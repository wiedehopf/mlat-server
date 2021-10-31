# -*- mode: python; indent-tabs-mode: nil -*-

# Part of mlat-server: a Mode S multilateration server
# Copyright (C) 2015  Oliver Jowett <oliver@mutability.co.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Maintains state for all aircraft known to some client.
Works out the set of aircraft we want the clients to send traffic for.
"""

import random
import asyncio
import time
from mlat import profile
from mlat.server import kalman, config

import logging
glogger = logging.getLogger("tracker")


class TrackedAircraft(object):
    """A single tracked aircraft."""

    def __init__(self, icao, allow_mlat):
        # ICAO address of this aircraft
        self.icao = icao

        # Allow mlat of this aircraft?
        self.allow_mlat = allow_mlat

        # set of receivers that can see this aircraft.
        # invariant: r.tracking.contains(a) iff a.tracking.contains(r)
        self.tracking = set()

        # set of receivers who use this aircraft for synchronization.
        # this aircraft is interesting if this set is non-empty.
        # invariant: r.sync_interest.contains(a) iff a.sync_interest.contains(r)
        self.sync_interest = set()

        # set of receivers who have seen ADS-B from this aircraft
        self.adsb_seen = set()

        # timestamp of when the last sync point was created using this aircraft
        # set this to 3 min in the past, hacky
        self.last_syncpoint_time = time.monotonic() - 180

        # set of receivers who want to use this aircraft for multilateration.
        # this aircraft is interesting if this set is non-empty.
        # invariant: r.mlat_interest.contains(a) iff a.mlat_interest.contains(r)
        self.mlat_interest = set()

        # set of receivers that have contributed to at least one multilateration
        # result. This is used to decide who to forward results to.
        self.successful_mlat = set()

        # number of mlat message resolves attempted
        self.mlat_message_count = 0
        # number of mlat messages that produced valid least-squares results
        self.mlat_result_count = 0
        # number of mlat messages that produced valid kalman state updates
        self.mlat_kalman_count = 0

        # last reported altitude (for multilaterated aircraft)
        self.altitude = None
        # time of last altitude (time.monotonic())
        self.last_altitude_time = None
        # altitude time tuples
        self.alt_history = []
        # dervided vertical rate
        self.vrate = None
        self.vrate_time = None

        # last multilateration, time (monotonic)
        self.last_result_time = None
        # last multilateration, ECEF position
        self.last_result_position = None
        # last multilateration, variance
        self.last_result_var = None
        # last multilateration, distinct receivers
        self.last_result_distinct = None
        # kalman filter state
        self.kalman = kalman.KalmanStateCA(self.icao)

        self.last_crappy_output = 0

        self.callsign = None
        self.squawk = None

    @property
    def interesting(self):
        """Is this aircraft interesting, i.e. are we asking any station to transmit data for it?"""
        return bool(self.sync_interest or (self.allow_mlat and self.mlat_interest))

    def __lt__(self, other):
        return self.icao < other.icao


class Tracker(object):
    """Tracks which receivers can see which aircraft, and asks receivers to
    forward traffic accordingly."""

    def __init__(self, coordinator, partition):
        self.aircraft = {}
        self.partition_id = partition[0] - 1
        self.partition_count = partition[1]
        self.coordinator = coordinator

    def in_local_partition(self, icao):
        if self.partition_count == 1:
            return True

        # mix the address a bit
        h = icao
        h = (((h >> 16) ^ h) * 0x45d9f3b) & 0xFFFFFFFF
        h = (((h >> 16) ^ h) * 0x45d9f3b) & 0xFFFFFFFF
        h = ((h >> 16) ^ h)
        return bool((h % self.partition_count) == self.partition_id)

    def add(self, receiver, icao_set):
        for icao in icao_set:
            ac = self.aircraft.get(icao)
            if ac is None:
                ac = self.aircraft[icao] = TrackedAircraft(icao, self.in_local_partition(icao))

            ac.tracking.add(receiver)
            receiver.tracking.add(ac)

    def remove(self, receiver, icao_set):
        for icao in icao_set:
            ac = self.aircraft.get(icao)
            if not ac:
                continue

            ac.tracking.discard(receiver)
            ac.successful_mlat.discard(receiver)
            receiver.tracking.discard(ac)
            if not ac.tracking:
                del self.aircraft[icao]

    def remove_all(self, receiver):
        for ac in receiver.tracking:
            ac.tracking.discard(receiver)
            ac.successful_mlat.discard(receiver)
            ac.sync_interest.discard(receiver)
            ac.adsb_seen.discard(receiver)
            ac.mlat_interest.discard(receiver)
            if not ac.tracking:
                del self.aircraft[ac.icao]

        receiver.tracking.clear()
        receiver.adsb_seen.clear()
        receiver.sync_interest.clear()
        receiver.mlat_interest.clear()

    @profile.trackcpu
    def update_interest(self, receiver):
        """Update the interest sets of one receiver based on the
        latest tracking and rate report data."""

        new_adsb = set()
        now = time.monotonic()


        new_mlat = {ac for ac in receiver.tracking if ac.allow_mlat and len(ac.adsb_seen) < 3}
        if receiver.last_rate_report is None:
            # Legacy client, no rate report, we cannot be very selective.
            new_sync = {ac for ac in receiver.tracking}
            if len(new_sync) > config.MAX_SYNC_AC:
                new_sync = set(random.sample(new_sync, k=config.MAX_SYNC_AC))

            receiver.update_interest_sets(new_sync, new_mlat, new_adsb)
            asyncio.get_event_loop().call_soon(receiver.refresh_traffic_requests)
            return


        # Work out the aircraft that are transmitting ADS-B that this
        # receiver wants to use for synchronization.
        ac_to_ratepair_map = {}
        ratepair_list = []
        rate_report_set = set()
        for icao, rate in receiver.last_rate_report.items():
            ac = self.aircraft.get(icao)
            if not ac:
                self.coordinator.receiver_tracking_add(receiver, {int(icao, 16)})

            rate_report_set.add(ac)
            new_adsb.add(ac)

            altFactor = None
            if ac.altitude is not None and ac.altitude > 0:
                altFactor = (1 + (ac.altitude / 20000)**1.5)
                #if receiver.user.startswith(config.DEBUG_FOCUS):
                #    glogger.warn('altFactor:' + str(altFactor) + ' alt: ' + str(ac.altitude))
            ac_to_ratepair_map[ac] = l = []  # list of (rateproduct, receiver, ac) tuples for this aircraft
            for r1 in ac.tracking:
                if receiver is r1:
                    continue

                if r1.last_rate_report is None:
                    # Receiver that does not produce rate reports, just take a guess.
                    rate1 = 0.8
                else:
                    rate1 = r1.last_rate_report.get(icao, 0.0)

                rp = rate * rate1 / 2.25
                # favor higher flying aircraft
                if altFactor is not None:
                    rp = rp * altFactor
                if rp < 0.01:
                    continue

                ratepair = (rp, r1, ac, rate)
                l.append(ratepair)
                ratepair_list.append(ratepair)

        ratepair_list.sort(reverse=True)

        ntotal = {}
        new_sync = set()
        total_rate = 0

        # select SYNC aircraft round1
        for rp, r1, ac, rate in ratepair_list:
            if ac in new_sync:
                continue  # already added

            if total_rate > config.MAX_SYNC_RATE:
                break

            if ntotal.get(r1, 0.0) < 0.3:
                # use this aircraft for sync
                new_sync.add(ac)
                total_rate += rate
                #if receiver.user.startswith(config.DEBUG_FOCUS):
                #    glogger.warn('1st round: ' + str(rate))
                # update rate-product totals for all receivers that see this aircraft
                for rp2, r2, ac2, rate in ac_to_ratepair_map[ac]:
                    ntotal[r2] = ntotal.get(r2, 0.0) + rp2

        # select SYNC aircraft round2 < 2.0 instead of < 1.0 ntotal
        for rp, r1, ac, rate in ratepair_list:
            if ac in new_sync:
                continue  # already added

            if total_rate > config.MAX_SYNC_RATE:
                break

            if ntotal.get(r1, 0.0) < 3.5:
                # use this aircraft for sync
                new_sync.add(ac)
                total_rate += rate
                #if receiver.user.startswith(config.DEBUG_FOCUS):
                #    glogger.warn('2nd round: ' + str(rate))
                # update rate-product totals for all receivers that see this aircraft
                for rp2, r2, ac2, rate in ac_to_ratepair_map[ac]:
                    ntotal[r2] = ntotal.get(r2, 0.0) + rp2

        if now - receiver.connectedSince < 45 and len(new_sync) < int(config.MAX_SYNC_AC / 2):
            acAvailable = ac_to_ratepair_map.keys()
            new_sync |= set(random.sample(acAvailable, k=min(len(acAvailable), int(config.MAX_SYNC_AC / 2))))

        addSome = int(config.MAX_SYNC_AC / 4) - len(new_sync)
        if addSome > 0:
            acAvailable = set(ac_to_ratepair_map.keys()).difference(new_sync)
            new_sync |= set(random.sample(acAvailable, k=min(len(acAvailable), addSome)))

        #if receiver.user.startswith(config.DEBUG_FOCUS):
        #    glogger.warn('new_sync:' + str([format(a.icao, '06x') for a in new_sync]))

        receiver.update_interest_sets(new_sync, new_mlat, new_adsb)
        asyncio.get_event_loop().call_soon(receiver.refresh_traffic_requests)
