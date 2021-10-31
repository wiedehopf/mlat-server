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
Manages the clock synchronization pairs between all receivers based on
DF17 position messages received by more than one receiver.
"""

__all__ = ('SyncPoint', 'ClockTracker')

import asyncio
import functools
import time
import logging
import math

import modes.message

from mlat import geodesy, constants, profile
from mlat.server import clocksync, config


class SyncPoint(object):
    """A potential clock synchronization point.
    Clock synchronization points are a pair of DF17 messages,
    and associated timing info from all receivers that see
    that pair.
    """

    def __init__(self, address, posA, posB, interval):
        """Construct a new sync point.

        address: the ICAO address of the sync aircraft
        posA: the ECEF position of the earlier message
        posB: the ECEF position of the later message
        interval: the nominal interval (in seconds)
          between the two messages; this is as measured by
          the first receiver to report the pair, and is used
          to distinguish cases where the same message is
          transmitted more than once.
        """

        self.address = address
        self.posA = posA
        self.posB = posB
        self.interval = interval
        self.receivers = []  # a list of (receiver, timestampA, timestampB) values


class ClockTracker(object):
    """Maintains clock pairings between receivers, and matches up incoming sync messages
    from receivers to update the parameters of the pairings."""

    def __init__(self, coordinator):
        # map of (sync key) -> list of sync points
        #
        # sync key is a pair of bytearrays: (msgA, msgB)
        # where msgA and msgB are the contents of the
        # earlier and later message of the pair respectively.
        self.sync_points = {}

        # map of (pair key) -> pairing
        #
        # pair key is (receiver 0, receiver 1) where receiver 0
        # is always less than receiver 1.
        self.clock_pairs = {}

        self.coordinator = coordinator

        # schedule periodic cleanup
        asyncio.get_event_loop().call_later(1.0, self._cleanup)

    def _cleanup(self):
        """Called periodically to clean up clock pairings that have expired and update pairing.valid"""

        asyncio.get_event_loop().call_later(5.0, self._cleanup)

        now = time.monotonic()
        prune = set()
        for k, pairing in self.clock_pairs.items():
            if pairing.expiry <= now or k[0].bad_syncs > 2 or k[1].bad_syncs > 2:
                prune.add((k, pairing))
            else:
                pairing.valid = pairing.check_valid(now)

        for k, pairing in prune:
            k[0].sync_peers[pairing.cat] -= 1
            k[1].sync_peers[pairing.cat] -= 1
            del self.clock_pairs[k]

    @profile.trackcpu
    def receiver_clock_reset(self, receiver):
        """
        Called by the coordinator when we should drop our clock sync
        state for a given receiver. This happens on input disconnect/
        reconnect.

        (This is actually the same work as receiver_disconnect for the moment)
        """
        for k, pairing in list(self.clock_pairs.items()):
            if k[0] is receiver or k[1] is receiver:
                k[0].sync_peers[pairing.cat] -= 1
                k[1].sync_peers[pairing.cat] -= 1
                del self.clock_pairs[k]

    @profile.trackcpu
    def receiver_disconnect(self, receiver):
        """
        Called by the coordinator when a receiver disconnects.

        Clears up any clock pairing involving the receiver immediately,
        as it's very likely that any existing sync data will be invalid
        if/when the receiver later reconnects.

        Sync points involving the receiver are not cleaned up immediately.
        It's assumed that the disconnected receiver has the "dead" flag
        set; this flag is tested before sync happens.
        """

        # Clean up clock_pairs immediately.
        # Any membership in a pending sync point is noticed when we try to sync more receivers with it.
        for k, pairing in list(self.clock_pairs.items()):
            if k[0] is receiver or k[1] is receiver:
                k[0].sync_peers[pairing.cat] -= 1
                k[1].sync_peers[pairing.cat] -= 1
                del self.clock_pairs[k]

    @profile.trackcpu
    def receiver_sync(self, receiver,
                      even_time, odd_time,
                      even_message, odd_message):
        """
        Called by the coordinator to handle a sync message from a receiver.

        Looks for a suitable existing sync point and, if there is one, does
        synchronization between this receiver and the existing receivers
        associated with the sync point.

        Otherwise, validates the message pair and, if it is suitable, creates a
        new sync point for it.

        receiver: the receiver reporting the sync message
        even_message: a DF17 airborne position message with F=0
        odd_message: a DF17 airborne position message with F=1
        even_time: the time of arrival of even_message, as seen by receiver.clock
        odd_time: the time of arrival of odd_message, as seen by receiver.clock
        """

        # Do sanity checks.

        if receiver.bad_syncs > 2:
            return

        # compute key and interval
        if even_time < odd_time:
            tA = even_time
            tB = odd_time
            key = (even_message, odd_message)
        else:
            tA = odd_time
            tB = even_time
            key = (odd_message, even_message)

        interval = (tB - tA) / receiver.clock.freq

        # Messages must be within 5 seconds of each other.
        if interval > 5.0:
            return

        # do we have a suitable existing match?
        syncpointlist = self.sync_points.get(key)
        if syncpointlist:
            for candidate in syncpointlist:
                if abs(candidate.interval - interval) < 1e-3:
                    # interval matches within 1ms, close enough.
                    self._add_to_existing_syncpoint(candidate, receiver, tA, tB)
                    return

        # No existing match. Validate the messages and maybe create a new sync point

        if receiver.sync_range_exceeded:
            return

        # basic validity
        even_message = modes.message.decode(even_message)
        if ((not even_message or
             even_message.DF != 17 or
             not even_message.crc_ok or
             even_message.estype != modes.message.ESType.airborne_position or
             even_message.F)):
            return

        odd_message = modes.message.decode(odd_message)
        if ((not odd_message or
             odd_message.DF != 17 or
             not odd_message.crc_ok or
             odd_message.estype != modes.message.ESType.airborne_position or
             not odd_message.F)):
            return

        if even_message.address != odd_message.address:
            return

        # quality checks
        if even_message.nuc < 6 or even_message.altitude is None:
            return

        if odd_message.nuc < 6 or odd_message.altitude is None:
            return

        if abs(even_message.altitude - odd_message.altitude) > 5000:
            return

        # find global positions
        try:
            even_lat, even_lon, odd_lat, odd_lon = modes.cpr.decode(even_message.LAT,
                                                                    even_message.LON,
                                                                    odd_message.LAT,
                                                                    odd_message.LON)
        except ValueError:
            # CPR failed
            return

        # convert to ECEF, do range checks
        even_ecef = geodesy.llh2ecef((even_lat,
                                      even_lon,
                                      even_message.altitude * constants.FTOM))
        if geodesy.ecef_distance(even_ecef, receiver.position) > config.MAX_RANGE:
            # suppress this spam, can't help if ppl give a wrong location
            # logging.info("{a:06X}: receiver range check (even) failed".format(a=even_message.address))
            receiver.sync_range_exceeded = 1
            return

        odd_ecef = geodesy.llh2ecef((odd_lat,
                                     odd_lon,
                                     odd_message.altitude * constants.FTOM))

        # checking range for the even position and intermessage distance is sufficient
        #if geodesy.ecef_distance(odd_ecef, receiver.position) > config.MAX_RANGE:
            #logging.info("{a:06X}: receiver range check (odd) failed".format(a=odd_message.address))
            #return

        if geodesy.ecef_distance(even_ecef, odd_ecef) > config.MAX_INTERMESSAGE_RANGE:
            #logging.info("{a:06X}: intermessage range check failed".format(a=even_message.address))
            return

        #valid, do some extra bookkeeping before sync stuff

        ac = self.coordinator.tracker.aircraft.get(even_message.address)
        if ac:
            now = time.monotonic()
            ac.last_syncpoint_time = now
            ac.last_altitude_time = now
            ac.altitude = even_message.altitude


        # valid. Create a new sync point.
        if even_time < odd_time:
            syncpoint = SyncPoint(even_message.address, even_ecef, odd_ecef, interval)
        else:
            syncpoint = SyncPoint(even_message.address, odd_ecef, even_ecef, interval)

        self._add_to_existing_syncpoint(syncpoint, receiver, tA, tB)

        if not syncpointlist:
            syncpointlist = self.sync_points[key] = []
        syncpointlist.append(syncpoint)

        # schedule cleanup of the syncpoint after 2 seconds -
        # we should have seen all copies of those messages by
        # then.
        asyncio.get_event_loop().call_later(
            2.0,
            functools.partial(self._cleanup_syncpoint,
                              key=key,
                              syncpoint=syncpoint))

    def _add_to_existing_syncpoint(self, syncpoint, r0, t0A, t0B):
        # add a new receiver and timestamps to an existing syncpoint

        # new state for the syncpoint: receiver, timestamp A, timestamp B,
        # and a flag indicating if this receiver actually managed to sync
        # with another receiver using this syncpoint (used for stats)

        receiverDistA = geodesy.ecef_distance(syncpoint.posA, r0.position)
        receiverDistB = geodesy.ecef_distance(syncpoint.posB, r0.position)

        # add receiver distance check here
        if receiverDistA > config.MAX_RANGE or receiverDistB > config.MAX_RANGE:
            r0.sync_range_exceeded = 1
            return

        r0.sync_range_exceeded = 0

        # propagation delays, in clock units
        delay0A = receiverDistA * r0.clock.freq / constants.Cair
        delay0B = receiverDistB * r0.clock.freq / constants.Cair

        td0A = t0A - delay0A
        td0B = t0B - delay0B

        # compute interval, adjusted for transmitter motion
        i0 = td0B - td0A

        r0l = [r0, td0B, i0, False]

        now = time.monotonic()

        # try to sync the new receiver with all receivers that previously
        # saw the same pair
        for r1l in syncpoint.receivers:
            r1, td1B, i1, r1sync = r1l

            if r1.dead:
                # receiver went away before we started resolving this
                continue

            if r0 is r1:
                # odd, but could happen
                continue

            # order the clockpair so that the receiver that sorts lower is the base clock

            if r0 < r1:
                k = (r0, r1)
            else:
                k = (r1, r0)

            pairing = self.clock_pairs.get(k)
            if pairing is None:
                receiver_distance = r0.distance[r1.uid]
                cat = math.floor(receiver_distance / config.DISTANCE_CATEGORY_STEP)
                if cat > 4:
                    cat = 4

                p0 = r0.sync_peers[cat]
                p1 = r1.sync_peers[cat]
                limit = config.MAX_PEERS[cat]

                if p0 > limit or p1 > limit:
                    if p0 > 0.4 * limit and p1 > 0.4 * limit:
                        #if r0.user.startswith("euerdorf") or r1.user.startswith("euerdorf"):
                        if r0.user.startswith("Kirby") or r1.user.startswith("Kirby"):
                            logging.warning("rejected new sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                        continue

                self.clock_pairs[k] = pairing = clocksync.ClockPairing(r0, r1, cat)

                r0.sync_peers[pairing.cat] += 1
                r1.sync_peers[pairing.cat] += 1
            else:
                if pairing.n > 10 and now < pairing.updated + 0.8:
                    continue

                cat = pairing.cat
                p0 = r0.sync_peers[cat]
                p1 = r1.sync_peers[cat]
                limit = config.MAX_PEERS[cat] * 1.2

                if p0 > limit or p1 > limit:
                    if p0 > 0.7 * limit and p1 > 0.7 * limit:
                        #if r0.user.startswith("euerdorf") or r1.user.startswith("euerdorf"):
                        if r0.user.startswith("Kirby") or r1.user.startswith("Kirby"):
                            logging.warning("rejected existing sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                        r0.sync_peers[pairing.cat] -= 1
                        r1.sync_peers[pairing.cat] -= 1
                        del self.clock_pairs[k]
                        continue

            if r0 < r1:
                if not pairing.update(syncpoint.address, td0B, td1B, i0, i1, now):
                    continue
            else:
                if not pairing.update(syncpoint.address, td1B, td0B, i1, i0, now):
                    continue

            # sync worked, note it for stats
            r0l[3] = r1l[3] = True

        # update syncpoint with the new receiver and we're done
        syncpoint.receivers.append(r0l)

    @profile.trackcpu
    def _cleanup_syncpoint(self, key, syncpoint):
        """Expire a syncpoint. This happens ~2 seconds after the first copy
        of a message pair is received.

        key: the key of the syncpoint
        syncpoint: the syncpoint itself
        """

        # remove syncpoint from self.sync_points, clean up empty entries
        l = self.sync_points[key]
        l.remove(syncpoint)
        if not l:
            del self.sync_points[key]

        # stats update
        for r, _, _, synced in syncpoint.receivers:
            if synced:
                r.sync_count += 1

    def dump_receiver_state(self):
        state = {}
        for (r0, r1), pairing in self.clock_pairs.items():
            if pairing.n < 2:
                continue

            state.setdefault(r0.user, {})[r1.user] = [pairing.n,
                              round(pairing.error * 1e6, 1),
                              round(pairing.drift * 1e6),
                              round(r1.bad_syncs, 2),
                              pairing.jumped]
                    #removed: #pairing.ts_peer[-1] - pairing.ts_base[-1]]
            state.setdefault(r1.user, {})[r0.user] = [pairing.n,
                              round(pairing.error * 1e6, 1),
                              round(pairing.i_drift * 1e6),
                              round(r0.bad_syncs, 2),
                              pairing.jumped]
                    #removed: #pairing.ts_base[-1] - pairing.ts_peer[-1]]
            # reset jumped indicator
            pairing.jumped = 0
        return state
