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


import asyncio
import functools
import time
import logging
import math

from cpython cimport array
import array

import modes_cython.message

from mlat import geodesy, constants, profile
from mlat import clocksync, config
from libc.math cimport sqrt

__all__ = ('SyncPoint', 'ClockTracker')


cdef class SyncPoint(object):
    """A potential clock synchronization point.
    Clock synchronization points are a pair of DF17 messages,
    and associated timing info from all receivers that see
    that pair.
    """
    cdef public int address
    cdef public array.array posA
    cdef public array.array posB
    cdef public double interval
    cdef public list receivers

    def __init__(self, message_details, interval):
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
        address, posA, posB = message_details

        self.address = address
        self.posA = posA
        self.posB = posB
        self.interval = interval
        self.receivers = []  # a list of (receiver, timestampA, timestampB) values

cdef double ecef_distance(p0, p1):
    """Returns the straight-line distance in metres between two ECEF points."""
    return sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)

cdef get_limit(int cat):
    if cat == 0:
        return 30
    if cat == 1:
        return 16
    if cat == 2:
        return 6
    if cat == 3:
        return 6

cdef _add_to_existing_syncpoint(clock_pairs, syncpoint, r0, double t0A, double t0B):
    # add a new receiver and timestamps to an existing syncpoint

    # new state for the syncpoint: receiver, timestamp A, timestamp B,
    # and a flag indicating if this receiver actually managed to sync
    # with another receiver using this syncpoint (used for stats)

    cdef double receiverDistA = ecef_distance(syncpoint.posA, r0.position)
    cdef double receiverDistB = ecef_distance(syncpoint.posB, r0.position)

    cdef double now = time.time()

    # add receiver distance check here
    if receiverDistA > config.MAX_RANGE or receiverDistB > config.MAX_RANGE:
        r0.sync_range_exceeded = now
        return

    # propagation delays, in clock units
    #cdef double delayFactor = r0.clock.freq / constants.Cair
    cdef double delayFactor = r0.clock.delayFactor
    delay0A = receiverDistA * delayFactor
    delay0B = receiverDistB * delayFactor

    cdef double td0A = t0A - delay0A
    cdef double td0B = t0B - delay0B

    # compute interval, adjusted for transmitter motion
    cdef double i0 = td0B - td0A

    r0l = (r0, td0B, i0)

    # try to sync the new receiver with all receivers that previously
    # saw the same pair
    cdef double td1B, i1
    cdef int cat
    cdef double p0, p1, limit
    cdef dict distances = r0.distance
    for r1l in syncpoint.receivers:
        r1, td1B, i1 = r1l

        # order the clockpair so that the receiver that sorts lower is the base clock

        if r0 < r1:
            k = (r0, r1)
        else:
            k = (r1, r0)

        pairing = clock_pairs.get(k)

        if pairing is None:
            if r1.dead or r0 is r1:
                # receiver went away before we started resolving this
                # odd, but could happen
                continue

            cat = (int) (distances[r1.uid] / 50e3)
            if cat > 3:
                cat = 3

            limit = 0.7 * get_limit(cat)
            p0 = r0.sync_peers[cat]
            p1 = r1.sync_peers[cat]

            if p0 > limit and p1 > limit:
                #if r0.user.startswith(config.DEBUG_FOCUS) or r1.user.startswith(config.DEBUG_FOCUS):
                #    logging.warning("rejected new sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                continue

            clock_pairs[k] = pairing = clocksync.ClockPairing(r0, r1, cat)
            r0.sync_peers[cat] += 1
            r1.sync_peers[cat] += 1

        else:
            if now - pairing.updated < 0.5:
                continue
            cat = pairing.cat

            limit = 1.2 * get_limit(cat)
            p0 = r0.sync_peers[cat]
            p1 = r1.sync_peers[cat]


            if p0 > limit and p1 > limit:
                #if r0.user.startswith(config.DEBUG_FOCUS) or r1.user.startswith(config.DEBUG_FOCUS):
                #    logging.warning("rejected existing sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                r0.sync_peers[cat] -= 1
                r1.sync_peers[cat] -= 1
                del clock_pairs[k]
                continue

        if r0 < r1:
            pairing.update(syncpoint.address, td0B, td1B, i0, i1, now)
        else:
            pairing.update(syncpoint.address, td1B, td0B, i1, i0, now)

    # update syncpoint with the new receiver and we're done
    syncpoint.receivers.append(r0l)

class ClockTracker(object):
    """Maintains clock pairings between receivers, and matches up incoming sync messages
    from receivers to update the parameters of the pairings."""

    def __init__(self, coordinator, loop):
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
        self.loop = loop

        # schedule periodic cleanup
        self.loop.call_later(1.0, self._cleanup)

    def _cleanup(self):
        """Called periodically to clean up clock pairings that have expired and update pairing.valid"""

        self.loop.call_later(10.0, self._cleanup)

        now = time.time()
        prune = set()
        for k, pairing in self.clock_pairs.items():
            pairing.check_valid(now)
            if now - pairing.updated > 45:
                prune.add((k, pairing))
            if not pairing.valid and now - pairing.updated > 30:
                prune.add((k, pairing))

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
                      double even_time, double odd_time,
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

        # compute key and interval
        cdef double tA, tB
        if even_time < odd_time:
            tA = even_time
            tB = odd_time
            key = (even_message, odd_message)
        else:
            tA = odd_time
            tB = even_time
            key = (odd_message, even_message)

        # do we have a suitable existing match?
        message_details, syncpointlist = self.sync_points.get(key, (None, None))

        # check if sync point is invalid
        if syncpointlist == False:
            return

        cdef double freq = receiver.clock.freq
        cdef double interval = (tB - tA) / freq

        # Messages must be within 5 seconds of each other.
        if interval > 5.0:
            return

        if isinstance(syncpointlist, list):
            for candidate in syncpointlist:
                if abs(candidate.interval - interval) < 1e-3:
                    # interval matches within 1ms, close enough.
                    _add_to_existing_syncpoint(self.clock_pairs, candidate, receiver, tA, tB)
                    return
            # no matching syncpoint in syncpointlist, add one
            syncpoint = SyncPoint(message_details, interval)

            syncpointlist.append(syncpoint)

            _add_to_existing_syncpoint(self.clock_pairs, syncpoint, receiver, tA, tB)

            return

        now = time.time()

        # this receiver isn't allowed to create new sync points due to exceeding the range
        if now - receiver.sync_range_exceeded < 15.0:
            return

        message_details = (None, None, None)
        if syncpointlist is None:
            # No existing match. Create an invalid sync point, if it pans out, replace it.
            self.sync_points[key] = (message_details, False)
            # schedule cleanup of the syncpoint after 3 seconds -
            # we should have seen all copies of those messages by then.
            #self.loop.call_later(3.0, functools.partial(self._cleanup_syncpointlist,key=key))
            # change to centralized cleanup every 15 seconds, maybe save some CPU

        # Validate the messages and maybe create a real sync point list

        # basic validity
        even_message = modes_cython.message.decode(even_message)
        odd_message = modes_cython.message.decode(odd_message)

        ac = self.coordinator.tracker.aircraft.get(even_message.address)
        if ac:
            ac.seen = now
        if (ac
                and even_message.estype == modes_cython.message.ESType.surface_position
                and odd_message.estype == modes_cython.message.ESType.surface_position
                ):
            ac.last_adsb_time = now

        if ((not even_message or
             even_message.DF != 17 or
             not even_message.crc_ok or
             even_message.estype != modes_cython.message.ESType.airborne_position or
             even_message.F)):
            return

        if ((not odd_message or
             odd_message.DF != 17 or
             not odd_message.crc_ok or
             odd_message.estype != modes_cython.message.ESType.airborne_position or
             not odd_message.F)):
            return

        if even_message.address != odd_message.address:
            return

        # find global positions
        try:
            even_lat, even_lon, odd_lat, odd_lon = modes_cython.message.decode_cpr(even_message.LAT,
                                                                    even_message.LON,
                                                                    odd_message.LAT,
                                                                    odd_message.LON)
        except ValueError:
            # CPR failed
            return

        # quality checks
        if even_message.altitude is None or odd_message.altitude is None:
            return

        if abs(even_message.altitude - odd_message.altitude) > 5000:
            return

        # sort out some bad transponders / strange positions
        if even_lat > 85 or even_lat < -85:
            return
        if even_lat == 0 or even_lon == 0:
            return

        # convert to ECEF, do range checks
        even_ecef = geodesy.llh2ecef((even_lat,
                                      even_lon,
                                      even_message.altitude * constants.FTOM))


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

        #do some extra bookkeeping now we know this message is legit
        if ac:
            ac.last_adsb_time = now
            ac.last_altitude_time = now
            ac.altitude = even_message.altitude

        # more quality checking
        if even_message.nuc < 6 or odd_message.nuc < 6:
            return

        # valid message, set the message details for use by the SyncPoint
        if even_time < odd_time:
            message_details = (even_message.address, even_ecef, odd_ecef)
        else:
            message_details = (even_message.address, odd_ecef, even_ecef)

        if geodesy.ecef_distance(even_ecef, receiver.position) > config.MAX_RANGE:
            # suppress this spam, can't help if ppl give a wrong location
            # logging.info("{a:06X}: receiver range check (even) failed".format(a=even_message.address))
            receiver.sync_range_exceeded = now
            self.sync_points[key] = (message_details, 'dummy')
            return

        # valid. Create a new Sync point, add to it and creat the sync point list

        syncpoint = SyncPoint(message_details, interval)

        _add_to_existing_syncpoint(self.clock_pairs, syncpoint, receiver, tA, tB)

        self.sync_points[key] = (message_details, [ syncpoint ])



    @profile.trackcpu
    def _cleanup_syncpointlist(self, key):
        """Expire a syncpoint list. This happens ~3 seconds after the first copy
        of a message pair is received.

        key: the key of the syncpoint
        """

        try:
            del self.sync_points[key]
        except KeyError:
            # something is not right, reset the dict
            logging.warn("len(sync_points): {a}".format(a=len(self.sync_points)))
            self.sync_points.clear()

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
