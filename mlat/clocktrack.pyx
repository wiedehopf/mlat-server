#!python
#cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, nonecheck=False

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

Maintains clock synchronization between individual pairs of receivers.
"""

import asyncio
import functools
import time
import logging
import math
import bisect

from cpython cimport array
import array

import modes_cython.message

from mlat import geodesy, constants, profile, config
# cython stuff:
from libc.math cimport sqrt
from libc.string cimport memmove

import pygraph.classes.graph
import pygraph.algorithms.minmax
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np


__all__ = ('SyncPoint', 'ClockTracker', 'Clock', 'ClockPairing', 'make_clock')

glogger = logging.getLogger("clocksync")

cdef double MAX_RANGE = 500e3

cdef class SyncShard:
    """part of a sync point, data for that receiver in relation to the sync point
    """
    cdef receiver
    cdef double td
    cdef double i

cdef class SyncPoint:
    """A potential clock synchronization point.
    Clock synchronization points are a pair of DF17 messages,
    and associated timing info from all receivers that see
    that pair.
    """
    cdef int address
    cdef array.array posA
    cdef array.array posB
    cdef double interval
    cdef list receivers
    cdef ac

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
        self.address, self.posA, self.posB, self.ac, do_mlat = message_details

        self.interval = interval
        self.receivers = []  # a list of (receiver, timestampA, timestampB) values

cdef double ecef_distance(p0, p1):
    """Returns the straight-line distance in metres between two ECEF points."""
    return sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)

cdef get_limit(int cat):
    if cat == 0:
        return 32
    if cat == 1:
        return 32
    if cat == 2:
        return 16
    if cat == 3:
        return 16

cdef _add_to_existing_syncpoint(dict clock_pairs, SyncPoint syncpoint, r0, double t0A, double t0B):
    # add a new receiver and timestamps to an existing syncpoint

    # new state for the syncpoint: receiver, timestamp A, timestamp B,
    # and a flag indicating if this receiver actually managed to sync
    # with another receiver using this syncpoint (used for stats)

    cdef double receiverDistA = ecef_distance(syncpoint.posA, r0.position)
    cdef double receiverDistB = ecef_distance(syncpoint.posB, r0.position)

    cdef double now = time.time()

    # receiver distance check, checking one range is sufficient
    if receiverDistA > MAX_RANGE:
        r0.sync_range_exceeded = now
        return

    # propagation delays, in clock units
    #cdef double delayFactor = r0.clock.freq / constants.Cair
    cdef double delayFactor = r0.clock.delayFactor
    cdef double delay0A = receiverDistA * delayFactor
    cdef double delay0B = receiverDistB * delayFactor

    cdef double td0A = t0A - delay0A
    cdef double td0B = t0B - delay0B

    # compute interval, adjusted for transmitter motion
    cdef double i0 = td0B - td0A

    cdef SyncShard r0l = SyncShard()
    cdef SyncShard r1l
    r0l.receiver = r0
    r0l.td = td0B
    r0l.i = i0

    # try to sync the new receiver with all receivers that previously
    # saw the same pair
    cdef double td1B, i1
    cdef int cat
    cdef double p0, p1, limit
    cdef dict distances = r0.distance
    cdef ClockPairing pairing

    for r1l in syncpoint.receivers:
        r1 = r1l.receiver
        td1B = r1l.td
        i1 = r1l.i

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
                #if r0.focus or r1.focus:
                #    logging.warning("rejected new sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                continue

            clock_pairs[k] = pairing = ClockPairing(r0, r1, cat)
            r0.sync_peers[cat] += 1
            r1.sync_peers[cat] += 1

        else:
            if pairing.n > 8 and now - pairing.update_attempted < 0.2:
                continue
            cat = pairing.cat

            limit = 1.2 * get_limit(cat)
            p0 = r0.sync_peers[cat]
            p1 = r1.sync_peers[cat]


            if p0 > limit and p1 > limit:
                #if r0.focus or r1.focus:
                #    logging.warning("rejected existing sync: %06x cat: %d p0: %d p1: %d limit: %d", syncpoint.address, cat, p0, p1, limit)
                r0.sync_peers[cat] -= 1
                r1.sync_peers[cat] -= 1
                del clock_pairs[k]
                continue

        if r0 < r1:
            pairing.update(syncpoint.address, td0B, td1B, i0, i1, now, syncpoint.ac)
        else:
            pairing.update(syncpoint.address, td1B, td0B, i1, i0, now, syncpoint.ac)

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
        cdef ClockPairing pairing
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

        Only reset the offsets, the drift shouldn't be affected
        """
        cdef ClockPairing pairing
        for k, pairing in list(self.clock_pairs.items()):
            if k[0] is receiver or k[1] is receiver:
                pairing.reset_offsets()

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

        cdef ClockPairing pairing
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

        cdef SyncPoint syncpoint

        self.coordinator.stats_sync_msgs += 1
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
        message_details, syncpointlist = self.sync_points.get(key, (0, 0))

        if message_details == -1:
            return

        now = time.time()

        if message_details and message_details[4]:
            # do MLAT on this sync message if do_mlat was set in the message details
            self.coordinator.receiver_mlat(receiver, even_time, even_message, now)

        # check if sync point is invalid
        if syncpointlist == -1:
            return

        cdef double freq = receiver.clock.freq
        cdef double interval = (tB - tA) / freq

        # Messages must be within 5 seconds of each other.
        if interval > 5.0:
            return

        if syncpointlist:
            for syncpoint in syncpointlist:
                if abs(syncpoint.interval - interval) < 0.75e-3:
                    # interval matches within 0.75ms, close enough.
                    _add_to_existing_syncpoint(self.clock_pairs, syncpoint, receiver, tA, tB)
                    return
            # no matching syncpoint in syncpointlist, add one
            syncpoint = SyncPoint(message_details, interval)

            syncpointlist.append(syncpoint)

            _add_to_existing_syncpoint(self.clock_pairs, syncpoint, receiver, tA, tB)

            return

        # this receiver isn't allowed to create new sync points due to exceeding the range
        if now - receiver.sync_range_exceeded < 15.0:
            return

        if syncpointlist:
            # No existing match. Create an invalid sync point, if it pans out, replace it.
            self.sync_points[key] = (-1, -1)
            # schedule cleanup of the syncpoint after 3 seconds -
            # we should have seen all copies of those messages by then.
            #self.loop.call_later(3.0, functools.partial(self._cleanup_syncpointlist,key=key))
            # change to centralized cleanup every 15 seconds, maybe save some CPU

        # Validate the messages and maybe create a real sync point list

        # basic validity
        even_message = modes_cython.message.decode(even_message)
        odd_message = modes_cython.message.decode(odd_message)
        if (not even_message or not odd_message):
            return

        ac = self.coordinator.tracker.aircraft.get(even_message.address)
        if not ac:
            return

        ac.seen = now

        if ac and (ac.do_mlat or ac.force_mlat) and len(ac.tracking) >= 3:
            do_mlat = True
            self.coordinator.receiver_mlat(receiver, even_time, even_message, now)
        else:
            do_mlat = False

        # set dummy message details for potential use for MLATing bad ADS-B data
        message_details = (even_message.address, None, None, ac, do_mlat)
        self.sync_points[key] = (message_details, -1)

        if ((even_message.DF != 17 or
             not even_message.crc_ok or
             even_message.F)):
            return

        if ((odd_message.DF != 17 or
             not odd_message.crc_ok or
             not odd_message.F)):
            return

        if even_message.address != odd_message.address:
            return

        if (odd_message.estype != modes_cython.message.ESType.airborne_position or
            even_message.estype != modes_cython.message.ESType.airborne_position):
            if (even_message.estype == modes_cython.message.ESType.surface_position
                    and odd_message.estype == modes_cython.message.ESType.surface_position):
                ac.last_adsb_time = now

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

        if ecef_distance(even_ecef, odd_ecef) > config.MAX_INTERMESSAGE_RANGE:
            #logging.info("{a:06X}: intermessage range check failed".format(a=even_message.address))
            return

        #do some extra bookkeeping now we know this message is legit
        ac.last_altitude_time = now
        ac.altitude = even_message.altitude

        # more quality checking
        if even_message.nuc < 6 or odd_message.nuc < 6:
            return

        # valid message, set the message details for use by the SyncPoint
        if even_time < odd_time:
            message_details = (even_message.address, even_ecef, odd_ecef, ac, do_mlat)
        else:
            message_details = (even_message.address, odd_ecef, even_ecef, ac, do_mlat)

        # valid. Create a new Sync point, add to it and create the sync point list

        self.coordinator.stats_sync_points += 1

        syncpoint = SyncPoint(message_details, interval)

        self.sync_points[key] = (message_details, [ syncpoint ])

        if ecef_distance(even_ecef, receiver.position) > MAX_RANGE:
            # suppress this spam, can't help if ppl give a wrong location
            # logging.info("{a:06X}: receiver range check (even) failed".format(a=even_message.address))
            receiver.sync_range_exceeded = now
            return

        _add_to_existing_syncpoint(self.clock_pairs, syncpoint, receiver, tA, tB)

    def clear_all_sync_points(self):
        """Clear all syncpoint lists from sync_points.
            called every 15 seconds by coordinator
        """
        self.sync_points.clear()

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
            self.clear_all_sync_points()

    def dump_receiver_state(self):
        state = {}
        cdef double outlier_percent
        cdef ClockPairing pairing
        for (r0, r1), pairing in self.clock_pairs.items():
            if pairing.n < 2:
                continue

            if pairing.update_total < 4:
                outlier_percent = 50.0 * pairing.outlier_total / pairing.update_total
            else:
                outlier_percent = 100.0 * pairing.outlier_total / pairing.update_total


            pairing.outlier_total /= 2
            pairing.update_total /= 2


            state.setdefault(r0.user, {})[r1.user] = [pairing.n,
                              round(pairing.error * 1e6, 1),
                              round(pairing.drift * 1e6),
                              round(r1.bad_syncs, 2),
                              pairing.jumped,
                              round(outlier_percent, 1),
                              pairing.updated]
                    #removed: #pairing.ts_peer[-1] - pairing.ts_base[-1]]
            state.setdefault(r1.user, {})[r0.user] = [pairing.n,
                              round(pairing.error * 1e6, 1),
                              round(pairing.i_drift * 1e6),
                              round(r0.bad_syncs, 2),
                              pairing.jumped,
                              round(outlier_percent, 1),
                              pairing.updated]
                    #removed: #pairing.ts_base[-1] - pairing.ts_peer[-1]]
            # reset jumped indicator
            pairing.jumped = 0
        return state

cdef class Clock(object):
    """A particular clock. Stores characteristics of a clock,
    and acts as part of the key in the clock pairing map.
    """

    cdef readonly double freq
    cdef readonly double max_freq_error
    cdef readonly double jitter
    cdef readonly double delayFactor

    def __init__(self, freq, max_freq_error, jitter):
        """Create a new clock representation.

        freq: the clock frequency in Hz (float)
        max_freq_error: the maximum expected relative frequency error (i.e. 1e-6 is 1PPM) (float)
        jitter: the expected jitter of a typical reading, in seconds, standard deviation  (float)
        """
        self.freq = freq
        self.max_freq_error = max_freq_error
        self.jitter = jitter
        self.delayFactor = freq / constants.Cair


def make_clock(clock_type):
    """Return a new Clock instance for the given clock type."""

    if clock_type == 'radarcape_gps':
        return Clock(freq=1e9, max_freq_error=1e-6, jitter=15e-9)
    if clock_type == 'beast' or clock_type == 'radarcape_12mhz':
        return Clock(freq=12e6, max_freq_error=5e-6, jitter=83e-9)
    if clock_type == 'sbs':
        return Clock(freq=20e6, max_freq_error=100e-6, jitter=500e-9)
    if clock_type == 'dump1090' or clock_type == 'unknown':
        return Clock(freq=12e6, max_freq_error=100e-6, jitter=500e-9)
    raise NotImplementedError("{ct}".format(ct=clock_type))

cdef int cp_size = 32

cdef int drift_n_stable = 12

cdef class ClockPairing(object):
    """Describes the current relative characteristics of a pair of clocks."""

    cdef readonly bint valid
    cdef readonly double updated
    cdef readonly double variance
    cdef readonly int n
    cdef double update_attempted
    cdef base
    cdef peer
    cdef int cat

    cdef public double factor
    cdef public double i_factor
    cdef public double base_avg
    cdef public double peer_avg

    cdef base_clock
    cdef peer_clock
    cdef double base_freq
    cdef double peer_freq
    cdef double raw_drift
    cdef double drift
    cdef double i_drift

    cdef int drift_n
    cdef int drift_outliers
    cdef int outlier_reset_cooldown
    cdef double outlier_total
    cdef double update_total
    # needs to be cp_size big, can't use it here though
    cdef double ts_base[32]
    cdef double ts_peer[32]
    cdef double var[32]
    cdef double var_sum
    cdef int outliers
    cdef double cumulative_error
    cdef double error

    cdef int jumped

    cdef double relative_freq
    cdef double i_relative_freq
    cdef double drift_max
    cdef double drift_max_delta
    cdef double outlier_threshold
    cdef double update_last_sync


    def __init__(self, base, peer, cat):
        self.base = base
        self.peer = peer
        self.cat = cat
        self.base_clock = base.clock
        self.peer_clock = peer.clock
        self.base_freq = base.clock.freq
        self.peer_freq = peer.clock.freq

        self.relative_freq = peer.clock.freq / base.clock.freq
        self.i_relative_freq = base.clock.freq / peer.clock.freq
        self.drift_max = 0.75 * (base.clock.max_freq_error + peer.clock.max_freq_error)
        self.drift_max_delta = self.drift_max / 20.0
        # self.outlier_threshold = 4 * sqrt(peer.clock.jitter ** 2 + base.clock.jitter ** 2) # 4 sigma
        # this was about 2.5 us for rtl-sdr receivers
        self.outlier_threshold = 1.5 * 1e-6 # 1e-6 -> 1 us

        self.updated = 0
        self.update_attempted = 0

        self.raw_drift = 0
        self.drift = 0
        self.i_drift = 0
        self.drift_outliers = 0

        self.outliers = 0
        self.jumped = 0

        self.outlier_reset_cooldown = 5 # number of sync pair updates before this sync pair can be trusted

        self.valid = False
        self.n = 0
        self.var_sum = 0.0
        self.cumulative_error = 0.0
        self.error = -1e-6
        self.variance = -1e-6

        self.outlier_total = 0
        self.update_total = 1e-3

        self.update_last_sync = 0

    cdef bint check_valid(self, double now):
        if self.n < 2 or self.drift_n < 2 or self.outlier_total / self.update_total > 0.5:
            self.variance = -1e-6
            self.error = -1e-6
            self.valid = False
            return False

        """Variance of recent predictions of the sync point versus the actual sync point."""
        self.variance = self.var_sum / self.n
        """Standard error of recent predictions."""
        self.error = sqrt(self.variance)

        """True if this pairing is usable for clock syncronization."""
        self.valid = (self.outlier_reset_cooldown < 1
                and self.n > 4
                and self.drift_n > 4
                and self.variance < 16e-12
                and now - self.updated < 35.0)

        if self.valid and now > self.update_last_sync:
            self.update_last_sync = now + 5
            self.base.last_sync = now
            self.peer.last_sync = now

        return self.valid

    cdef update(self, address, double base_ts, double peer_ts, double base_interval, double peer_interval, double now, ac):
        """Update the relative drift and offset of this pairing given:

        address: the ICAO address of the sync aircraft, for logging purposes
        base_ts: the timestamp of a recent point in time measured by the base clock
        peer_ts: the timestamp of the same point in time measured by the peer clock
        base_interval: the duration of a recent interval measured by the base clock
        peer_interval: the duration of the same interval measured by the peer clock

        Returns True if the update was used, False if it was an outlier.
        """
        cdef double prediction = 0
        cdef double prediction_error = 0
        cdef bint outlier = False
        cdef bint do_reset = False
        cdef double outlier_threshold

        # clean old data
        if self.n > cp_size - 1 or base_ts - self.ts_base[0] > 50.0 * self.base_freq:
            self._prune_old_data(now)

        if not ac.sync_dont_use:
            self.update_total += 1
        self.update_attempted = now

        if self.n > 0 and not outlier:
            # ts_base and ts_peer define a function constructed by linearly
            # interpolating between each pair of values.
            #
            # This function must be monotonically increasing or one of our clocks
            # has effectively gone backwards. If this happens, give up and start
            # again.

            if peer_ts <= self.ts_peer[self.n - 1] or base_ts <= self.ts_base[self.n - 1]:
                if peer_ts < self.ts_peer[self.n - 1] and base_ts < self.ts_base[self.n - 1]:
                    return False
                if peer_ts == self.ts_peer[self.n - 1] or base_ts == self.ts_base[self.n - 1]:
                    return False

                # just in case, make this pair invalid for the moment
                # the next update will set it to valid again
                self.valid = False

                self.outliers += 10
                outlier = True
                self.outlier_total += 1

                if self.outliers <= 10:
                    # don't reset quite yet, maybe something strange was unique
                    return False

        cdef double abs_error
        # predict from existing data, compare to actual value
        if self.n > 0 and not outlier:
            prediction = self.predict_peer(base_ts)
            prediction_error = (prediction - peer_ts) / self.peer_freq

            #if abs(prediction_error) > self.outlier_threshold and abs(prediction_error) > self.error * 4 : # 4 sigma

            if self.n >= 4:
                outlier_threshold = self.outlier_threshold
            else:
                outlier_threshold = 2.0 * self.outlier_threshold

            abs_error = abs(prediction_error)
            self.base.num_syncs += 1
            self.peer.num_syncs += 1
            if abs_error > outlier_threshold:
                if self.peer.bad_syncs < 0.01 and self.base.bad_syncs < 0.01:
                    ac.sync_bad += 1

                if ac.sync_dont_use:
                    return False

                # disable this for the moment
                if 0:
                    if self.peer.bad_syncs < 0.01:
                        self.base.num_outliers += 1
                    if self.base.bad_syncs < 0.01:
                        self.peer.num_outliers += 1

                outlier = True
                self.outlier_total += 1
                if abs_error > 3 * outlier_threshold:
                    self.outliers += 20
                else:
                    self.outliers += 6

                if self.outliers <= 77 and abs_error < 5 * outlier_threshold:
                    return False

                if abs_error > 3 * outlier_threshold:
                    do_reset = True
                    if not self.jumped:
                        if self.peer.bad_syncs < 0.01:
                            self.base.incrementJumps()
                        if self.base.bad_syncs < 0.01:
                            self.peer.incrementJumps()

                    self.jumped = 1
            else:
                ac.sync_good += 1
                ac.last_adsb_time = now

            # disable this
            if 0 and self.n >= 2:
                # wiedehopf: add hacky sync averaging
                # modify new base_ts and peer_ts towards the geometric mean between predition and actual value
                # changing the prediction functions to take into account more past values would likely be the cleaner approach
                # but this modification is significantly easier in regards to the code required
                # so far it seems to be working quite well
                # note that using weight 1/2 so the exact geometric mean seems to be unstable
                # weights 1/4 and 1/3 seem to work well though
                prediction_base = self.predict_base(peer_ts)
                if self.n >= 6 and self.drift_n > drift_n_stable:
                    peer_ts += (prediction - peer_ts) * 0.25
                    base_ts += (prediction_base - base_ts) * 0.25
                else:
                    peer_ts += (prediction - peer_ts) * 0.10
                    base_ts += (prediction_base - base_ts) * 0.10

        if ac.sync_dont_use:
            return False

        cdef double outlier_percent
        if outlier and do_reset:
            if (self.peer.focus and self.base.bad_syncs < 0.01) or (self.base.focus and self.peer.bad_syncs < 0.01):
                outlier_percent = 100.0 * self.outlier_total / self.update_total
                glogger.warning(f'ac {address:06X} step_us {(prediction_error*1e6):.1f} drift_ppm {(self.drift*1e6):.1f} outlier_percent {outlier_percent:.3f} {round(self.outlier_total)}/{round(self.update_total)} pair: {self}')
            #if self.peer.bad_syncs < 0.1 and self.base.bad_syncs < 0.1:
            #   glogger.warning("{r}: {a:06X}: step by {e:.1f}us".format(r=self, a=address, e=prediction_error*1e6))

            # outlier .. we need to reset this clock pair
            self.reset_offsets()
            self.outlier_reset_cooldown = 15 # number of sync pair updates before this sync pair can be trusted
            # as we just reset everything, this is the first point and the prediction error is zero
            prediction_error = 0

        self.outliers = max(0, self.outliers - 18)

        if not outlier:
            self.cumulative_error = max(-50e-6, min(50e-6, self.cumulative_error + prediction_error))  # limit to 50us

        self.outlier_reset_cooldown = max(0, self.outlier_reset_cooldown - 1)

        # update clock drift based on interval ratio
        # this might reject the update
        if not self._update_drift(base_interval, peer_interval):
            self.check_valid(now)
            return False

        self.factor = self.relative_freq * (1 + self.drift)
        self.i_factor = self.i_relative_freq * (1 + self.i_drift)

        # update clock offset based on the actual clock values
        self._update_offset(base_ts, peer_ts, prediction_error)

        self.updated = now
        self.check_valid(now)
        return True

    cdef void _prune_old_data(self, double now):
        cdef int i = 0

        cdef int new_max = cp_size - 12
        if self.n > new_max:
            i = self.n - new_max

        cdef double latest_base_ts = self.ts_base[self.n - 1]
        cdef double limit = 45.0 * self.base_freq
        while i < self.n and (latest_base_ts - self.ts_base[i]) > limit:
            i += 1

        if i > 0:
            self.n -= i
            memmove(self.ts_base, self.ts_base + i, self.n * sizeof(double))
            memmove(self.ts_peer, self.ts_peer + i, self.n * sizeof(double))
            memmove(self.var, self.var + i, self.n * sizeof(double))
            self.var_sum = 0
            for k in range(self.n):
                self.var_sum += self.var[k]
            self.check_valid(now)

    cdef bint _update_drift(self, double base_interval, double peer_interval):
        # try to reduce the effects of catastropic cancellation here:
        #new_drift = (peer_interval / base_interval) / self.relative_freq - 1.0
        cdef double adjusted_base_interval = base_interval * self.relative_freq
        cdef double new_drift = (peer_interval - adjusted_base_interval) / adjusted_base_interval

        if abs(new_drift) > self.drift_max:
            # Bad data, ignore entirely
            #glogger.warn("{0}: drift_max".format(self))
            return False

        if self.drift_n <= 0 or self.drift_outliers > 15:
            # First sample, just trust it outright
            self.raw_drift = self.drift = new_drift
            self.i_drift = -1 * self.drift / (1.0 + self.drift)
            self.cumulative_error = 0.0
            self.drift_outliers = 0
            self.drift_n = 1
            return True

        cdef double drift_error = new_drift - self.raw_drift
        if abs(drift_error) > self.drift_max_delta:
            # Too far away from the value we expect, discard
            #glogger.warn("{0}: drift_max_delta".format(self))
            if self.peer.focus or self.base.focus:
                glogger.warn("{r}: drift_error_ppm out of limits: {de:.1f}".format(r=self, de=1e6*drift_error))
            self.drift_outliers += 1
            return False

        self.drift_outliers = max(0, self.drift_outliers - 2)

        cdef double KP = 0.06
        cdef double KI = 0.008

        # for relatively new pairs allow quicker adjustment of relative drift
        cdef double adjustment_factor
        cdef int unstable = drift_n_stable - self.drift_n
        if unstable > 0:
            adjustment_factor = 1 + (0.2 / KP) * (unstable / drift_n_stable)
            KP *= adjustment_factor

        self.drift_n += 1

        # move towards the new value
        self.raw_drift += drift_error * KP
        self.drift = self.raw_drift - KI * self.cumulative_error
        self.i_drift = -1 * self.drift / (1.0 + self.drift)
        return True

    cpdef void reset_offsets(self):
        self.valid = False
        self.n = 0
        self.var_sum = 0.0
        self.error = -1e-6
        self.variance = -1e-6
        self.outliers = 0
        self.cumulative_error = 0.0

        self.base_avg = 0
        self.peer_avg = 0

    cdef void _update_offset(self, double base_ts, double peer_ts, double prediction_error):
        # insert this into self.ts_base / self.ts_peer / self.var in the right place

        cdef double p_var = prediction_error * prediction_error

        self.ts_base[self.n] = base_ts
        self.ts_peer[self.n] = peer_ts
        self.var[self.n] = p_var

        cdef double elapsed
        cdef double max_elapsed = 20
        cdef double keep

        if self.n == 0:
            self.base_avg = base_ts
            self.peer_avg = peer_ts
        else:
            elapsed = base_ts - self.ts_base[self.n]
            if elapsed < max_elapsed:
                keep = ((1 - elapsed / max_elapsed) ** 2) * 0.75
                self.base_avg = self.base_avg * keep + base_ts * (1 - keep)
                self.peer_avg = self.peer_avg * keep + peer_ts * (1 - keep)
            else:
                self.base_avg = base_ts
                self.peer_avg = peer_ts

        self.n += 1

        self.var_sum += p_var


    cpdef double predict_peer(self, double base_ts):
        """
        Given a time from the base clock, predict the time of the peer clock.
        """

        cdef int n = self.n
        if n == 0:
            raise ValueError("predict_peer called on n == 0 clock pair")

        return self.peer_avg + (base_ts - self.base_avg) * self.factor

    cpdef double predict_base(self, double peer_ts):
        """
        Given a time from the peer clock, predict the time of the base
        clock.
        """

        cdef int n = self.n
        if n == 0:
            raise ValueError("predict_base called on n == 0 clock pair")

        return self.base_avg + (peer_ts - self.peer_avg) * self.i_factor

    def __str__(self):
        return self.base.user + ':' + self.peer.user



"""
Clock normalization routines.
"""

#cdef double _predict(double target_ref, double source_ref, double factor, double source_ts):
#    return target_ref + (source_ts - source_ref) * factor

cdef class _Predictor(object):
    """Simple object for holding prediction state"""

    cdef public double target_ref
    cdef public double source_ref
    cdef public double factor
    cdef public double variance

    def __init__(self, double target_ref, double source_ref, double factor, double variance):
        self.target_ref = target_ref
        self.source_ref = source_ref
        self.factor = factor
        self.variance = variance


#cdef _identity_predict(x):
#    return x


cdef _make_predictors(clock_pairs, station0, station1, double now):
    """Return a tuple of predictors (p_01, p_10) where:

    p_01 will predict a station1 timestamp given a station0 timestamp
    p_10 will predict a station0 timestamp given a station1 timestamp

    Returns None if no suitable clock sync model is available for
    this pair of stations.
    """

    #if station0.epoch is not None and station0.epoch == station1.epoch:
    #    # Assume clocks are closely synchronized to the epoch (and therefore to each other)
    #    predictor = _Predictor(_identity_predict, station0.clock.jitter ** 2 + station1.clock.jitter ** 2)
    #    return (predictor, predictor)

    if station0 < station1:
        pair = clock_pairs.get((station0, station1))
    else:
        pair = clock_pairs.get((station1, station0))

    if pair is None or not pair.valid:
        return None

    cdef ClockPairing pairing = pair

    cdef double variance = pairing.variance
    cdef double stale = now - pairing.updated

    # increase variance for stale pairings
    variance *= 1.0 + stale / 60.0

    # increase variance for pairing with fewer sync points
    if pairing.n < 10:
        variance *= 1 + (10 - pairing.n) * 0.05

    cdef _Predictor base_predictor = _Predictor(pairing.base_avg, pairing.peer_avg, pairing.i_factor, variance)
    cdef _Predictor peer_predictor = _Predictor(pairing.peer_avg, pairing.base_avg, pairing.factor, variance)

    if station0 < station1:
        return (peer_predictor, base_predictor)
    else:
        return (base_predictor, peer_predictor)

cdef _label_heights(g, node, heights):
    """Label each node in the tree with a root of 'node'
    with its height, filling the map 'heights' which
    should be initially empty."""

    # we use heights as a visited-map too.
    heights[node] = 0
    for each in g.neighbors(node):
        if each not in heights:
            _label_heights(g, each, heights)
            mn = heights[each] + g.edge_weight((node, each))
            if mn > heights[node]:
                heights[node] = mn


cdef _tallest_branch(g, node, heights, ignore=None):
    """Find the edge in the tree rooted at 'node' that is part of
    the tallest branch. If ignore is not None, ignore that neighbour.
    Returns (pathlen,node)"""
    tallest = (0, None)

    for each in g.neighbors(node):
        if each is ignore:
            continue

        eh = heights[each] + g.edge_weight((node, each))
        if eh > tallest[0]:
            tallest = (eh, each)

    return tallest


cdef _convert_timestamps(g, timestamp_map, predictor_map, node, results, conversion_chain, variance):
    """Rewrite node and all unvisited nodes reachable from node using the
    chain of clocksync objects in conversion_chain, populating the results dict.

    node: the root node to convert
    timestamp_map: dict of node -> [(timestamp, utc), ...] to convert
    results: dict of node -> (variance, [(converted timestamp, utc), ...])
    conversion_chain: list of predictor tuples to apply to node, in order
    variance: the total error introduced by chain: sum([p.variance for p in chain])
    """

    # convert our own timestamp using the provided chain
    r = []
    results[node] = (variance, r)   # also used as a visited-map
    for ts, utc in timestamp_map[node]:
        for predictor in conversion_chain:
            ts = predictor.target_ref + (ts - predictor.source_ref) * predictor.factor
        r.append((ts, utc))

    # convert all reachable unvisited nodes using a conversion to our timestamp
    # followed by the provided chain
    for neighbor in g.neighbors(node):
        if neighbor not in results:
            predictor = predictor_map[(neighbor, node)]
            _convert_timestamps(g, timestamp_map, predictor_map,
                                neighbor,
                                results,
                                [predictor] + conversion_chain, variance + predictor.variance)


@profile.trackcpu
def normalize(clocktracker, timestamp_map):
    """
    Given {receiver: [(timestamp, utc), ...]}

    return [{receiver: (variance, [(timestamp, utc), ...])}, ...]
    where timestamps are normalized to some arbitrary base timescale within each map;
    one map is returned per connected subgraph."""

    # Represent the stations as a weighted graph where there
    # is an edge between S0 and S1 with weight W if we have a
    # sufficiently recent clock correlation between S0 and S1 with
    # estimated variance W.
    #
    # This graph may have multiple disconnected components. Treat
    # each separately and do this:
    #
    # Find the minimal spanning tree of the component. This will
    # give us the edges to use to convert between timestamps with
    # the lowest total error.
    #
    # Pick a central node of the MST to use as the the timestamp
    # basis, where a central node is a node that minimizes the maximum
    # path cost from the central node to any other node in the spanning
    # tree.
    #
    # Finally, convert all timestamps in the tree to the basis of the
    # central node.

    # populate initial graph
    g = pygraph.classes.graph.graph()
    g.add_nodes(timestamp_map.keys())

    # build a weighted graph where edges represent usable clock
    # synchronization paths, and the weight of each edge represents
    # the estimated variance introducted by converting a timestamp
    # across that clock synchronization.

    # also build a map of predictor objects corresponding to the
    # edges for later use

    now = time.time()

    predictor_map = {}
    for si in timestamp_map.keys():
        for sj in timestamp_map.keys():
            if si < sj:
                predictors = _make_predictors(clocktracker.clock_pairs, si, sj, now)
                if predictors:
                    predictor_map[(si, sj)] = predictors[0]
                    predictor_map[(sj, si)] = predictors[1]
                    g.add_edge((si, sj), wt=predictors[0].variance)

    # find a minimal spanning tree for each component of the graph
    mst_forest = pygraph.algorithms.minmax.minimal_spanning_tree(g)

    # rebuild the graph with only the spanning edges, retaining weights
    # also note the roots of each tree as we go
    g = pygraph.classes.graph.graph()
    g.add_nodes(mst_forest.keys())
    roots = []
    for edge in mst_forest.items():
        if edge[1] is None:
            roots.append(edge[0])
        else:
            g.add_edge(edge, wt=predictor_map[edge].variance)

    # for each spanning tree, find a central node and convert timestamps
    resultComponents = []
    for root in roots:
        # label heights of nodes, where the height of a node is
        # the length of the most expensive path to a child of the node
        heights = {}
        _label_heights(g, root, heights)

        # Find the longest path in the spanning tree; we want to
        # resolve starting at the center of this path, as this minimizes
        # the maximum path length to any node

        # find the two tallest branches leading from the root
        tall1 = _tallest_branch(g, root, heights)
        tall2 = _tallest_branch(g, root, heights, ignore=tall1[1])

        # Longest path is TALL1 - ROOT - TALL2
        # We want to move along the path into TALL1 until the distances to the two
        # tips of the path are equal length. This is the same as finding a node on
        # the path within TALL1 with a height of about half the longest path.
        target = (tall1[0] + tall2[0]) / 2
        central = root
        step = tall1[1]
        while step and abs(heights[central] - target) > abs(heights[step] - target):
            central = step
            _, step = _tallest_branch(g, central, heights, ignore=central)

        # Convert timestamps so they are using the clock units of "central"
        # by walking the spanning tree edges. Then finally convert to wallclock
        # times as the last step by dividing by the final clock's frequency
        results = {}

        factor = 1 / central.clock.freq
        variance = central.clock.jitter**2
        source_ref = 0
        target_ref = 0
        conversion_chain = [_Predictor(target_ref, source_ref, factor, variance)]

        _convert_timestamps(g, timestamp_map, predictor_map, central, results,
                            conversion_chain, central.clock.jitter**2)

        resultComponents.append(results)

    return resultComponents

class _my_component(object):
    """Simple object for holding a graph component"""
    def __init__(self, label, receivers, size):
        self.label = label
        self.receivers = receivers
        self.size = size


@profile.trackcpu
def normalize2(clocktracker, timestamp_map):
    """
    Given {receiver: [(timestamp, utc), ...]}

    return [{receiver: (variance, [(timestamp, utc), ...])}, ...]
    where timestamps are normalized to some arbitrary base timescale within each map;
    one map is returned per connected subgraph."""

    # normalize2: use scipy.sparse.csgraph instead of pythongraph

    # Represent the stations as a weighted graph where there
    # is an edge between S0 and S1 with weight W if we have a
    # sufficiently recent clock correlation between S0 and S1 with
    # estimated variance W.
    #
    # This graph may have multiple disconnected components. Treat
    # each separately and do this:
    #
    # Find the minimal spanning tree of the component. This will
    # give us the edges to use to convert between timestamps with
    # the lowest total error.
    #
    # Pick a central node of the MST to use as the the timestamp
    # basis, where a central node is a node that minimizes the maximum
    # path cost from the central node to any other node in the spanning
    # tree.
    #
    # Finally, convert all timestamps in the tree to the basis of the
    # central node.

    receivers = list(timestamp_map.keys())
    receivers.sort() # to get a matrix with only entries above the diagonal when doing si < sj

    # populate initial graph
    g = pygraph.classes.graph.graph()
    g.add_nodes(receivers)

    # build a weighted graph where edges represent usable clock
    # synchronization paths, and the weight of each edge represents
    # the estimated variance introducted by converting a timestamp
    # across that clock synchronization.

    # also build a map of predictor objects corresponding to the
    # edges for later use

    cdef double now = time.time()

    reclen = len(receivers)
    predictor_count = 0
    predictor_map = {}

    row = []
    col = []
    data = []
    ri = 0
    for si in receivers:
        ci = 0
        for sj in receivers:
            if si < sj:
                predictors = _make_predictors(clocktracker.clock_pairs, si, sj, now)
                if predictors:
                    predictor_map[(si, sj)] = predictors[0]
                    predictor_map[(sj, si)] = predictors[1]
                    g.add_edge((si, sj), wt=predictors[0].variance)
                    predictor_count += 1
                    data.append(predictors[0].variance)
                    row.append(ri)
                    col.append(ci)
            ci += 1
        ri += 1

    if predictor_count < 2:
        return []

    cm = csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(reclen, reclen))

    #for row in cm.toarray():
    #    print(row)

    mst = minimum_spanning_tree(csgraph=cm, overwrite=True)

    n_components, labels = connected_components(csgraph=mst, directed=False, return_labels=True)
    #print('labels: ' + str(labels))

    comps = {}
    for label in labels:
        if label not in comps:
            comps[label] = (_my_component(label=label, receivers=[], size=0))

    index = 0
    for label in labels:
        comps[label].size += 1 # increment component size
        comps[label].receivers.append(receivers[index])
        index += 1

    # make our dict a list so we can sort it
    comps = list(comps.values())
    comps.sort(key=lambda x: x.size, reverse=True)

    if len(comps) == 0:
        return [] # no results, return empty list

    bigComp = comps[0] # biggest component
    roots = [bigComp.receivers[0]] # let's just stay with a list for a moment ... doesn't hurt even if we only do one entry
    #print(bigComp.size)

    if bigComp.size < 3:
        return [] # too small, don't continue

    # rebuild the graph with only the spanning edges, retaining weights
    g = pygraph.classes.graph.graph()
    g.add_nodes(bigComp.receivers)

    coo = mst.tocoo(copy=False)
    for index in range(len(coo.data)):
        weight = coo.data[index]
        si = receivers[coo.row[index]]
        sj = receivers[coo.col[index]]
        if si in bigComp.receivers:
            if si < sj:
                edge = (si,sj)
            else:
                edge = (sj,si)
            #g.add_edge((si, sj), wt=predictor_map[edge].variance)
            g.add_edge((si, sj), wt=weight)

    # for each spanning tree, find a central node and convert timestamps
    # actually we're only searching the biggest spanning tree now
    resultComponents = []
    for root in roots:
        # label heights of nodes, where the height of a node is
        # the length of the most expensive path to a child of the node
        heights = {}
        _label_heights(g, root, heights)

        # Find the longest path in the spanning tree; we want to
        # resolve starting at the center of this path, as this minimizes
        # the maximum path length to any node

        # find the two tallest branches leading from the root
        tall1 = _tallest_branch(g, root, heights)
        tall2 = _tallest_branch(g, root, heights, ignore=tall1[1])

        # Longest path is TALL1 - ROOT - TALL2
        # We want to move along the path into TALL1 until the distances to the two
        # tips of the path are equal length. This is the same as finding a node on
        # the path within TALL1 with a height of about half the longest path.
        target = (tall1[0] + tall2[0]) / 2
        central = root
        step = tall1[1]
        while step and abs(heights[central] - target) > abs(heights[step] - target):
            central = step
            _, step = _tallest_branch(g, central, heights, ignore=central)

        # Convert timestamps so they are using the clock units of "central"
        # by walking the spanning tree edges. Then finally convert to wallclock
        # times as the last step by dividing by the final clock's frequency
        results = {}

        factor = 1 / central.clock.freq
        variance = central.clock.jitter**2
        source_ref = 0
        target_ref = 0
        conversion_chain = [_Predictor(target_ref, source_ref, factor, variance)]

        _convert_timestamps(g, timestamp_map, predictor_map, central, results,
                            conversion_chain, central.clock.jitter**2)

        resultComponents.append(results)

    return resultComponents
