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
Top level glue that knows about all receivers and moves data between
the various sub-objects that make up the server.
"""

import signal
import asyncio
import json
import logging
import time
import os
from contextlib import closing

from mlat import geodesy, profile, constants
from mlat.server import tracker, clocksync, clocktrack, mlattrack, util

glogger = logging.getLogger("coordinator")


class Receiver(object):
    """Represents a particular connected receiver and the associated
    connection that manages it."""

    def __init__(self, uuid, user, connection, clock, position_llh, privacy, connection_info):
        self.uuid = uuid
        self.user = user
        self.connection = connection
        self.clock = clock
        self.position_llh = position_llh
        self.position = geodesy.llh2ecef(position_llh)
        self.privacy = privacy
        self.connection_info = connection_info
        self.dead = False

        self.sync_count = 0
        self.last_rate_report = None
        self.tracking = set()
        self.sync_interest = set()
        self.mlat_interest = set()
        self.requested = set()

        self.distance = {}

        # Receivers with bad_syncs>0 are not used to calculate positions
        self.bad_syncs = 0

    def update_interest_sets(self, new_sync, new_mlat):
        for added in new_sync.difference(self.sync_interest):
            added.sync_interest.add(self)

        for removed in self.sync_interest.difference(new_sync):
            removed.sync_interest.discard(self)

        for added in new_mlat.difference(self.mlat_interest):
            added.mlat_interest.add(self)

        for removed in self.mlat_interest.difference(new_mlat):
            removed.mlat_interest.discard(self)

        self.sync_interest = new_sync
        self.mlat_interest = new_mlat

    @profile.trackcpu
    def refresh_traffic_requests(self):
        self.requested = {x for x in self.tracking if x.interesting}
        self.connection.request_traffic(self, {x.icao for x in self.requested})

    def __lt__(self, other):
        return self.uuid < other.uuid

    def __str__(self):
        return self.uuid

    def __repr__(self):
        return 'Receiver({0!r},{0!r},{1!r})@{2}'.format(self.uuid,
                                                        self.user,
                                                        self.connection,
                                                        id(self))


class Coordinator(object):
    """Master coordinator. Receives all messages from receivers and dispatches
    them to clock sync / multilateration / tracking as needed."""

    def __init__(self, work_dir, partition=(1, 1), tag="mlat", authenticator=None, pseudorange_filename=None):
        """If authenticator is not None, it should be a callable that takes two arguments:
        the newly created Receiver, plus the 'auth' argument provided by the connection.
        The authenticator may modify the receiver if needed. The authenticator should either
        return silently on success, or raise an exception (propagated to the caller) on
        failure.
        """

        self.work_dir = work_dir
        self.receivers = {}    # keyed by uuid
        self.sighup_handlers = []
        self.authenticator = authenticator
        self.partition = partition
        self.tag = tag
        self.tracker = tracker.Tracker(partition)
        self.clock_tracker = clocktrack.ClockTracker()
        self.mlat_tracker = mlattrack.MlatTracker(self,
                                                  blacklist_filename=work_dir + '/blacklist.txt',
                                                  pseudorange_filename=pseudorange_filename)
        self.output_handlers = [self.forward_results]

        self.receiver_mlat = self.mlat_tracker.receiver_mlat
        self.receiver_sync = self.clock_tracker.receiver_sync

    def start(self):
        self._write_state_task = asyncio.async(self.write_state())
        if profile.enabled:
            self._write_profile_task = asyncio.async(self.write_profile())
        else:
            self._write_profile_task = None
        return util.completed_future

    def add_output_handler(self, handler):
        self.output_handlers.append(handler)

    def remove_output_handler(self, handler):
        self.output_handlers.remove(handler)

    # it's a pity that asyncio's add_signal_handler doesn't let you have
    # multiple handlers per signal. so wire up a multiple-handler here.
    def add_sighup_handler(self, handler):
        if not self.sighup_handlers:
            asyncio.get_event_loop().add_signal_handler(signal.SIGHUP, self.sighup)
        self.sighup_handlers.append(handler)

    def remove_sighup_handler(self, handler):
        self.sighup_handlers.remove(handler)
        if not self.sighup_handlers:
            asyncio.get_event_loop().remove_signal_handler(signal.SIGHUP)

    def sighup(self):
        for handler in self.sighup_handlers[:]:
            handler()

    @profile.trackcpu
    def _really_write_state(self):
        aircraft_state = {}
        mlat_count = 0
        sync_count = 0
        now = time.time()
        for ac in self.tracker.aircraft.values():
            s = aircraft_state['{0:06X}'.format(ac.icao)] = {}
            s['interesting'] = 1 if ac.interesting else 0
            s['allow_mlat'] = 1 if ac.allow_mlat else 0
            s['tracking'] = len(ac.tracking)
            s['sync_interest'] = len(ac.sync_interest)
            s['mlat_interest'] = len(ac.mlat_interest)
            s['mlat_message_count'] = ac.mlat_message_count
            s['mlat_result_count'] = ac.mlat_result_count
            s['mlat_kalman_count'] = ac.mlat_kalman_count

            if ac.last_result_time is not None and ac.kalman.valid:
                s['last_result'] = round(now - ac.last_result_time, 1)
                lat, lon, alt = ac.kalman.position_llh
                s['lat'] = round(lat, 3)
                s['lon'] = round(lon, 3)
                s['alt'] = round(alt * constants.MTOF, 0)
                s['heading'] = round(ac.kalman.heading, 0)
                s['speed'] = round(ac.kalman.ground_speed, 0)

            if ac.interesting:
                if ac.sync_interest:
                    sync_count += 1
                if ac.mlat_interest:
                    mlat_count += 1

        if self.partition[1] > 1:
            util.setproctitle('{tag} {i}/{n} ({r} clients) ({m} mlat {s} sync {t} tracked)'.format(
                tag=self.tag,
                i=self.partition[0],
                n=self.partition[1],
                r=len(self.receivers),
                m=mlat_count,
                s=sync_count,
                t=len(self.tracker.aircraft)))
        else:
            util.setproctitle('{tag} ({r} clients) ({m} mlat {s} sync {t} tracked)'.format(
                tag=self.tag,
                r=len(self.receivers),
                m=mlat_count,
                s=sync_count,
                t=len(self.tracker.aircraft)))

        sync = {}
        locations = {}

        for r in self.receivers.values():
            sync[r.uuid] = {
                'peers': self.clock_tracker.dump_receiver_state(r),
                'bad_syncs': r.bad_syncs
            }
            locations[r.uuid] = {
                'user': r.user,
                'lat': r.position_llh[0],
                'lon': r.position_llh[1],
                'alt': r.position_llh[2],
                'privacy': r.privacy,
                'connection': r.connection_info
            }

        # The sync matrix json can be large.  This means it might take a little time to write out.
        # This therefore means someone could start reading it before it has completed writing...
        # So, write out to a temp file first, and then call os.rename(), which is ATOMIC, to overwrite the real file.
        # (Do this for each file, because why not?)
        syncfile = self.work_dir + '/sync.json'
        locationsfile = self.work_dir + '/locations.json'
        aircraftfile = self.work_dir + '/aircraft.json'

        # This random bit can be used for each file
        tmprand = str(int(time.time()))

        # sync.json
        tmpfile = syncfile + '.tmp.' + tmprand
        with closing(open(tmpfile, 'w')) as f:
            json.dump(sync, fp=f, indent=True)
        # We should probably check for errors here, but let's fire-and-forget, instead...
        os.rename(tmpfile, syncfile)

        # locations.json
        tmpfile = locationsfile + '.tmp.' + tmprand
        with closing(open(tmpfile, 'w')) as f:
            json.dump(locations, fp=f, indent=True)
        os.rename(tmpfile, locationsfile)

        # aircraft.json
        tmpfile = aircraftfile + '.tmp.' + tmprand
        with closing(open(tmpfile, 'w')) as f:
            json.dump(aircraft_state, fp=f, indent=True)
        os.rename(tmpfile, aircraftfile)

        # blacklist receivers with bad clock
        # note this section of code runs every 30 seconds
        for r in self.receivers.values():
            bad_peers = 0
            # count how many peers we have bad sync with
            # don't count peers who have been timed out (state[4] > 0)
            # 3 microseconds error or more are considered a bad sync (state[1] > 3)
            num_peers = 10
            # start with 10 peers extra, so low peer receivers
            # aren't timed out by the percentage threshold
            # of bad_peers as easily.

            # iterate over sync state with all peers
            # state = [ 0: pairing sync count, 1: offset, 2: drift,
            #           3: timestamp?, 4: bad_syncs ]
            for state in sync[r.uuid].peers.values():
                if state[4] > 0:
                    continue
                num_peers += 1
                if state[1] > 3:
                    bad_peers += 1

            # If your sync with 5 receivers or more than 10 percent of peers is bad,
            # it's likely you are the reason.
            # You get 0.2 to 1 to your bad_sync score and timed out.

            if bad_peers > 5 or bad_peers/num_peers > 0.1:
                r.bad_syncs += min(1, 2*bad_peers/num_peers)
            else:
                r.bad_syncs -= 0.1

            # If your sync mostly looks good, your bad_sync score is decreased.
            # If you had a score before, once it goes down to zero you are
            # no longer timed out

            # Limit bad_sync score to the range of 0 to 6

            r.bad_syncs = max(0, min(6, r.bad_syncs))


    @asyncio.coroutine
    def write_state(self):
        while True:
            try:
                self._really_write_state()
            except Exception:
                glogger.exception("Failed to write state files")

            yield from asyncio.sleep(30.0)

    @asyncio.coroutine
    def write_profile(self):
        while True:
            yield from asyncio.sleep(60.0)

            try:
                with closing(open(self.work_dir + '/cpuprofile.txt', 'w')) as f:
                    profile.dump_cpu_profiles(f)
            except Exception:
                glogger.exception("Failed to write CPU profile")

    def close(self):
        self._write_state_task.cancel()
        if self._write_profile_task:
            self._write_profile_task.cancel()

    @asyncio.coroutine
    def wait_closed(self):
        util.safe_wait([self._write_state_task, self._write_profile_task])

    @profile.trackcpu
    def new_receiver(self, connection, uuid, user, auth, position_llh, clock_type, privacy, connection_info):
        """Assigns a new receiver ID for a given user.
        Returns the new receiver ID.

        May raise ValueError to disallow this receiver."""

        if uuid in self.receivers:
            raise ValueError('User {uuid}/{user} is already connected'.format(uuid=uuid, user=user))

        clock = clocksync.make_clock(clock_type)
        receiver = Receiver(uuid, user, connection, clock,
                            position_llh=position_llh,
                            privacy=privacy,
                            connection_info=connection_info)

        if self.authenticator is not None:
            self.authenticator(receiver, auth)  # may raise ValueError if authentication fails

        self._compute_interstation_distances(receiver)

        self.receivers[receiver.uuid] = receiver
        return receiver

    def _compute_interstation_distances(self, receiver):
        """compute inter-station distances for a receiver"""

        for other_receiver in self.receivers.values():
            if other_receiver is receiver:
                distance = 0
            else:
                distance = geodesy.ecef_distance(receiver.position, other_receiver.position)
            receiver.distance[other_receiver] = distance
            other_receiver.distance[receiver] = distance

    @profile.trackcpu
    def receiver_location_update(self, receiver, position_llh):
        """Note that a given receiver has moved."""
        receiver.position_llh = position_llh
        receiver.position = geodesy.llh2ecef(position_llh)

        self._compute_interstation_distances(receiver)

    @profile.trackcpu
    def receiver_disconnect(self, receiver):
        """Notes that the given receiver has disconnected."""

        receiver.dead = True
        self.tracker.remove_all(receiver)
        self.clock_tracker.receiver_disconnect(receiver)
        self.receivers.pop(receiver.uuid)

        # clean up old distance entries
        for other_receiver in self.receivers.values():
            other_receiver.distance.pop(receiver, None)

    @profile.trackcpu
    def receiver_tracking_add(self, receiver, icao_set):
        """Update a receiver's tracking set by adding some aircraft."""
        self.tracker.add(receiver, icao_set)
        if receiver.last_rate_report is None:
            # not receiving rate reports for this receiver
            self.tracker.update_interest(receiver)

    @profile.trackcpu
    def receiver_tracking_remove(self, receiver, icao_set):
        """Update a receiver's tracking set by removing some aircraft."""
        self.tracker.remove(receiver, icao_set)
        if receiver.last_rate_report is None:
            # not receiving rate reports for this receiver
            self.tracker.update_interest(receiver)

    @profile.trackcpu
    def receiver_clock_reset(self, receiver):
        """Reset current clock synchronization for a receiver."""
        self.clock_tracker.receiver_clock_reset(receiver)

    @profile.trackcpu
    def receiver_rate_report(self, receiver, report):
        """Process an ADS-B position rate report for a receiver."""
        receiver.last_rate_report = report
        self.tracker.update_interest(receiver)

    @profile.trackcpu
    def forward_results(self, receive_timestamp, address, ecef, ecef_cov, receivers, distinct, dof, kalman_state):
        broadcast = receivers
        ac = self.tracker.aircraft.get(address)
        if ac:
            ac.successful_mlat.update(receivers)
            broadcast = ac.successful_mlat
        for receiver in broadcast:
            try:
                receiver.connection.report_mlat_position(receiver,
                                                         receive_timestamp, address,
                                                         ecef, ecef_cov, receivers, distinct,
                                                         dof, kalman_state)
            except Exception:
                glogger.exception("Failed to forward result to receiver {r}".format(r=receiver.uuid))
                # eat the exception so it doesn't break our caller
