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

import random
import signal
import asyncio
import ujson
import logging
import logging.handlers
import time
import os
from contextlib import closing
import array

from mlat import geodesy, profile, constants
from mlat import tracker, clocktrack, mlattrack, util, config

glogger = logging.getLogger("coordinator")
random.seed()


class Receiver(object):
    """Represents a particular connected receiver and the associated
    connection that manages it."""

    def __init__(self, uid, user, connection, clock_type, position_llh, privacy, connection_info, uuid, coordinator, clock_tracker):
        self.uid = uid
        self.uuid = uuid
        self.user = user
        self.connection = connection
        self.coordinator = coordinator
        self.clock_tracker = clock_tracker
        self.clock = clocktrack.make_clock(clock_type)
        self.epoch = None
        if clock_type == 'radarcape_gps':
            self.epoch = 'gps_midnight'
        self.last_clock_reset = time.time()
        self.clock_reset_counter = 0

        self.position_llh = position_llh
        self.position = geodesy.llh2ecef(position_llh)

        self.privacy = privacy
        self.connection_info = connection_info
        self.dead = False
        self.connectedSince = time.time()

        self.num_outliers = 0
        self.num_syncs = 0

        self.outlier_percent_rolling = 0

        self.sync_peers = array.array('i', [0, 0, 0, 0, 0]) # number of peers per distance category
        self.peer_count = 0 # only updated when dumping state
        self.last_rate_report = None
        self.tracking = set()
        self.adsb_seen = set()
        self.sync_interest = set()
        self.mlat_interest = set()
        self.requested = set()
        self.mapLat = 0
        self.mapLon = 0
        self.mapAlt = 0

        self.distance = {}

        # timestamp this receiver last synced with the result being a valid clock pair
        self.last_sync = 0

        # Receivers with bad_syncs>0 are not used to calculate positions
        self.bad_syncs = 0
        self.sync_range_exceeded = 0

        self.recent_pair_jumps = 0

        self.focus = False

        self.coordinator.receiver_location_update(self, position_llh)

    def update_interest_sets(self, new_sync, new_mlat, new_adsb):

        if self.bad_syncs > 3 and len(new_sync) > 3:
            new_sync = set(random.sample(list(new_sync), k=3))

        if self.bad_syncs > 0:
            new_mlat = set()

        for added in new_adsb.difference(self.adsb_seen):
            added.adsb_seen.add(self)

        for removed in self.adsb_seen.difference(new_adsb):
            removed.adsb_seen.discard(self)


        for added in new_sync.difference(self.sync_interest):
            added.sync_interest.add(self)

        for removed in self.sync_interest.difference(new_sync):
            removed.sync_interest.discard(self)

        for added in new_mlat.difference(self.mlat_interest):
            added.mlat_interest.add(self)

        for removed in self.mlat_interest.difference(new_mlat):
            removed.mlat_interest.discard(self)

        self.adsb_seen = new_adsb
        self.sync_interest = new_sync
        self.mlat_interest = new_mlat

    def incrementJumps(self):
        self.recent_pair_jumps += 1
        total_peers = sum(self.sync_peers)
        if total_peers > 0 and self.recent_pair_jumps / total_peers > 0.25 and self.recent_pair_jumps > 3:
            self.recent_pair_jumps = 0
            self.clock_reset('sync detected clock jump')
            if self.focus:
                glogger.warning("{r}: detected clockjump, bad_syncs: {b}".format(r=self.user, b=round(self.bad_syncs, 1)))
            self.bad_syncs += 0.1

    def clock_reset(self, reason):
        """Reset current clock synchronization for this receiver."""
        self.clock_tracker.receiver_clock_reset(self)
        self.last_clock_reset = time.time()
        self.clock_reset_counter += 1
        if self.focus:
            if not reason:
                reason = ''
            glogger.warning("Clock reset: {r} count: {c} reason: {s}".format(r=self.user, c=self.clock_reset_counter, s=reason))

    @profile.trackcpu
    def refresh_traffic_requests(self):
        self.requested = self.sync_interest | self.mlat_interest
        self.connection.request_traffic(self, {x.icao for x in self.requested})

    def __lt__(self, other):
        return self.uid < other.uid

    def __str__(self):
        return self.user

    def __repr__(self):
        return 'Receiver({0!r},{0!r},{1!r})@{2}'.format(self.uid,
                                                        self.user,
                                                        self.connection,
                                                        id(self))


class Coordinator(object):
    """Master coordinator. Receives all messages from receivers and dispatches
    them to clock sync / multilateration / tracking as needed."""

    def __init__(self, work_dir, loop, status_interval, partition=(1, 1), tag="mlat", authenticator=None, pseudorange_filename=None):
        """If authenticator is not None, it should be a callable that takes two arguments:
        the newly created Receiver, plus the 'auth' argument provided by the connection.
        The authenticator may modify the receiver if needed. The authenticator should either
        return silently on success, or raise an exception (propagated to the caller) on
        failure.
        """

        self.work_dir = work_dir
        self.loop = loop

        self.uidCounter = 0
        # receivers:
        self.receivers = {} # keyed by uid
        self.usernames = {} # keyed by usernames

        self.sighup_handlers = []
        self.authenticator = authenticator
        self.partition = partition
        self.tag = tag
        self.tracker = tracker.Tracker(self, partition, loop)
        self.clock_tracker = clocktrack.ClockTracker(self, loop)
        self.mlat_tracker = mlattrack.MlatTracker(self,
                                                  blacklist_filename=work_dir + '/blacklist.txt',
                                                  pseudorange_filename=pseudorange_filename)
        self.output_handlers = []

        self.receiver_mlat = self.mlat_tracker.receiver_mlat
        self.receiver_sync = self.clock_tracker.receiver_sync


        self.handshake_logger = logging.getLogger("handshake")
        self.handshake_logger.setLevel(logging.DEBUG)

        self.handshake_handler = logging.handlers.RotatingFileHandler(
                (self.work_dir + '/handshakes.log'),
                maxBytes=(1*1024*1024), backupCount=2)

        self.handshake_logger.addHandler(self.handshake_handler)

        self.main_interval = 15.0

        self.last_cpu_time = 0
        self.stats_sync_points = 0
        self.stats_sync_msgs = 0
        self.stats_mlat_msgs = 0
        self.stats_valid_groups = 0
        self.stats_normalize = 0
        self.stats_solve_attempt = 0
        self.stats_solve_success = 0
        self.stats_solve_used = 0

        if status_interval is None:
            status_interval = 15
        status_interval = float(status_interval)
        if status_interval < 0:
            self.status_interval = 1e15 # a really long time
        else:
            self.status_interval = status_interval * 0.95

        self.next_status = time.time() + self.status_interval

    def start(self):
        self._every_15_task = asyncio.ensure_future(self.every_15())
        if profile.enabled:
            self._write_profile_task = asyncio.ensure_future(self.write_profile())
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
            self.loop.add_signal_handler(signal.SIGHUP, self.sighup)
        self.sighup_handlers.append(handler)

    def remove_sighup_handler(self, handler):
        self.sighup_handlers.remove(handler)
        if not self.sighup_handlers:
            self.loop.remove_signal_handler(signal.SIGHUP)

    def sighup(self):
        for handler in self.sighup_handlers[:]:
            handler()

    @profile.trackcpu
    def _write_state(self):
        aircraft_state = {}
        ac_count_mlat = len(self.tracker.mlat_wanted)
        ac_count_sync = 0
        now = time.time()
        for ac in self.tracker.aircraft.values():
            elapsed_seen = now - ac.seen
            #if elapsed_seen > 3 * 3600:
            #    continue
            icao_string = '{0:06X}'.format(ac.icao)
            s = {}
            aircraft_state[icao_string] = s
            s['icao'] = icao_string
            s['elapsed_seen'] = round(elapsed_seen, 1)
            s['interesting'] = 1 if ac.interesting else 0
            s['allow_mlat'] = 1 if ac.allow_mlat else 0
            s['tracking'] = len(ac.tracking)
            s['sync_interest'] = len(ac.sync_interest)
            s['mlat_interest'] = len(ac.mlat_interest)
            s['adsb_seen'] = len(ac.adsb_seen)
            s['mlat_message_count'] = ac.mlat_message_count
            s['mlat_result_count'] = ac.mlat_result_count
            s['mlat_kalman_count'] = ac.mlat_kalman_count


            sync_count = round(ac.sync_good + ac.sync_bad)
            s['sync_count_1min'] = sync_count

            sync_bad_percent = round(100 * ac.sync_bad / (sync_count + 0.01), 1)
            s['sync_bad_percent'] = sync_bad_percent
            ac.sync_bad_percent = sync_bad_percent

            if ac.sync_bad > 3 and sync_bad_percent > 15:
                ac.sync_dont_use = 1
            else:
                ac.sync_dont_use = 0


            ac.sync_good *= 0.8
            ac.sync_bad *= 0.8

            if ac.last_result_time is not None:
                s['last_result'] = round(now - ac.last_result_time, 1)
                lat, lon, alt = geodesy.ecef2llh(ac.last_result_position)
                s['lat'] = round(lat, 4)
                s['lon'] = round(lon, 4)
                alt = ac.altitude
                if alt is not None:
                    s['alt'] = round(alt)
                if ac.kalman.valid:
                    s['heading'] = round(ac.kalman.heading, 0)
                    s['speed'] = round(ac.kalman.ground_speed, 0)

            if elapsed_seen > 600:
                s['tracking_receivers'] = [receiver.uid for receiver in ac.tracking]

            if ac.sync_interest:
                ac_count_sync += 1

        sync = {}
        clients = {}

        receiver_states = self.clock_tracker.dump_receiver_state()

        bad_receivers = 0

        outlier_sum = 0
        sync_sum = 0

        # blacklist receivers with bad clock
        # note this section of code runs every 15 seconds
        for r in self.receivers.values():
            bad_peers = 0
            # count how many peers we have bad sync with
            # don't count peers who have been timed out (state[3] > 0)

            num_peers = 8

            # iterate over sync state with all peers
            # state
            # 0: pairing sync count
            # 1: offset
            # 2: drift
            # 3: bad_syncs
            # 4: pairing.jumped
            # 5: pairing.outlier_percent
            # 6: time since last pairing update
            # 7: how many synced peers does the peer have

            peers = receiver_states.get(r.user, {})
            bad_peer_list = []
            sum_outlier_percent = 0
            for username, state in peers.items():
                if state[3] > 0 or state[6] > 20 or state[7] < 8:
                    # skip peers which have bad sync
                    # skip peers which haven't updated in a bit
                    continue
                sum_outlier_percent += state[5]
                num_peers += 1
                if state[4] or state[5] > 35 or (state[0] > 10 and state[1] > 1.2) or (state[0] > 3 and state[1] > 1.8) or state[1] > 2.4:
                    bad_peers += 1
                    bad_peer_list.append(username)


            outlier_sum += sum_outlier_percent / 100
            sync_sum += num_peers
            #outlier_percent = 100 * r.num_outliers / (r.num_syncs + 0.1)
            outlier_percent = sum_outlier_percent / num_peers

            # running average for the outlier percent
            r.outlier_percent_rolling -= 0.1 * (r.outlier_percent_rolling - outlier_percent)

            # If your sync with more than 10 percent of peers is bad,
            # it's likely you are the reason.
            # You get 0.5 to 2 to your bad_sync score and timed out.

            if bad_peers / num_peers > 0.15 and bad_peers > 3:
                r.bad_syncs += min(0.5, 2*bad_peers/num_peers) + 0.1

            r.bad_syncs -= 0.1

            # If your sync mostly looks good, your bad_sync score is decreased.
            # If you had a score before, once it goes down to zero you are
            # no longer timed out

            # Limit bad_sync score to the range of 0 to 6

            r.bad_syncs = max(0, min(6, r.bad_syncs))

            if r.bad_syncs > 0:
                bad_receivers += 1

            if r.focus:
                glogger.warning("{u}: bad_syncs: {bs:0.1f} outlier percent: {pe:0.1f} bad peers: {bp} ratio: {r} list: {l}".format(
                    u=r.user, bs=r.bad_syncs, pe=r.outlier_percent_rolling,
                    bp=bad_peers, r=round(bad_peers/num_peers, 2), l=str(bad_peer_list)))


            r.recent_pair_jumps = 0

            # almost reset num_outlier / num_syns for each receiver, keep a bit of the last iteration
            r.num_outliers *= 0.25
            r.num_syncs *= 0.25

            # r.mapLat / r.mapLon is a fudged position for privacy
            sync[r.user] = {
                'peers': peers,
                'bad_syncs': r.bad_syncs,
                'lat': r.mapLat,
                'lon': r.mapLon
            }

            r.peer_count = len(sync[r.user]['peers'])

            clients[r.user] = {
                'user': r.user,
                'uid': r.uid,
                'uuid': r.uuid,
                'coords': "{0:.6f},{1:.6f}".format(r.position_llh[0], r.position_llh[1]),
                'lat': r.position_llh[0],
                'lon': r.position_llh[1],
                'alt': r.position_llh[2],
                'privacy': r.privacy,
                'connection': r.connection_info,
                'source_ip': r.connection.source_ip,
                'source_port': r.connection.source_port,
                'message_rate': round(r.connection.message_counter / 15.0),
                'peer_count': len(peers),
                'bad_sync_timeout': round(r.bad_syncs * 15 / 0.1),
                'outlier_percent': round(r.outlier_percent_rolling, 1),
                'bad_peer_list': str(bad_peer_list),
                'sync_interest': [format(a.icao, '06x') for a in r.sync_interest],
                'mlat_interest': [format(a.icao, '06x') for a in r.mlat_interest]
            }

            statistics = {
                'peer_count': len(peers),
                'bad_sync_timeout': round(r.bad_syncs * 15 / 0.1),
                'outlier_percent': round(r.outlier_percent_rolling, 1)
            }

            r.connection.send_stats(statistics)

            # reset message counter
            r.connection.message_counter = 0

        # The sync matrix json can be large.  This means it might take a little time to write out.
        # This therefore means someone could start reading it before it has completed writing...
        # So, write out to a temp file first, and then call os.replace(), which is ATOMIC, to overwrite the real file.
        # (Do this for each file, because why not?)
        syncfile = self.work_dir + '/sync.json'
        clientsfile = self.work_dir + '/clients.json'
        aircraftfile = self.work_dir + '/aircraft.json'

        # sync.json
        tmpfile = syncfile + '.tmp'
        with closing(open(tmpfile, 'w')) as f:
            ujson.dump(sync, f)
        # We should probably check for errors here, but let's fire-and-forget, instead...
        os.replace(tmpfile, syncfile)

        # clients.json
        tmpfile = clientsfile + '.tmp'
        with closing(open(tmpfile, 'w')) as f:
            ujson.dump(clients, f)
        os.replace(tmpfile, clientsfile)

        # aircraft.json
        tmpfile = aircraftfile + '.tmp'
        with closing(open(tmpfile, 'w')) as f:
            ujson.dump(aircraft_state, f)
        os.replace(tmpfile, aircraftfile)

        total_outlier_percent = 100 * outlier_sum / (sync_sum + 0.1)

        cpu_time = time.clock_gettime(time.CLOCK_PROCESS_CPUTIME_ID)
        cpu_time_us_per_sec = round((cpu_time - self.last_cpu_time) * (1e6 / 15))
        self.last_cpu_time = cpu_time
        try:
            with open('/run/node_exporter/mlat-server.prom', 'w', encoding='utf-8') as f:
                out = ''
                out += 'mlat_server_cpu_ppm ' + str(cpu_time_us_per_sec) + '\n'
                out += 'mlat_server_receivers ' + str(len(self.receivers)) + '\n'
                out += 'mlat_server_ac_mlat ' + str(ac_count_mlat) + '\n'
                out += 'mlat_server_ac_sync ' + str(ac_count_sync) + '\n'
                out += 'mlat_server_ac_total ' + str(len(self.tracker.aircraft)) + '\n'
                out += 'mlat_server_outlier_ppm ' + "{0:.0f}".format(total_outlier_percent * 1000) + '\n'
                out += 'mlat_server_sync_points ' + "{0:.0f}".format(self.stats_sync_points / self.main_interval) + '\n'
                out += 'mlat_server_sync_msgs ' + "{0:.0f}".format(self.stats_sync_msgs / self.main_interval) + '\n'
                out += 'mlat_server_mlat_msgs ' + "{0:.0f}".format(self.stats_mlat_msgs / self.main_interval) + '\n'
                out += 'mlat_server_valid_groups ' + "{0:.0f}".format(self.stats_valid_groups / self.main_interval) + '\n'
                out += 'mlat_server_normalize_called ' + "{0:.0f}".format(self.stats_normalize / self.main_interval) + '\n'
                out += 'mlat_server_solve_attempt ' + "{0:.0f}".format(self.stats_solve_attempt / self.main_interval) + '\n'
                out += 'mlat_server_solve_success ' + "{0:.0f}".format(self.stats_solve_success / self.main_interval) + '\n'
                out += 'mlat_server_solve_used ' + "{0:.0f}".format(self.stats_solve_used / self.main_interval) + '\n'

                f.write(out)
        except OSError:
            pass
        except:
            glogger.exception("prom stats")

        # reset stats
        self.stats_sync_points = 0
        self.stats_sync_msgs = 0
        self.stats_mlat_msgs = 0
        self.stats_valid_groups = 0
        self.stats_normalize = 0
        self.stats_solve_attempt = 0
        self.stats_solve_success = 0
        self.stats_solve_used = 0


        if self.partition[1] > 1:
            title_string = 'Status: {i}/{n} ({r} clients) ({m} mlat {s} sync {t} tracked)'.format(
                i=self.partition[0],
                n=self.partition[1],
                r=len(self.receivers),
                m=ac_count_mlat,
                s=ac_count_sync,
                t=len(self.tracker.aircraft))
        else:
            title_string = 'Status: ({r} clients {b} bad sync) ({o:.2f} outlier_percentage) ({m} mlat {s} sync {t} tracked)'.format(
                r=len(self.receivers),
                b=bad_receivers,
                o=total_outlier_percent,
                m=ac_count_mlat,
                s=ac_count_sync,
                t=len(self.tracker.aircraft))
        util.setproctitle(title_string)

        if now > self.next_status:
            self.next_status = now + self.status_interval
            glogger.warning(title_string)


    async def every_15(self):
        while True:
            sleep = asyncio.create_task(asyncio.sleep(self.main_interval))
            try:
                self._write_state()
                self.clock_tracker.clear_all_sync_points()
            except Exception:
                glogger.exception("Failed to write state files")

            await sleep

    async def write_profile(self):
        while True:
            await asyncio.sleep(60.0)

            try:
                with closing(open(self.work_dir + '/cpuprofile.txt', 'w')) as f:
                    profile.dump_cpu_profiles(f)
            except Exception:
                glogger.exception("Failed to write CPU profile")

    def close(self):
        self._every_15_task.cancel()
        if self._write_profile_task:
            self._write_profile_task.cancel()

    async def wait_closed(self):
        await util.safe_wait([self._every_15_task, self._write_profile_task])

    @profile.trackcpu
    def new_receiver(self, connection, uuid, user, auth, position_llh, clock_type, privacy, connection_info):
        """Assigns a new receiver ID for a given user.
        Returns the new receiver.

        May raise ValueError to disallow this receiver."""

        if user in self.usernames:
            raise ValueError('User {user} is already connected'.format(user=user))

        if self.uidCounter > 4611686018427387904:
            self.uidCounter = 0
        uid = self.uidCounter
        while uid in self.receivers:
            self.uidCounter += 1
            uid = self.uidCounter

        receiver = Receiver(uid, user, connection, clock_type,
                            position_llh=position_llh,
                            privacy=privacy,
                            connection_info=connection_info,
                            uuid=uuid,
                            coordinator=self,
                            clock_tracker=self.clock_tracker)

        if self.authenticator is not None:
            self.authenticator(receiver, auth)  # may raise ValueError if authentication fails

        self._compute_interstation_distances(receiver)

        self.receivers[receiver.uid] = receiver
        self.usernames[receiver.user] = receiver
        if receiver.user.startswith(config.DEBUG_FOCUS):
            receiver.focus = True
        return receiver

    def _compute_interstation_distances(self, receiver):
        """compute inter-station distances for a receiver"""

        for other_receiver in self.receivers.values():
            if other_receiver is receiver:
                distance = 0
            else:
                distance = geodesy.ecef_distance(receiver.position, other_receiver.position)
            receiver.distance[other_receiver.uid] = distance
            other_receiver.distance[receiver.uid] = distance

    @profile.trackcpu
    def receiver_location_update(self, receiver, position_llh):
        """Note that a given receiver has moved."""
        receiver.position_llh = position_llh
        receiver.position = geodesy.llh2ecef(position_llh)


        # fudge map position, set retained precision as a fraction of a degree:
        precision = 20
        r = receiver
        if r.privacy:
            r.mapLat = None
            r.mapLon = None
            r.mapAlt = None
        else:

            offX = -1/precision + 1/precision * random.random()
            offY = -1/precision + 1/precision * random.random()
            r.mapLat = round(round(r.position_llh[0] * precision) / precision + offX, 2)
            r.mapLon = round(round(r.position_llh[1] * precision) / precision + offY, 2)
            r.mapAlt = 50 * round(r.position_llh[2]/50)

        self._compute_interstation_distances(receiver)

    @profile.trackcpu
    def receiver_disconnect(self, receiver):
        """Notes that the given receiver has disconnected."""

        receiver.dead = True
        self.tracker.remove_all(receiver)
        self.clock_tracker.receiver_disconnect(receiver)
        self.receivers.pop(receiver.uid)
        self.usernames.pop(receiver.user)

        # clean up old distance entries
        for other_receiver in self.receivers.values():
            other_receiver.distance.pop(receiver.uid, None)

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
    def receiver_rate_report(self, receiver, report):
        """Process an ADS-B position rate report for a receiver."""
        receiver.last_rate_report = report
        self.tracker.update_interest(receiver)

    @profile.trackcpu
    def forward_results(self, receive_timestamp, address, ecef, ecef_cov, receivers, distinct, dof, kalman_state, error):

        # don't forward if kalman hasn't locked on and it's only 3 receivers
        if not kalman_state.valid and dof < 1:
            return

        broadcast = receivers
        # only send result to receivers who received this message
        result_new_old = [ None, None ]
        for receiver in broadcast:
            try:
                receiver.connection.report_mlat_position(receiver,
                                                         receive_timestamp, address,
                                                         ecef, ecef_cov, receivers, distinct,
                                                         dof, kalman_state, result_new_old)
            except Exception:
                glogger.exception("Failed to forward result to receiver {r}".format(r=receiver.user))
                # eat the exception so it doesn't break our caller
