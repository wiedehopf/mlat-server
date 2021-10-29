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
Maintains clock synchronization between individual pairs of receivers.
"""

import math
import time
import bisect
import logging

__all__ = ('Clock', 'ClockPairing', 'make_clock')

glogger = logging.getLogger("clocksync")


class Clock(object):
    """A particular clock. Stores characteristics of a clock,
    and acts as part of the key in the clock pairing map.
    """

    def __init__(self, epoch, freq, max_freq_error, jitter):
        """Create a new clock representation.

        epoch: a string indicating a fixed epoch, or None if freerunning
        freq: the clock frequency in Hz (float)
        max_freq_error: the maximum expected relative frequency error (i.e. 1e-6 is 1PPM) (float)
        jitter: the expected jitter of a typical reading, in seconds, standard deviation  (float)
        """
        self.epoch = epoch
        self.freq = freq
        self.max_freq_error = max_freq_error
        self.jitter = jitter


def make_clock(clock_type):
    """Return a new Clock instance for the given clock type."""

    if clock_type == 'radarcape_gps':
        return Clock(epoch='gps_midnight', freq=1e9, max_freq_error=1e-6, jitter=15e-9)
    if clock_type == 'beast' or clock_type == 'radarcape_12mhz':
        return Clock(epoch=None, freq=12e6, max_freq_error=5e-6, jitter=83e-9)
    if clock_type == 'sbs':
        return Clock(epoch=None, freq=20e6, max_freq_error=100e-6, jitter=500e-9)
    if clock_type == 'dump1090' or clock_type == 'unknown':
        return Clock(epoch=None, freq=12e6, max_freq_error=100e-6, jitter=500e-9)
    raise NotImplementedError("{ct}".format(ct=clock_type))


class ClockPairing(object):
    """Describes the current relative characteristics of a pair of clocks."""

    KP = 0.05
    KI = 0.01

    def __init__(self, base, peer, cat):
        self.base = base
        self.peer = peer
        self.cat = cat
        self.base_clock = base.clock
        self.peer_clock = peer.clock
        self.raw_drift = None
        self.drift = None
        self.i_drift = None
        self.n = 0
        self.ts_base = []
        self.ts_peer = []
        self.var = []
        self.var_sum = 0.0
        self.outliers = 0
        self.cumulative_error = 0.0
        self.error = None
        self.variance = None

        self.jumped = 0

        self.relative_freq = peer.clock.freq / base.clock.freq
        self.i_relative_freq = base.clock.freq / peer.clock.freq
        self.drift_max = base.clock.max_freq_error + peer.clock.max_freq_error
        self.drift_max_delta = self.drift_max / 10.0
        self.outlier_threshold = 5 * math.sqrt(peer.clock.jitter ** 2 + base.clock.jitter ** 2) # 5 sigma

        now = time.monotonic()
        self.expiry = now + 30.0
        self.validity = now + 25.0
        self.updated = now
        self.valid = False


    def updateVars(self):
        if self.n == 0:
            self.variance = None
            self.error = None
        else:
            """Variance of recent predictions of the sync point versus the actual sync point."""
            self.variance = self.var_sum / self.n

            """Standard error of recent predictions."""
            self.error = math.sqrt(self.variance)

    def check_valid(self, now):
        """True if this pairing is usable for clock syncronization."""
        return (self.n >= 3 and (self.var_sum / self.n) < 16e-12 and
                    self.outliers < 3 and self.validity > now)

    def update(self, address, base_ts, peer_ts, base_interval, peer_interval, now):
        """Update the relative drift and offset of this pairing given:

        address: the ICAO address of the sync aircraft, for logging purposes
        base_ts: the timestamp of a recent point in time measured by the base clock
        peer_ts: the timestamp of the same point in time measured by the peer clock
        base_interval: the duration of a recent interval measured by the base clock
        peer_interval: the duration of the same interval measured by the peer clock

        Returns True if the update was used, False if it was an outlier.
        """

        if self.n != 0 and base_ts <= self.ts_base[-1]:
            # timestamp is in the past or duplicated, don't use this
            return False

        # clean old data
        if self.n > 30 or (self.n > 1 and (base_ts - self.ts_base[0]) > 55 * self.base_clock.freq):
            self._prune_old_data(base_ts)

        outlier = False
        # predict from existing data, compare to actual value
        if self.n > 0:
            prediction = self.predict_peer(base_ts)
            prediction_error = (prediction - peer_ts) / self.peer_clock.freq

            if (abs(prediction_error) > self.outlier_threshold or self.n > 8) and abs(prediction_error) > self.error * 4 : # 4 sigma
                outlier = True
                self.outliers += 1
                self.outlier = min(7, self.outliers)
                if self.outliers < 5:
                    # don't accept this one
                    self.valid = self.check_valid(now)
                    return False
            else:
                # wiedehopf: add hacky sync averaging
                # modify new base_ts and peer_ts towards the geometric mean between predition and actual value
                # changing the prediction functions to take into account more past values would likely be the cleaner approach
                # but this modification is significantly easier in regards to the code required
                # so far it seems to be working quite well
                # note that using weight 1/2 so the exact geometric mean seems to be unstable
                # weights 1/4 and 1/3 seem to work well though
                prediction_base = self.predict_base(peer_ts)
                peer_ts += (prediction - peer_ts) / 3
                base_ts += (prediction_base - base_ts) / 3
        else:
            prediction_error = 0  # first sync point, no error

        # update clock drift based on interval ratio
        # this might reject the update
        if not self._update_drift(address, base_interval, peer_interval):
            self.valid = self.check_valid(now)
            return False

        # update clock offset based on the actual clock values
        self._update_offset(address, base_ts, peer_ts, prediction_error, outlier)

        self.expiry = now + 45.0
        self.validity = now + 30.0
        self.updated = now
        self.valid = self.check_valid(now)
        return True

    def _prune_old_data(self, latest_base_ts):
        i = 0

        if self.n > 20:
            i = self.n - 20

        while i < self.n and (latest_base_ts - self.ts_base[i]) > 45 * self.base_clock.freq:
            i += 1

        if i > 0:
            del self.ts_base[0:i]
            del self.ts_peer[0:i]
            del self.var[0:i]
            self.n -= i
            self.var_sum = sum(self.var)
            self.updateVars()

    def _update_drift(self, address, base_interval, peer_interval):
        # try to reduce the effects of catastropic cancellation here:
        #new_drift = (peer_interval / base_interval) / self.relative_freq - 1.0
        adjusted_base_interval = base_interval * self.relative_freq
        new_drift = (peer_interval - adjusted_base_interval) / adjusted_base_interval

        if abs(new_drift) > self.drift_max:
            # Bad data, ignore entirely
            return False

        if self.drift is None:
            # First sample, just trust it outright
            self.raw_drift = self.drift = new_drift
            self.i_drift = -1 * self.drift / (1.0 + self.drift)
            return True

        drift_error = new_drift - self.raw_drift
        if abs(drift_error) > self.drift_max_delta:
            # Too far away from the value we expect, discard
            return False

        # move towards the new value
        self.raw_drift += drift_error * self.KP
        self.drift = self.raw_drift - self.KI * self.cumulative_error
        self.i_drift = -1 * self.drift / (1.0 + self.drift)
        return True

    def _update_offset(self, address, base_ts, peer_ts, prediction_error, outlier):
        # insert this into self.ts_base / self.ts_peer / self.var in the right place
        if self.n != 0:
            assert base_ts > self.ts_base[-1]

            # ts_base and ts_peer define a function constructed by linearly
            # interpolating between each pair of values.
            #
            # This function must be monotonically increasing or one of our clocks
            # has effectively gone backwards. If this happens, give up and start
            # again.

            if peer_ts < self.ts_peer[-1]:
                self.ts_base = []
                self.ts_peer = []
                self.var = []
                self.var_sum = 0
                self.updateVars()
                self.cumulative_error = 0
                self.n = 0

                if not self.jumped:
                    self.jumped = 1
                    if self.peer.user.startswith("euerdorf") or self.base.user.startswith("euerdorf"):
                        glogger.warn("{0}: monotonicity broken, reset".format(self))
                    #if self.peer.bad_syncs < 0.1 and self.base.bad_syncs < 0.1:
                    #    glogger.warn("{0}: monotonicity broken, reset".format(self))

                    if self.peer.bad_syncs < 0.1:
                        self.base.incrementJumps()
                    if self.base.bad_syncs < 0.1:
                        self.peer.incrementJumps()

        self.n += 1
        self.ts_base.append(base_ts)
        self.ts_peer.append(peer_ts)

        p_var = prediction_error ** 2
        self.var.append(p_var)
        self.var_sum += p_var
        self.updateVars()

        # do not include outliers in our integral term
        if not outlier:
            self.cumulative_error = max(-50e-6, min(50e-6, self.cumulative_error + prediction_error))  # limit to 50us
            self.outliers = max(0, self.outliers - 1)

        if outlier and not self.jumped:
            self.jumped = 1
            if self.peer.user.startswith("euerdorf") or self.base.user.startswith("euerdorf"):
                glogger.warning("{r}: {a:06X}: step by {e:.1f}us".format(r=self, a=address, e=prediction_error*1e6))
            #if self.peer.bad_syncs < 0.1 and self.base.bad_syncs < 0.1:
            #    glogger.warning("{r}: {a:06X}: step by {e:.1f}us".format(r=self, a=address, e=prediction_error*1e6))
            if self.peer.bad_syncs < 0.1:
                self.base.incrementJumps()
            if self.base.bad_syncs < 0.1:
                self.peer.incrementJumps()

    def predict_peer(self, base_ts):
        """
        Given a time from the base clock, predict the time of the peer clock.
        """

        if self.n == 0:
            return None

        if base_ts < self.ts_base[0] or self.n == 1:
            # extrapolate before first point or if we only have one point
            elapsed = base_ts - self.ts_base[0]
            return (self.ts_peer[0] +
                    elapsed * self.relative_freq +
                    elapsed * self.relative_freq * self.drift)

        if base_ts > self.ts_base[-2]:
            # extrapolate after or before the last point
            elapsed = base_ts - self.ts_base[-1]
            result = (self.ts_peer[-1] +
                    elapsed * self.relative_freq +
                    elapsed * self.relative_freq * self.drift)

            if self.ts_base[-1] - self.ts_base[-2] > 10 * self.base_clock.freq and base_ts > self.ts_base[-1]:
                return result

            elapsed = base_ts - self.ts_base[-2]
            result += (self.ts_peer[-2] +
                    elapsed * self.relative_freq +
                    elapsed * self.relative_freq * self.drift)
            return result / 2

        i = bisect.bisect_left(self.ts_base, base_ts)
        # interpolate between two points
        return (self.ts_peer[i-1] +
                (self.ts_peer[i] - self.ts_peer[i-1]) *
                (base_ts - self.ts_base[i-1]) /
                (self.ts_base[i] - self.ts_base[i-1]))

    def predict_base(self, peer_ts):
        """
        Given a time from the peer clock, predict the time of the base
        clock.
        """

        if self.n == 0:
            return None

        if peer_ts < self.ts_peer[0] or self.n == 1:
            # extrapolate before first point or if we only have one point
            elapsed = peer_ts - self.ts_peer[0]
            return (self.ts_base[0] +
                    elapsed * self.i_relative_freq +
                    elapsed * self.i_relative_freq * self.i_drift)

        if peer_ts > self.ts_peer[-2]:
            # extrapolate after or before the last point
            elapsed = peer_ts - self.ts_peer[-1]
            result = (self.ts_base[-1] +
                    elapsed * self.i_relative_freq +
                    elapsed * self.i_relative_freq * self.i_drift)

            if self.ts_peer[-1] - self.ts_peer[-2] > 10 * self.peer_clock.freq and peer_ts > self.ts_peer[-1]:
                return result

            elapsed = peer_ts - self.ts_peer[-2]
            result += (self.ts_base[-2] +
                    elapsed * self.i_relative_freq +
                    elapsed * self.i_relative_freq * self.i_drift)
            return result / 2

        i = bisect.bisect_left(self.ts_peer, peer_ts)
        # interpolate between two points
        return (self.ts_base[i-1] +
                (self.ts_base[i] - self.ts_base[i-1]) *
                (peer_ts - self.ts_peer[i-1]) /
                (self.ts_peer[i] - self.ts_peer[i-1]))

    def __str__(self):
        return self.base.user + ':' + self.peer.user
