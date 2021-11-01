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

import asyncio
import logging
import logging.handlers
import time
import math
import functools
import socket
import numpy

from mlat import constants, geodesy
from mlat.server import util, net

"""
Various output methods for multilateration results.
"""


def format_time(timestamp):
    return time.strftime("%H:%M:%S", time.gmtime(timestamp)) + ".{0:03.0f}".format(math.modf(timestamp)[0] * 1000)


def format_date(timestamp):
    return time.strftime("%Y/%m/%d", time.gmtime(timestamp))


def csv_quote(s):
    if s is None:
        return ''
    if s.find('\n') == -1 and s.find('"') == -1 and s.find(',') == -1:
        return s
    else:
        return '"' + s.replace('"', '""') + '"'


class LocalCSVWriter(object):
    """Writes multilateration results to a local CSV file"""

    TEMPLATE = '{t:.3f},{address:06X},{callsign},{squawk},{lat:.4f},{lon:.4f},{alt},{err:.0f},{n},{d},{receivers},{dof},{vrate}'  # noqa
    KTEMPLATE = '{t:.3f},{address:06X},{callsign},{squawk},{lat:.4f},{lon:.4f},{alt},{err:.0f},{n},{d},{receivers},{dof},{klat:.4f},{klon:.4f},{kalt:.0f},{kheading:.0f},{kspeed:.0f},{vrate},{kerr:.0f}'  # noqa

    def __init__(self, coordinator, filename):
        self.logger = logging.getLogger("csv")
        self.coordinator = coordinator
        self.filename = filename
        self.coordinator.add_output_handler(self.write_result)

        self.pos_logger = logging.getLogger("positions")
        self.pos_logger.setLevel(logging.DEBUG)

        self.pos_handler = logging.handlers.RotatingFileHandler(
                self.filename, maxBytes=(1*1024*1024), backupCount=10)
        self.pos_logger.addHandler(self.pos_handler)

    def start(self):
        return util.completed_future

    def close(self):
        self.coordinator.remove_output_handler(self.write_result)

    def wait_closed(self):
        return util.completed_future

    def write_result(self, receive_timestamp, address, ecef, ecef_cov, receivers, distinct, dof, kalman_state, error):
        try:
            lat, lon, alt = geodesy.ecef2llh(ecef)

            ac = self.coordinator.tracker.aircraft[address]
            callsign = ac.callsign
            squawk = ac.squawk

            if ecef_cov is None:
                err_est = -1
            else:
                var_est = numpy.sum(numpy.diagonal(ecef_cov))
                if var_est >= 0:
                    err_est = math.sqrt(var_est)
                else:
                    err_est = -1
            # never use MLAT calculated altitude, most of the time it's so inaccurate that it's useless
            # better have the altitude be undefined
            if ac.last_altitude_time and receive_timestamp - ac.last_altitude_time < 5:
                # ft
                altitude = str(ac.altitude)
            else:
                altitude = ''
            if ac.vrate_time and receive_timestamp - ac.vrate_time < 5:
                # fpm
                vrate = str(ac.vrate)
            else:
                vrate = ''

            if kalman_state.valid and kalman_state.last_update >= receive_timestamp:
                line = self.KTEMPLATE.format(
                    t=receive_timestamp,
                    address=address,
                    callsign=csv_quote(callsign),
                    squawk=csv_quote(squawk),
                    lat=lat,
                    lon=lon,
                    alt=altitude,
                    err=err_est,
                    n=len(receivers),
                    d=distinct,
                    dof=dof,
                    receivers=csv_quote(','.join([receiver.user for receiver in receivers])),
                    klat=kalman_state.position_llh[0],
                    klon=kalman_state.position_llh[1],
                    kalt=kalman_state.position_llh[2] * constants.MTOF,
                    kheading=kalman_state.heading,
                    kspeed=kalman_state.ground_speed * constants.MS_TO_KTS,
                    vrate=vrate,
                    kerr=kalman_state.position_error)
            else:
                line = self.TEMPLATE.format(
                    t=receive_timestamp,
                    address=address,
                    callsign=csv_quote(callsign),
                    squawk=csv_quote(squawk),
                    lat=lat,
                    lon=lon,
                    alt=altitude,
                    err=err_est,
                    n=len(receivers),
                    d=distinct,
                    dof=dof,
                    receivers=csv_quote(','.join([receiver.user for receiver in receivers])),
                    vrate=vrate)


            self.pos_logger.debug(line)

        except Exception:
            self.logger.exception("Failed to write result")
            # swallow the exception so we don't affect our caller


class BasestationClient(object):
    """Writes results in Basestation port-30003 format to network clients."""

    TEMPLATE = 'MSG,{mtype},1,1,{addr:06X},1,{rcv_date},{rcv_time},{now_date},{now_time},{callsign},{altitude},{speed},{heading},{lat},{lon},{vrate},{squawk},{fs},{emerg},{ident},{aog}\n'  # noqa

    def __init__(self, reader, writer, *, coordinator, use_kalman_data, heartbeat_interval=30.0):
        peer = writer.get_extra_info('peername')
        self.host = peer[0]
        self.port = peer[1]
        self.logger = util.TaggingLogger(logging.getLogger("basestation"),
                                         {'tag': '{host}:{port}'.format(host=self.host,
                                                                        port=self.port)})
        self.reader = reader
        self.writer = writer
        self.coordinator = coordinator
        self.use_kalman_data = use_kalman_data
        self.heartbeat_interval = heartbeat_interval
        self.last_output = time.time()
        self.heartbeat_task = asyncio.ensure_future(self.send_heartbeats())
        self.reader_task = asyncio.ensure_future(self.read_until_eof())

        self.logger.warning("Connection established")
        self.coordinator.add_output_handler(self.write_result)

    def close(self):
        if not self.writer:
            return  # already closed

        self.logger.info("Connection lost")
        self.coordinator.remove_output_handler(self.write_result)
        self.heartbeat_task.cancel()
        self.writer.close()
        self.writer = None

    @asyncio.coroutine
    def wait_closed(self):
        yield from util.safe_wait([self.heartbeat_task, self.reader_task])

    @asyncio.coroutine
    def read_until_eof(self):
        try:
            while True:
                r = yield from self.reader.read(1024)
                if len(r) == 0:
                    self.logger.info("Client EOF")
                    # EOF
                    self.close()
                    return
        except socket.error:
            self.close()
            return

    @asyncio.coroutine
    def send_heartbeats(self):
        try:
            while True:
                now = time.time()
                delay = self.last_output + self.heartbeat_interval - now
                if delay > 0.1:
                    yield from asyncio.sleep(delay)
                    continue

                self.writer.write(b'\n')
                self.last_output = now

        except socket.error:
            self.close()
            return

    def write_result(self, receive_timestamp, address, ecef, ecef_cov, receivers, distinct, dof, kalman_data, error):
        try:
            ac = self.coordinator.tracker.aircraft[address]

            speed = ''
            heading = ''
            vrate = ''
            altitude = ''

            if self.use_kalman_data:
                if not kalman_data.valid and dof < 1:
                    if receive_timestamp - ac.last_crappy_output > 60:
                        # ignore the first crappy output
                        ac.last_crappy_output = receive_timestamp - 1
                        return
                    if receive_timestamp - ac.last_crappy_output > 30 or receive_timestamp == ac.last_crappy_output:
                        ac.last_crappy_output = receive_timestamp
                    else:
                        return

                if not kalman_data.valid or kalman_data.last_update < receive_timestamp:
                    lat, lon, alt = geodesy.ecef2llh(ecef)
                else:
                    # always use non kalman position, only speed, heading, vertical speed are used from kalman
                    # lat, lon, alt = kalman_data.position_llh
                    lat, lon, alt = geodesy.ecef2llh(ecef)
                    speed = int(round(kalman_data.ground_speed * constants.MS_TO_KTS))
                    heading = int(round(kalman_data.heading))
                    vrate = int(round(kalman_data.vertical_speed * constants.MS_TO_FPM))

            else:
                lat, lon, alt = geodesy.ecef2llh(ecef)

            callsign = ac.callsign
            squawk = ac.squawk
            send_timestamp = time.time()

            # never use MLAT calculated altitude, most of the time it's so inaccurate that it's useless
            # better have the altitude be undefined
            if ac.last_altitude_time and receive_timestamp - ac.last_altitude_time < 5:
                # ft
                altitude = ac.altitude
            else:
                altitude = ''
            if ac.vrate_time and receive_timestamp - ac.vrate_time < 5:
                # fpm
                vrate = ac.vrate
            else:
                vrate = ''

            line = self.TEMPLATE.format(mtype=3,
                                        addr=address,
                                        rcv_date=format_date(receive_timestamp),
                                        rcv_time=format_time(receive_timestamp),
                                        now_date=format_date(send_timestamp),
                                        now_time=format_time(send_timestamp),
                                        callsign=csv_quote(callsign),
                                        squawk=csv_quote(squawk),
                                        lat=round(lat, 6),
                                        lon=round(lon, 6),
                                        altitude=altitude,
                                        speed=speed,
                                        heading=heading,
                                        vrate=vrate,
                                        fs=len(receivers),
                                        emerg=(error if error is not None else ''),
                                        ident='',
                                        aog='')
            self.writer.write(line.encode('ascii'))
            self.last_output = time.time()

        except Exception:
            self.logger.exception("Failed to write result")
            # swallow the exception so we don't affect our caller


def make_basestation_listener(host, port, coordinator, use_kalman_data):
    factory = functools.partial(BasestationClient,
                                coordinator=coordinator,
                                use_kalman_data=use_kalman_data)
    return net.MonitoringListener(host, port, factory,
                                  logger=logging.getLogger('basestation'),
                                  description='Basestation output listener')


def make_basestation_connector(host, port, coordinator, use_kalman_data):
    factory = functools.partial(BasestationClient,
                                coordinator=coordinator,
                                use_kalman_data=use_kalman_data)
    return net.MonitoringConnector(host, port, 30.0, factory)
