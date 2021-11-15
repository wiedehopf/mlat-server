# -*- mode: python; indent-tabs-mode: nil -*-

# Part of mlat-server: a Mode S multilateration server
# Copyright (C) 2015  Oliver Jowett <oliver@mutability.co.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Utility functions to convert between coordinate systems and calculate distances.
"""

from libc.math cimport sqrt, sin, cos, atan2, acos

import math

from cpython cimport array
import array

# degrees to radians
cdef double DTOR = math.pi / 180.0
# radians to degrees
cdef double RTOD = 180.0 / math.pi

# WGS84 ellipsoid Earth parameters
cdef double WGS84_A = 6378137.0
cdef double WGS84_F = 1.0/298.257223563
cdef double WGS84_B = WGS84_A * (1 - WGS84_F)
cdef double WGS84_ECC_SQ = 1 - WGS84_B * WGS84_B / (WGS84_A * WGS84_A)
cdef double WGS84_ECC = sqrt(WGS84_ECC_SQ)

# Average radius for a spherical Earth
SPHERICAL_R = 6371e3

# Some derived values
_wgs84_ep = math.sqrt((WGS84_A**2 - WGS84_B**2) / WGS84_B**2)
_wgs84_ep2_b = _wgs84_ep**2 * WGS84_B
_wgs84_e2_a = WGS84_ECC_SQ * WGS84_A


def llh2ecef(llh):
    """Converts from WGS84 lat/lon/height to ellipsoid-earth ECEF"""

    cdef double lat = llh[0] * DTOR
    cdef double lng = llh[1] * DTOR
    cdef double alt = llh[2]

    cdef double slat = sin(lat)
    cdef double slng = sin(lng)
    cdef double clat = cos(lat)
    cdef double clng = cos(lng)

    cdef double d = sqrt(1 - (slat * slat * WGS84_ECC_SQ))
    cdef double rn = WGS84_A / d

    cdef double x = (rn + alt) * clat * clng
    cdef double y = (rn + alt) * clat * slng
    cdef double z = (rn * (1 - WGS84_ECC_SQ) + alt) * slat

    return array.array('d', (x, y, z))


def ecef2llh(ecef):
    "Converts from ECEF to WGS84 lat/lon/height"

    cdef double x, y, z
    x, y, z = ecef

    cdef double lon = atan2(y, x)

    cdef double p = sqrt(x**2 + y**2)
    cdef double th = atan2(WGS84_A * z, WGS84_B * p)
    cdef double lat = atan2(z + _wgs84_ep2_b * sin(th)**3,
                     p - _wgs84_e2_a * cos(th)**3)

    cdef double N = WGS84_A / sqrt(1 - WGS84_ECC_SQ * sin(lat)**2)
    cdef double alt = p / cos(lat) - N

    return array.array('d', (lat * RTOD, lon * RTOD, alt))


def greatcircle(p0, p1):
    """Returns a great-circle distance in metres between two LLH points,
    _assuming spherical earth_ and _ignoring altitude_. Don't use this if you
    need a distance accurate to better than 1%."""

    cdef double lat0 = p0[0] * DTOR
    cdef double lon0 = p0[1] * DTOR
    cdef double lat1 = p1[0] * DTOR
    cdef double lon1 = p1[1] * DTOR
    return SPHERICAL_R * acos(
        sin(lat0) * sin(lat1) +
        cos(lat0) * cos(lat1) * cos(abs(lon0 - lon1)))


# direct implementation here turns out to be _much_ faster (10-20x) compared to
# scipy.spatial.distance.euclidean or numpy-based approaches
def ecef_distance(p0, p1):
    """Returns the straight-line distance in metres between two ECEF points."""
    return sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)
