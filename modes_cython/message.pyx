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
Decoder for Mode S responses and ADS-B extended squitter messages.
"""

__all__ = ('ESType', 'decode', 'decode_cpr', 'DF0', 'DF4', 'DF5', 'DF11', 'DF16',
           'DF17', 'DF18', 'DF20', 'DF21', 'ExtendedSquitter', 'CommB')

from enum import Enum
import math
import bisect

ais_charset = " ABCDEFGHIJKLMNOPQRSTUVWXYZ????? ???????????????0123456789??????"

class ModeSMessage:
    """
    A decoded Mode S message.

    All subclasses have the following fields present, though some may be
    set to None:

      DF: downlink format
      address: ICAO address of transmitting aircraft. For some message types
        this is derived from the CRC field and may be unreliable.
      altitude: decoded altitude in feet, or None if not present / not available
      callsign: decoded callsign, or None if not present
      squawk: decoded squawk, or None if not present
      crc_ok: True if the CRC is OK. False if it is bad. None if the correctness
        of the CRC cannot be checked (e.g. the messages uses AP or PI)
    """


class DF0(ModeSMessage):
    """
    DF0 (Short air-air surveillance / ACAS) message.

    Fields: DF, VS, CC, SL, RI, AC, altitude, address
    """

    def __init__(self, frombuf):
        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.VS = (frombuf[0] & 0x04) >> 2  # 1 bit
        self.CC = (frombuf[0] & 0x02) >> 1  # 1 bit
        # 1 bit pad
        self.SL = (frombuf[1] & 0xe0) >> 5  # 3 bits
        # 2 bits pad
        self.RI = ((frombuf[1] & 0x03) << 1) | ((frombuf[2] & 0x80) >> 7)  # 4 bits
        # 2 bits pad
        self.AC = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        # 24 bits A/P

        self.squawk = self.callsign = None
        self.altitude = decode_ac13(self.AC)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class DF4(ModeSMessage):
    """
    DF4 (Surveillance, altitude reply) message.

    Fields: DF, FS, DR, UM, AC, altitude, address
    """

    def __init__(self, frombuf):
        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.FS = (frombuf[0] & 0x07)       # 3 bits
        self.DR = (frombuf[1] & 0xf8) >> 3  # 5 bits
        self.UM = ((frombuf[1] & 0x07) << 3) | ((frombuf[2] & 0xe0) >> 5)  # 6 bits
        self.AC = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        # 24 bits A/P

        self.squawk = self.callsign = None
        self.altitude = decode_ac13(self.AC)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class DF5(ModeSMessage):
    """
    DF5 (Surveillance, identity reply) message.

    Fields: DF, FS, DR, UM, ID, squawk, address
    """

    def __init__(self, frombuf):
        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.FS = (frombuf[0] & 0x07)       # 3 bits
        self.DR = (frombuf[1] & 0xf8) >> 3  # 5 bits
        self.UM = ((frombuf[1] & 0x07) << 3) | ((frombuf[2] & 0xe0) >> 5)  # 6 bits
        self.ID = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        # 24 bits A/P

        self.altitude = self.callsign = None
        self.squawk = decode_id13(self.ID)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class DF11(ModeSMessage):
    """
    DF11 (All-call reply) message.

    Fields: DF, CA, AA, address, crc_ok
    """

    def __init__(self, frombuf):
        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.CA = (frombuf[0] & 0x07)       # 3 bits
        self.AA = (frombuf[1] << 16) | (frombuf[2] << 8) | frombuf[3]  # 24 bits
        # 24 bits P/I

        self.squawk = self.callsign = self.altitude = None

        r = crc_residual(frombuf)
        if r == 0:
            self.crc_ok = True
        elif (r & ~0x7f) == 0:
            self.crc_ok = None
        else:
            self.crc_ok = False
        self.address = self.AA


class DF16(ModeSMessage):
    """
    DF16 (Long air-air surveillance / ACAS) message.

    Fields: DF, VS, SL, RI, AC, altitude, address
    """

    def __init__(self, frombuf):
        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.VS = (frombuf[0] & 0x04) >> 2  # 1 bit
        # 2 bits pad
        self.SL = (frombuf[1] & 0xe0) >> 5  # 3 bits
        # 2 bits pad
        self.RI = ((frombuf[1] & 0x03) << 1) | ((frombuf[2] & 0x80) >> 7)  # 4 bits
        # 2 bits pad
        self.AC = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        self.MV = frombuf[4:11]  # 56 bits
        # 24 bits A/P

        self.squawk = self.callsign = None
        self.altitude = decode_ac13(self.AC)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class CommB(ModeSMessage):
    """A message containing a Comm-B reply.

    Fields: MB, callsign
    """

    def __init__(self, frombuf):
        self.MB = frombuf[4:11]  # 56 bits

        if frombuf[4] != 0x20:
            self.callsign = None
        else:
            callsign = (
                ais_charset[(frombuf[5] & 0xfc) >> 2] +
                ais_charset[((frombuf[5] & 0x03) << 4) | ((frombuf[6] & 0xf0) >> 4)] +
                ais_charset[((frombuf[6] & 0x0f) << 2) | ((frombuf[7] & 0xc0) >> 6)] +
                ais_charset[frombuf[7] & 0x3f] +
                ais_charset[(frombuf[8] & 0xfc) >> 2] +
                ais_charset[((frombuf[8] & 0x03) << 4) | ((frombuf[9] & 0xf0) >> 4)] +
                ais_charset[((frombuf[9] & 0x0f) << 2) | ((frombuf[10] & 0xc0) >> 6)] +
                ais_charset[frombuf[10] & 0x3f]
            )

            if callsign != '        ' and callsign.find('?') == -1:
                self.callsign = callsign
            else:
                self.callsign = None


class DF20(CommB):
    """
    DF20 (Comm-B, altitude reply) message.

    Fields: DF, FS, DR, UM, AC, altitude, address, MB, callsign
    """

    def __init__(self, frombuf):
        CommB.__init__(self, frombuf)

        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.FS = (frombuf[0] & 0x07)       # 3 bits
        self.DR = (frombuf[1] & 0xf8) >> 3  # 5 bits
        self.UM = ((frombuf[1] & 0x07) << 3) | ((frombuf[2] & 0xe0) >> 5)  # 6 bits
        self.AC = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        # 56 bits MB
        # 24 bits A/P

        self.squawk = None
        self.altitude = decode_ac13(self.AC)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class DF21(CommB):
    """
    DF21 (Comm-B, identity reply) message.

    Fields: DF, FS, DR, UM, ID, squawk, address, MB, callsign
    """

    def __init__(self, frombuf):
        CommB.__init__(self, frombuf)

        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.FS = (frombuf[0] & 0x07)       # 3 bits
        self.DR = (frombuf[1] & 0xf8) >> 3  # 5 bits
        self.UM = ((frombuf[1] & 0x07) << 3) | ((frombuf[2] & 0xe0) >> 5)  # 6 bits
        self.ID = ((frombuf[2] & 0x1f) << 8) | frombuf[3]  # 13 bits
        # 56 bits MB
        # 24 bits A/P

        self.altitude = None
        self.squawk = decode_id13(self.ID)
        self.crc_ok = None
        self.address = crc_residual(frombuf)


class ESType(Enum):
    """Identifies the type of an Extended Squitter message."""
    id_and_category = 1
    airborne_position = 2
    surface_position = 3
    airborne_velocity = 4
    other = 5

es_types = {
    0: (ESType.airborne_position, 0),
    1: (ESType.id_and_category, None),
    2: (ESType.id_and_category, None),
    3: (ESType.id_and_category, None),
    4: (ESType.id_and_category, None),
    5: (ESType.surface_position, 9),
    6: (ESType.surface_position, 8),
    7: (ESType.surface_position, 7),
    8: (ESType.surface_position, 6),
    9: (ESType.airborne_position, 9),
    10: (ESType.airborne_position, 8),
    11: (ESType.airborne_position, 7),
    12: (ESType.airborne_position, 6),
    13: (ESType.airborne_position, 5),
    14: (ESType.airborne_position, 4),
    15: (ESType.airborne_position, 3),
    16: (ESType.airborne_position, 2),
    17: (ESType.airborne_position, 1),
    18: (ESType.airborne_position, 0),
    19: (ESType.airborne_velocity, None),
    20: (ESType.airborne_position, 9),
    21: (ESType.airborne_position, 8),
    22: (ESType.airborne_position, 0)
}


class ExtendedSquitter(ModeSMessage):
    """A message that carries an Extended Squitter message.

    Fields: estype, nuc

    For airborne positions: SS, SAF, AC12, T, F, LAN, LON, altitude
    For id and category: CATEGORY, callsign
    """

    def __init__(self, frombuf):
        metype = (frombuf[4] & 0xf8) >> 3
        self.estype, self.nuc = es_types.get(metype, (ESType.other, None))

        if self.estype is ESType.airborne_position:
            self.SS = (frombuf[4] & 0x06) >> 1
            self.SAF = frombuf[4] & 0x01
            self.AC12 = (frombuf[5] << 4) | ((frombuf[6] & 0xf0) >> 4)
            self.T = (frombuf[6] & 0x08) >> 3
            self.F = (frombuf[6] & 0x04) >> 2
            self.LAT = (((frombuf[6] & 0x03) << 15) |
                        (frombuf[7] << 7) |
                        ((frombuf[8] & 0xfe) >> 1))
            self.LON = (((frombuf[8] & 0x01) << 16) |
                        (frombuf[9] << 8) |
                        frombuf[10])
            self.altitude = decode_ac12(self.AC12)
            self.callsign = None

        elif self.estype is ESType.id_and_category:
            self.CATEGORY = frombuf[4] & 0x07
            self.altitude = None
            self.callsign = (
                ais_charset[(frombuf[5] & 0xfc) >> 2] +
                ais_charset[((frombuf[5] & 0x03) << 4) | ((frombuf[6] & 0xf0) >> 4)] +
                ais_charset[((frombuf[6] & 0x0f) << 2) | ((frombuf[7] & 0xc0) >> 6)] +
                ais_charset[frombuf[7] & 0x3f] +
                ais_charset[(frombuf[8] & 0xfc) >> 2] +
                ais_charset[((frombuf[8] & 0x03) << 4) | ((frombuf[9] & 0xf0) >> 4)] +
                ais_charset[((frombuf[9] & 0x0f) << 2) | ((frombuf[10] & 0xc0) >> 6)] +
                ais_charset[frombuf[10] & 0x3f]
            )

        else:
            self.altitude = None
            self.callsign = None


class DF17(ExtendedSquitter):
    """DF17 (Extended Squitter) message.

    Fields: DF, CA, AA, address, crc_ok; plus those of ExtendedSquitter.
    """

    def __init__(self, frombuf):
        ExtendedSquitter.__init__(self, frombuf)

        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.CA = (frombuf[0] & 0x07)       # 3 bits
        self.AA = (frombuf[1] << 16) | (frombuf[2] << 8) | frombuf[3]  # 24 bits
        # 56 bits ME
        # 24 bits CRC

        self.squawk = None
        self.crc_ok = (crc_residual(frombuf) == 0)
        self.address = self.AA


class DF18(ExtendedSquitter):
    """DF18 (Extended Squitter / Non-Transponder) message.

    Fields: DF, CF, AA, address, crc_ok; plus those of ExtendedSquitter.
    """

    def __init__(self, frombuf):
        ExtendedSquitter.__init__(self, frombuf)

        self.DF = (frombuf[0] & 0xf8) >> 3  # 5 bits
        self.CF = (frombuf[0] & 0x07)       # 3 bits
        self.AA = (frombuf[1] << 16) | (frombuf[2] << 8) | frombuf[3]  # 24 bits
        # 56 bits ME
        # 24 bits CRC

        self.squawk = None
        self.crc_ok = (crc_residual(frombuf) == 0)
        self.address = self.AA


message_types = {
    0: DF0,
    4: DF4,
    5: DF5,
    11: DF11,
    16: DF16,
    17: DF17,
    18: DF18,
    20: DF20,
    21: DF21
}


def decode(frombuf):
    """
    Decode a Mode S message.

      frombuf: a 7-byte or 14-byte message containing the encoded Mode S message

    Returns a suitable message object, or None if the message type is not
    handled.
    """

    df = (frombuf[0] & 0xf8) >> 3
    try:
        return message_types[df](frombuf)
    except KeyError:
        return None
    except IndexError:
        return None

"""
Decoders for the 12- and 13-bit altitude encodings used in Mode S responses
and ADS-B extended squitter messages.
"""

def _decode_ac13(ac13):
    if ac13 is None or ac13 == 0:    # no data
        return None
    if ac13 & 0x0040:                # M bit set
        return None
    if ac13 & 0x0010:                # Q bit set
        n = ((ac13 & 0x1f80) >> 2) | ((ac13 & 0x0020) >> 1) | (ac13 & 0x000f)
        return n * 25 - 1000

    # convert from Gillham code
    if not (ac13 & 0x1500):
        return None  # illegal C bits

    h = 0
    if ac13 & 0x1000:
        h ^= 7  # C1
    if ac13 & 0x0400:
        h ^= 3  # C2
    if ac13 & 0x0100:
        h ^= 1  # C4
    if h & 5:
        h ^= 5
    if h > 5:
        return None  # illegal C bits

    f = 0
    if ac13 & 0x0010:
        f ^= 0x1ff  # D1
    if ac13 & 0x0004:
        f ^= 0x0ff  # D2
    if ac13 & 0x0001:
        f ^= 0x07f  # D4
    if ac13 & 0x0800:
        f ^= 0x03f  # A1
    if ac13 & 0x0200:
        f ^= 0x01f  # A2
    if ac13 & 0x0080:
        f ^= 0x00f  # A4
    if ac13 & 0x0020:
        f ^= 0x007  # B1
    if ac13 & 0x0008:
        f ^= 0x003  # B2
    if ac13 & 0x0002:
        f ^= 0x001  # B4

    if f & 1:
        h = (6 - h)

    a = 500 * f + 100 * h - 1300
    if a < -1200:
        return None  # illegal value

    return a


def decode_ac13(ac13):
    """Decodes a Mode S 13-bit altitude field.

    The expected ordering is as specified in §3.1.2.6.5.4 of Annex 10:

      C1, A1, C2, A2, C4, A4, (M), B1, (Q), B2, D2, B4, D4

    Returns signed altitude in feet, or None if not decodable.
    """

    if ac13 is None:
        return None
    return _alt_table[ac13 & 0x1fff]


def decode_ac12(ac12):
    """Decode a 12-bit AC altitude field from an extended squitter.

    The expected ordering is as specified in Doc 9871 Table A-2-5:

     the altitude code (AC) as specified in §3.1.2.6.5.4 of Annex 10,
     but with the M-bit removed

    Returns signed altitude in feet, or None if not a valid altitude."""

    if ac12 is None:
        return None
    return _alt_table[((ac12 & 0x0fc0) << 1) | (ac12 & 0x003f)]


# precompute the lookup table
_alt_table = [_decode_ac13(i) for i in range(2**13)]

"""
Decoder for 12-bit squawk (identity) fields contained in some Mode S messages.
"""

def _make_upper_table():
    ut = []
    for i in range(64):
        v = 0
        id13 = i << 7
        if id13 & 0x1000:
            v |= 0x0010  # C1
        if id13 & 0x0800:
            v |= 0x1000  # A1
        if id13 & 0x0400:
            v |= 0x0020  # C2
        if id13 & 0x0200:
            v |= 0x2000  # A2
        if id13 & 0x0100:
            v |= 0x0040  # C4
        if id13 & 0x0080:
            v |= 0x4000  # A4
        ut.append(v)
    return ut


def _make_lower_table():
    lt = []
    for id13 in range(64):
        v = 0
        # 0040 unused (M/X)
        if id13 & 0x0020:
            v |= 0x0100  # B1
        if id13 & 0x0010:
            v |= 0x0001  # D1/Q
        if id13 & 0x0008:
            v |= 0x0200  # B2
        if id13 & 0x0004:
            v |= 0x0002  # D2
        if id13 & 0x0002:
            v |= 0x0400  # B4
        if id13 & 0x0001:
            v |= 0x0004  # D4
        lt.append(v)

    return lt


def decode_id13(id13):
    """Decode a 13-bit Mode A squawk.

    The expected ordering is that from Annex 10 vol 4 3.1.2.6.7.1:

      C1, A1, C2, A2, C4, A4, ZERO, B1, D1, B2, D2, B4, D4

    Returns the squawk as a 4-character string."""

    return '{0:04x}'.format(_id13_lt[id13 & 63] | _id13_ut[id13 >> 7])

_id13_lt = _make_lower_table()
_id13_ut = _make_upper_table()

"""
Calculates the 24-bit CRC used in Mode S messages.
"""

def crc_residual(payload):
    """Computes the 24-bit Mode S CRC residual for a message.

    The CRC residual is the CRC computed across the first 4 or 11 bytes,
    XOR-ed with the CRC value stored in the final 3 bytes.

    For a message using Address/Parity, the expected residual is the
    transmitter's address.

    For a message using Parity/Interrogator, the expected residual is
    the interrogator ID.

    For an extended squitter message or a DF11 acquisition squitter, the
    expected residual is zero.

    Errors in the message or in the CRC value itself will appear as errors
    in the residual value.
    """

    t = _crc_table
    rem = t[payload[0]]
    for b in payload[1:-3]:
        rem = ((rem & 0xFFFF) << 8) ^ t[b ^ (rem >> 16)]

    rem = rem ^ (payload[-3] << 16) ^ (payload[-2] << 8) ^ (payload[-1])
    return rem


def _make_crc_table():
    # precompute the CRC table
    t = []

    poly = 0xfff409
    for i in range(256):
        c = i << 16
        for j in range(8):
            if c & 0x800000:
                c = (c << 1) ^ poly
            else:
                c = (c << 1)

        t.append(c & 0xffffff)

    return t

_crc_table = _make_crc_table()

"""
Decoder for the Compact Position Reporting (CPR) position encoding used in
ADS-B extended squitter messages.
"""

nl_table = (
    (10.47047130, 59),
    (14.82817437, 58),
    (18.18626357, 57),
    (21.02939493, 56),
    (23.54504487, 55),
    (25.82924707, 54),
    (27.93898710, 53),
    (29.91135686, 52),
    (31.77209708, 51),
    (33.53993436, 50),
    (35.22899598, 49),
    (36.85025108, 48),
    (38.41241892, 47),
    (39.92256684, 46),
    (41.38651832, 45),
    (42.80914012, 44),
    (44.19454951, 43),
    (45.54626723, 42),
    (46.86733252, 41),
    (48.16039128, 40),
    (49.42776439, 39),
    (50.67150166, 38),
    (51.89342469, 37),
    (53.09516153, 36),
    (54.27817472, 35),
    (55.44378444, 34),
    (56.59318756, 33),
    (57.72747354, 32),
    (58.84763776, 31),
    (59.95459277, 30),
    (61.04917774, 29),
    (62.13216659, 28),
    (63.20427479, 27),
    (64.26616523, 26),
    (65.31845310, 25),
    (66.36171008, 24),
    (67.39646774, 23),
    (68.42322022, 22),
    (69.44242631, 21),
    (70.45451075, 20),
    (71.45986473, 19),
    (72.45884545, 18),
    (73.45177442, 17),
    (74.43893416, 16),
    (75.42056257, 15),
    (76.39684391, 14),
    (77.36789461, 13),
    (78.33374083, 12),
    (79.29428225, 11),
    (80.24923213, 10),
    (81.19801349, 9),
    (82.13956981, 8),
    (83.07199445, 7),
    (83.99173563, 6),
    (84.89166191, 5),
    (85.75541621, 4),
    (86.53536998, 3),
    (87.00000000, 2),
    (90.00000000, 1)
)

nl_lats = [x[0] for x in nl_table]
nl_vals = [x[1] for x in nl_table]


def NL(lat):
    if lat < 0:
        lat = -lat

    nl = nl_vals[bisect.bisect_left(nl_lats, lat)]
    return nl


def MOD(a, b):
    r = a % b
    if r < 0:
        r += b
    return r


def decode_cpr(latE, lonE, latO, lonO):
    """Perform globally unambiguous position decoding for a pair of
    airborne CPR messages.

    latE, lonE: the raw latitude and longitude values of the even message
    latO, lonO: the raw latitude and longitude values of the odd message

    Return a tuple of (even latitude, even longitude, odd latitude, odd longitude)

    Raises ValueError if the messages do not produce a useful position."""

    # Compute the Latitude Index "j"
    j = math.floor(((59 * latE - 60 * latO) / 131072.0) + 0.5)
    rlatE = (360.0 / 60.0) * (MOD(j, 60) + latE / 131072.0)
    rlatO = (360.0 / 59.0) * (MOD(j, 59) + latO / 131072.0)

    # adjust for southern hemisphere values, which are in the range (270,360)
    if rlatE >= 270:
        rlatE -= 360
    if rlatO >= 270:
        rlatO -= 360

    # Check to see that the latitude is in range: -90 .. +90
    if rlatE < -90 or rlatE > 90 or rlatO < -90 or rlatO > 90:
        raise ValueError('latitude out of range')

    # Find latitude zone, abort if the two positions are not in the same zone
    nl = NL(rlatE)
    if nl != NL(rlatO):
        raise ValueError('messages lie in different latitude zones')

    # Compute n(i)
    nE = nl
    nO = max(1, nl - 1)

    # Compute the Longitude Index "m"
    m = math.floor((((lonE * (nl - 1)) - (lonO * nl)) / 131072.0) + 0.5)

    # Compute global longitudes
    rlonE = (360.0 / nE) * (MOD(m, nE) + lonE / 131072.0)
    rlonO = (360.0 / nO) * (MOD(m, nO) + lonO / 131072.0)

    # Renormalize to -180 .. +180
    rlonE -= math.floor((rlonE + 180) / 360) * 360
    rlonO -= math.floor((rlonO + 180) / 360) * 360

    return (rlatE, rlonE, rlatO, rlonO)
