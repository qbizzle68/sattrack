import json
from enum import Enum
from functools import total_ordering
from math import sin, cos, pi, radians, tan, asin, sqrt, degrees, atan2, acos, atan

from pyevspace import EVector

from sattrack.exceptions import SunRiseSetException
from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.structures.coordinates import GeoPosition
from sattrack.topos import toTopocentric
from sattrack.util.constants import AU, TWOPI
from sattrack.util.conversions import atan3

DELTAT = 72.6  # 142.8 #157.2  # seconds from https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm
_EMPTY_GEO = GeoPosition(0, 0, 0)


@total_ordering
class TwilightType(Enum):
    """
    Enumerated value to distinguish between the various sun angles relative to a horizon.
    The classifications for sun angle to twilight type is:
    Sun < -18 degrees --> Night time
    Sun < -12 degrees --> Astronomical twilight
    Sun < -6 degrees --> Nautical twilight
    Sun < -50 arc-minutes --> Civil twilight
    Otherwise --> Day time
    Most object can't be seen until their position is in civil twilight.
    """

    Day = 0
    Civil = 1
    Nautical = 2
    Astronomical = 3
    Night = 4

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented


def getSunPosition(time: JulianDate) -> EVector:
    """
    Very simple algorithm for computing the position of the Sun in a geocentric equitorial reference frame.

    Args:
        time: Time to find the Sun's position.

    Returns:
        The Sun's position vector in kilometers.
    """

    # #   time since noon TT on Jan 1, 2000
    # n = time.difference(J2000)
    # #   mean longitude of the Sun
    # L = 4.89495042 + 0.0172027924 * n
    # #   mean anomaly of the Sun
    # g = 6.240040768 + 0.0172019703 * n
    # #   ecliptic longitude of the Sun
    # lmbda = L + 0.0334230552 * sin(g) + 0.00034906585 * sin(2 * g)
    # #   distance of the Sun from the Earth
    # R = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2 * g)
    # #   obliquity of the ecliptic (axial tilt)
    # e = 0.4090877234 - 6.981317008e-9 * n
    #
    # return EVector(
    #     R * AU * cos(lmbda),
    #     R * AU * cos(e) * sin(lmbda),
    #     R * AU * sin(e) * sin(lmbda)
    # )
    spc = SunPositionController2(time, _EMPTY_GEO)
    return spc.getSunPosition()


def getSunCoordinates(time: JulianDate) -> EVector:
    spc = SunPositionController2(time, _EMPTY_GEO)
    return spc.getSunCoordinates()


def getSunRiseSetTimes(time: JulianDate, geo: GeoPosition):
    """Returns the solar rise and set times for the day given."""
    spc = SunPositionController2(time, geo)
    vals = spc.getAngleTimes()
    return vals[0], vals[1]

def getSunRiseSetTimes2(time: JulianDate, geo: GeoPosition):
    """Returns the solar rise and set times, but with special rules. If the time is after sunrise, the previous sunrise
    is returned, otherwise the next sunrise will be returned. Therefore, if the time is during the day the returned
    times are that day's sunrise and sunset, and if the time is during the night the returned times are the next day's
    sunrise and sunset. The next sunset is always returned."""
    spc = SunPositionController2(time, geo)
    riseDayOf, setDayOf, transitDayOf = spc.getAngleTimes()
    if time.value() <= setDayOf.value():
        return riseDayOf, setDayOf
    else:
        time_p1 = JulianDate.fromNumber(time.number(), 0.51 - time.getTimeZone() / 24.0, time.getTimeZone())
        spc.setTime(time_p1)
        riseTime, setTime, UNUSED = spc.getAngleTimes()
        return riseTime, setTime

def getSunTransitTime(time: JulianDate, geo: GeoPosition):
    spc = SunPositionController2(time, geo)
    vals = spc.getAngleTimes()
    return vals[2]


def getSunAltAz(time: JulianDate, geo: GeoPosition):
    spc = SunPositionController2(time, geo)
    sunPos = spc.getSunPosition()
    sunPosSez = toTopocentric(sunPos, time, geo)
    return degrees(asin(sunPosSez[2] / sunPosSez.mag())), degrees(atan3(sunPosSez[1], -sunPosSez[0]))


class SunInfo:

    def __init__(self, day: JulianDate, geo: GeoPosition):
        spc = SunPositionController2(day, geo)
        self._riseTime, self._setTime, self._transitTime = spc.getAngleTimes()
        self._dayLength = self._setTime.value() - self._riseTime.value()
        self._distance = spc.getSunDistance()
        self._transitAlt = spc.getTransitAltitude()
        self._civilTwilightStart, self._civilTwilightEnd, UNUSED = spc.getAngleTimes(-6.0)
        self._nauticalTwilightStart, self._nauticalTwilightEnd, UNUSED = spc.getAngleTimes(-12.0)
        self._astronomicalTwilightStart, self._astronomicalTwilightEnd, UNUSED = spc.getAngleTimes(-18.0)
        UNUSED, self._riseAzimuth = getSunAltAz(self._riseTime, geo)
        UNUSED, self._setAzimuth = getSunAltAz(self._setTime, geo)

    def __iter__(self):
        yield from {
            'riseTime': self._riseTime,
            'riseAzimuth': self._riseAzimuth,
            'setTime': self._setTime,
            'setAzimuth': self._setAzimuth,
            'transitTime': self._transitTime,
            'transitAltitude': self._transitAlt,
            'dayLength': self._dayLength,
            'distance': self._distance,
            'civilTwilightStart': self._civilTwilightStart,
            'civilTwilightEnd': self._civilTwilightEnd,
            'nauticalTwilightStart': self._nauticalTwilightStart,
            'nauticalTwilightEnd': self._nauticalTwilightEnd,
            'astronomicalTwilightStart': self._astronomicalTwilightStart,
            'astronomicalTwilightEnd': self._astronomicalTwilightEnd
        }

    def __str__(self):
        return f'Rise Time:  {SunInfo.__truncateTime(self._riseTime)} -- {"%.0f" % self._riseAzimuth}°\n' \
               f'Set Time:  {SunInfo.__truncateTime(self._setTime)} -- {"%.0f" % self._setAzimuth}°\n' \
               f'Day Length: {SunInfo.__convertDayLength(self._dayLength)}\n' \
               f'Dawn:\n' \
               f'Astronomical Twilight Start:  {SunInfo.__truncateTime(self._astronomicalTwilightStart)}\n' \
               f'Nautical Twilight Start:  {SunInfo.__truncateTime(self._nauticalTwilightStart)}\n' \
               f'Civil Twilight Start:  {SunInfo.__truncateTime(self._civilTwilightStart)}\n' \
               f'Dusk:\n' \
               f'Civil Twilight End:  {SunInfo.__truncateTime(self._civilTwilightEnd)}\n' \
               f'Nautical Twilight End:  {SunInfo.__truncateTime(self._nauticalTwilightEnd)}\n' \
               f'Astronomical Twilight End:  {SunInfo.__truncateTime(self._astronomicalTwilightEnd)}\n' \
               f'Transit Time:  {SunInfo.__truncateTime(self._transitTime)} -- {"%.0f" % self._transitAlt}°\n' \
               f'Distance:  {round(self._distance):,} km'

    def __repr__(self):
        return json.dumps(self.toJson(), default=lambda o: o.toJson())

    def toJson(self):
        rtn = {'riseAzimuth': self._riseAzimuth, 'setAzimuth': self._setAzimuth, 'transitAltitude': self._transitAlt,
               'dayLength': self._dayLength, 'distance': self._distance, 'riseTime': dict(self._riseTime),
               'setTime': dict(self._setTime), 'transitTime': dict(self._transitTime),
               'civilTwilightStart': dict(self._civilTwilightStart), 'civilTwilightEnd': dict(self._civilTwilightEnd),
               'nauticalTwilightStart': dict(self._nauticalTwilightStart),
               'nauticalTwilightEnd': dict(self._nauticalTwilightEnd),
               'astronomicalTwilightStart': dict(self._astronomicalTwilightStart),
               'astronomicalTwilightEnd': dict(self._astronomicalTwilightEnd)}
        return rtn

    @classmethod
    def __convertDayLength(cls, length: float) -> str:
        """Convert a solar day fraction into a time based string."""
        seconds = round(length * 86400)
        hours = int(seconds / 3600.0)
        seconds -= hours * 3600.0
        minutes = int(seconds / 60.0)
        seconds = int(seconds - minutes * 60.0)  # to remove trailing zeros in the string

        # hours = int(length * 24.0)
        # minutes = int((length - hours / 24.0) * 60.0)
        # seconds = int((length - (hours / 24.0) - (minutes / 60.0)) * 3500.0)
        hrStr = str(hours) if hours >= 10 else f'0{hours}'
        minStr = str(minutes) if minutes >= 10 else f'0{minutes}'
        secStr = str(seconds) if seconds >= 10 else f'0{seconds}'
        return f'{hrStr}:{minStr}:{secStr}'

    @classmethod
    def __truncateTime(cls, time: JulianDate) -> str:
        """Convert a julian date time to a string with the seconds rounded to an integer"""
        timeStr = time.time()
        sec = round(float(timeStr[6:]))
        secStr = str(sec) if sec >= 10 else f'0{sec}'
        return f'{timeStr[:6]}{secStr}'

    def getRiseTime(self):
        return self._riseTime

    def getSetTime(self):
        return self._setTime

    def getDayLength(self):
        return self._dayLength

    def getTransitTime(self):
        return self._transitTime

    def getTransitAltitude(self):
        return self._transitAlt

    def getDistance(self):
        return self._distance

    def getCivilTwilightStart(self):
        return self._civilTwilightStart

    def getCivilTwilightEnd(self):
        return self._civilTwilightEnd

    def getNauticalTwilightStart(self):
        return self._nauticalTwilightStart

    def getNauticalTwilightEnd(self):
        return self._nauticalTwilightEnd

    def getAstronomicalTwilightStart(self):
        return self._astronomicalTwilightStart

    def getAstronomicalTwilightEnd(self):
        return self._astronomicalTwilightEnd


# def getSunPosition2(time: JulianDate):
#     """
#     A more in depth algorithm for finding the Sun's position vector.
#
#     Args:
#         time: Time to find the Sun's position.
#
#     Returns:
#         The Sun's position vector in kilometers.
#     """
#
#     JDE = time.future(DELTAT / 86400)
#     # JC = time.difference(J2000) / 36525.0
#     JCE = JDE.difference(J2000) / 36525.0
#     JME = JCE / 10.0
#
#     L = expandTableValues(LTable, JME) % TWOPI
#     B = expandTableValues(BTable, JME) % TWOPI
#     R = expandTableValues(RTable, JME)
#
#     theta = (L + pi) % TWOPI
#     beta = -B
#     dPsi, dEpsilon = getNutationDeltas(JCE)
#     epsilon0 = 0.0
#     for i in range(11):
#         epsilon0 += UTable[i] * ((JME/10.0) ** i)
#     epsilon = radians(epsilon0 / 3600.0) + dEpsilon
#     dTau = -radians(20.4898 / (3600 * R))
#     lmbda = theta + dPsi + dTau
#     alpha = atan3(sin(lmbda) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lmbda))
#     delta = asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lmbda))
#
#     zComp = sin(delta)
#     xComp = sqrt((cos(delta) ** 2) / (1 + (tan(alpha) ** 2)))
#     if pi / 2 < alpha < 3*pi/2:
#         xComp = -abs(xComp)
#     else:
#         xComp = abs(xComp)
#
#     yComp = xComp * tan(alpha)
#     if alpha < pi:
#         yComp = abs(yComp)
#     else:
#         yComp = -abs(yComp)
#
#     '''
#     This code works, I only need the right ascension and declination for my calculations for now.
#     phi = radians(geo.getLatitude())
#     nu0Raw = 280.46061837 + 360.98564736629 * time.difference(J2000) + 0.000387933 * JC*JC - (JC*JC*JC / 38710000)
#     nu0 = radians(nu0Raw) % TWOPI
#     nu = nu0 + dPsi * cos(epsilon)
#     H = (nu + radians(geo.getLongitude()) - alpha) % TWOPI
#     ksi = radians(8.794 / (3600 * R))
#     u = atan(0.99664719 * tan(phi))
#     x = cos(u) + geo.getElevation() * cos(phi) / 6378140
#     y = 0.99664719 * sin(u) + geo.getElevation() * sin(phi) / 6378140
#     dAlpha = atan2(-x * sin(ksi) * sin(H), cos(delta) - x * sin(ksi) * cos(H))
#     alphaP = alpha + dAlpha
#     deltaP = atan2((sin(delta) - y*sin(ksi))*cos(dAlpha), cos(delta) - x*sin(ksi)*cos(H))
#     HP = H - dAlpha
#     e0 = asin(sin(phi)*sin(deltaP) + cos(phi)*cos(deltaP)*cos(HP))
#     # put atmospheric refraction correction here if temp and pressure daa is known
#     e = e0
#     zenith = pi/2 - e
#     azimuth = (atan2(sin(HP), cos(HP)*sin(phi) - tan(deltaP)*cos(phi)) + pi) % TWOPI
#     return azimuth, zenith'''
#     # todo: find correct way to compute R
#     n = time.difference(J2000)
#     g = 6.240040768 + 0.0172019703 * n
#     R = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2 * g)
#     return EVector(xComp, yComp, zComp) * R * AU


def getTwilightType(time: JulianDate, geo: GeoPosition) -> TwilightType:
    """
    Computes the twilight-type from the angle of the sun.

    Args:
        time: Time to find the twilight-type.
        geo: GeoPosition whose twilight-type is to be computed.

    Returns:
        One of the TwilightType enumerations.
    """

    spc = SunPositionController2(time, geo)
    # sunSEZPos = toTopocentric(getSunPosition2(time), time, geo)
    sunSEZPos = toTopocentric(spc.getSunPosition(), time, geo)
    sunAngle = degrees(asin(sunSEZPos[2] / sunSEZPos.mag()))
    if sunAngle < -18:
        return TwilightType.Night
    elif sunAngle < -12:
        return TwilightType.Astronomical
    elif sunAngle < -6:
        return TwilightType.Nautical
    elif sunAngle < -(5.0 / 6.0):
        return TwilightType.Civil
    else:
        return TwilightType.Day


class SunPositionController2:

    def __init__(self, time: JulianDate, geo: GeoPosition):
        self._DELTAT = 72.6  # seconds from https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm
        self._JD = time
        self._JDE = time.future(self._DELTAT / 86400.0)
        self._JC = time.difference(J2000) / 36525.0
        self._JC = time.difference(J2000) / 36525.0
        self._JCE = self._JDE.difference(J2000) / 36525.0
        self._JME = self._JCE / 10.0
        self._geo = geo
        self._internal = {'phi': radians(geo.getLatitude()), 'sigma': radians(geo.getLongitude())}
        self._LTable = [
            [
                [175347046, 0, 0],
                [3341656, 4.6692568, 6283.07585],
                [34894, 4.6261, 12566.1517],
                [3497, 2.7441, 5753.3849],
                [3418, 2.8289, 3.5231],
                [3136, 3.6277, 77713.7715],
                [2676, 4.4181, 7860.4194],
                [2343, 6.1352, 3930.2097],
                [1324, 0.7425, 11506.7698],
                [1273, 2.0371, 529.691],
                [1199, 1.1096, 1577.3435],
                [990, 5.233, 5884.927],
                [902, 2.045, 26.298],
                [857, 3.508, 398.149],
                [780, 1.179, 5223.694],
                [753, 2.533, 5507.553],
                [505, 4.583, 18849.228],
                [492, 4.205, 775.523],
                [357, 2.92, 0.067],
                [317, 5.849, 11790.629],
                [284, 1.899, 796.298],
                [271, 0.315, 10977.079],
                [243, 0.345, 5486.778],
                [206, 4.806, 2544.314],
                [205, 1.869, 5573.143],
                [202, 2.458, 6069.777],
                [156, 0.833, 213.299],
                [132, 3.411, 2942.463],
                [126, 1.083, 20.775],
                [115, 0.645, 0.98],
                [103, 0.636, 4694.003],
                [102, 0.976, 15720.839],
                [102, 4.267, 7.114],
                [99, 6.21, 2146.17],
                [98, 0.68, 155.42],
                [86, 5.98, 161000.69],
                [85, 1.3, 6275.96],
                [85, 3.67, 71430.7],
                [80, 1.81, 17260.15],
                [79, 3.04, 12036.46],
                [75, 1.76, 5088.63],
                [74, 3.5, 3154.69],
                [74, 4.68, 801.82],
                [70, 0.83, 9437.76],
                [62, 3.98, 8827.39],
                [61, 1.82, 7084.9],
                [57, 2.78, 6286.6],
                [56, 4.39, 14143.5],
                [56, 3.47, 6279.55],
                [52, 0.19, 12139.55],
                [52, 1.33, 1748.02],
                [51, 0.28, 5856.48],
                [49, 0.49, 1194.45],
                [41, 5.37, 8429.24],
                [41, 2.4, 19651.05],
                [39, 6.17, 10447.39],
                [37, 6.04, 10213.29],
                [37, 2.57, 1059.38],
                [36, 1.71, 2352.87],
                [36, 1.78, 6812.77],
                [33, 0.59, 17789.85],
                [30, 0.44, 83996.85],
                [30, 2.74, 1349.87],
                [25, 3.16, 4690.48]
            ],
            [
                [628331966747.0, 0, 0],
                [206059, 2.678235, 6283.07585],
                [4303, 2.6351, 12566.1517],
                [425, 1.59, 3.523],
                [119, 5.769, 26.298],
                [109, 2.966, 1577.344],
                [93, 2.59, 18849.23],
                [72, 1.14, 529.69],
                [68, 1.87, 398.15],
                [67, 4.41, 5507.55],
                [59, 2.89, 5223.69],
                [56, 2.17, 155.42],
                [45, 0.4, 796.3],
                [36, 0.47, 775.52],
                [29, 2.65, 7.11],
                [21, 5.34, 0.98],
                [19, 1.85, 5486.78],
                [19, 4.97, 213.3],
                [17, 2.99, 6275.96],
                [16, 0.03, 2544.31],
                [16, 1.43, 2146.17],
                [15, 1.21, 10977.08],
                [12, 2.83, 1748.02],
                [12, 3.26, 5088.63],
                [12, 5.27, 1194.45],
                [12, 2.08, 4694],
                [11, 0.77, 553.57],
                [10, 1.3, 6286.6],
                [10, 4.24, 1349.87],
                [9, 2.7, 242.73],
                [9, 5.64, 951.72],
                [8, 5.3, 2352.87],
                [6, 2.65, 9437.76],
                [6, 4.67, 4690.48]
            ],
            [
                [52919, 0, 0],
                [8720, 1.0721, 6283.0758],
                [309, 0.867, 12566.152],
                [27, 0.05, 3.52],
                [16, 5.19, 26.3],
                [16, 3.68, 155.42],
                [10, 0.76, 18849.23],
                [9, 2.06, 77713.77],
                [7, 0.83, 775.52],
                [5, 4.66, 1577.34],
                [4, 1.03, 7.11],
                [4, 3.44, 5573.14],
                [3, 5.14, 796.3],
                [3, 6.05, 5507.55],
                [3, 1.19, 242.73],
                [3, 6.12, 529.69],
                [3, 0.31, 398.15],
                [3, 2.28, 553.57],
                [2, 4.38, 5223.69],
                [2, 3.75, 0.98]
            ],
            [
                [289, 5.844, 6283.076],
                [35, 0, 0],
                [17, 5.49, 12566.15],
                [3, 5.2, 155.42],
                [1, 4.72, 3.52],
                [1, 5.3, 18849.23],
                [1, 5.97, 242.73]
            ],
            [
                [114, 3.142, 0],
                [8, 4.13, 6283.08],
                [1, 3.84, 12566.15]
            ],
            [
                [1, 3.14, 0]
            ]
        ]

        self._BTable = [
            [
                [280, 3.199, 84334.662],
                [102, 5.422, 5507.553],
                [80, 3.88, 5223.69],
                [44, 3.7, 2352.87],
                [32, 4, 1577.34]
            ],
            [
                [9, 3.9, 5507.55],
                [6, 1.73, 5223.69]
            ]
        ]

        self._RTable = [
            [
                [100013989, 0, 0],
                [1670700, 3.0984635, 6283.07585],
                [13956, 3.05525, 12566.1517],
                [3084, 5.1985, 77713.7715],
                [1628, 1.1739, 5753.3849],
                [1576, 2.8469, 7860.4194],
                [925, 5.453, 11506.77],
                [542, 4.564, 3930.21],
                [472, 3.661, 5884.927],
                [346, 0.964, 5507.553],
                [329, 5.9, 5223.694],
                [307, 0.299, 5573.143],
                [243, 4.273, 11790.629],
                [212, 5.847, 1577.344],
                [186, 5.022, 10977.079],
                [175, 3.012, 18849.228],
                [110, 5.055, 5486.778],
                [98, 0.89, 6069.78],
                [86, 5.69, 15720.84],
                [86, 1.27, 161000.69],
                [65, 0.27, 7260.15],
                [63, 0.92, 529.69],
                [57, 2.01, 83996.85],
                [56, 5.24, 71430.7],
                [49, 3.25, 2544.31],
                [47, 2.58, 775.52],
                [45, 5.54, 9437.76],
                [43, 6.01, 6275.96],
                [39, 5.36, 4694],
                [38, 2.39, 8827.39],
                [37, 0.83, 19651.05],
                [37, 4.9, 12139.55],
                [36, 1.67, 12036.46],
                [35, 1.84, 2942.46],
                [33, 0.24, 7084.9],
                [32, 0.18, 5088.63],
                [32, 1.78, 398.15],
                [28, 1.21, 6286.6],
                [28, 1.9, 6279.55],
                [26, 4.59, 10447.39]
            ],
            [
                [103019, 1.10749, 6283.07585],
                [1721, 1.0644, 12566.1517],
                [702, 3.142, 0],
                [32, 1.02, 18849.23],
                [31, 2.84, 5507.55],
                [25, 1.32, 5223.69],
                [18, 1.42, 1577.34],
                [10, 5.91, 10977.08],
                [9, 1.42, 6275.96],
                [9, 0.27, 5486.78]
            ],
            [
                [4359, 5.7846, 6283.0758],
                [124, 5.579, 12566.152],
                [12, 3.14, 0],
                [9, 3.63, 77713.77],
                [6, 1.87, 5573.14],
                [3, 5.47, 18849.23]
            ],
            [
                [145, 4.273, 6283.076],
                [7, 3.92, 12566.15]
            ],
            [
                [4, 2.56, 6283.08]
            ]
        ]

        self._XTable = [
            [297.85036, 445267.111480, -0.0019142, 1.0 / 189474.0],
            [357.527772, 35999.050340, -0.0001603, 1.0 / -300000.0],
            [134.96298, 477198.867398, 0.0086972, 1.0 / 56250.0],
            [93.27191, 483202.017538, -0.0036825, 1.0 / 327270.0],
            [125.04452, -1934.136261, 0.0020708, 1.0 / 450000.0]
        ]

        self._NutationTable = [
            [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
            [-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1],
            [0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5],
            [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
            [0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1],
            [0, 0, 1, 0, 0, 712, 0.1, -7, 0],
            [-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6],
            [0, 0, 0, 2, 1, -386, -0.4, 200, 0],
            [0, 0, 1, 2, 2, -301, 0, 129, -0.1],
            [-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3],
            [-2, 0, 1, 0, 0, -158, 0, 0, 0],
            [-2, 0, 0, 2, 1, 129, 0.1, -70, 0],
            [0, 0, -1, 2, 2, 123, 0, -53, 0],
            [2, 0, 0, 0, 0, 63, 0, 0, 0],
            [0, 0, 1, 0, 1, 63, 0.1, -33, 0],
            [2, 0, -1, 2, 2, -59, 0, 26, 0],
            [0, 0, -1, 0, 1, -58, -0.1, 32, 0],
            [0, 0, 1, 2, 1, -51, 0, 27, 0],
            [-2, 0, 2, 0, 0, 48, 0, 0, 0],
            [0, 0, -2, 2, 1, 46, 0, -24, 0],
            [2, 0, 0, 2, 2, -38, 0, 16, 0],
            [0, 0, 2, 2, 2, -31, 0, 13, 0],
            [0, 0, 2, 0, 0, 29, 0, 0, 0],
            [-2, 0, 1, 2, 2, 29, 0, -12, 0],
            [0, 0, 0, 2, 0, 26, 0, 0, 0],
            [-2, 0, 0, 2, 0, -22, 0, 0, 0],
            [0, 0, -1, 2, 1, 21, 0, -10, 0],
            [0, 2, 0, 0, 0, 17, -0.1, 0, 0],
            [2, 0, -1, 0, 1, 16, 0, -8, 0],
            [-2, 2, 0, 2, 2, -16, 0.1, 7, 0],
            [0, 1, 0, 0, 1, -15, 0, 9, 0],
            [-2, 0, 1, 0, 1, -13, 0, 7, 0],
            [0, -1, 0, 0, 1, -12, 0, 6, 0],
            [0, 0, 2, -2, 0, 11, 0, 0, 0],
            [2, 0, -1, 2, 1, -10, 0, 5, 0],
            [2, 0, 1, 2, 2, -8, 0, 3, 0],
            [0, 1, 0, 2, 2, 7, 0, -3, 0],
            [-2, 1, 1, 0, 0, -7, 0, 0, 0],
            [0, -1, 0, 2, 2, -7, 0, 3, 0],
            [2, 0, 0, 2, 1, -7, 0, 3, 0],
            [2, 0, 1, 0, 0, 6, 0, 0, 0],
            [-2, 0, 2, 2, 2, 6, 0, -3, 0],
            [-2, 0, 1, 2, 1, 6, 0, -3, 0],
            [2, 0, -2, 0, 1, -6, 0, 3, 0],
            [2, 0, 0, 0, 1, -6, 0, 3, 0],
            [0, -1, 1, 0, 0, 5, 0, 0, 0],
            [-2, -1, 0, 2, 1, -5, 0, 3, 0],
            [-2, 0, 0, 0, 1, -5, 0, 3, 0],
            [0, 0, 2, 2, 1, -5, 0, 3, 0],
            [-2, 0, 2, 0, 1, 4, 0, 0, 0],
            [-2, 1, 0, 2, 1, 4, 0, 0, 0],
            [0, 0, 1, -2, 0, 4, 0, 0, 0],
            [-1, 0, 1, 0, 0, -4, 0, 0, 0],
            [-2, 1, 0, 0, 0, -4, 0, 0, 0],
            [1, 0, 0, 0, 0, -4, 0, 0, 0],
            [0, 0, 1, 2, 0, 3, 0, 0, 0],
            [0, 0, -2, 2, 2, -3, 0, 0, 0],
            [-1, -1, 1, 0, 0, -3, 0, 0, 0],
            [0, 1, 1, 0, 0, -3, 0, 0, 0],
            [0, -1, 1, 2, 2, -3, 0, 0, 0],
            [2, -1, -1, 2, 2, -3, 0, 0, 0],
            [0, 0, 3, 2, 2, -3, 0, 0, 0],
            [2, -1, 0, 2, 2, -3, 0, 0, 0]
        ]

        self._UTable = [
            84381.448, -4680.93, -1.55, 1999.25,
            51.38, -249.67, -39.05, 7.12,
            27.87, 5.79, 2.45
        ]

    def getTime(self) -> JulianDate:
        return self._JD

    def setTime(self, time: JulianDate):
        self._JD = time
        self._JDE = time.future(self._DELTAT / 86400.0)
        self._JC = time.difference(J2000) / 36525.0
        self._JCE = self._JDE.difference(J2000) / 36525.0
        self._JME = self._JCE / 10.0
        self._internal = {'phi': radians(self._geo.getLatitude()), 'sigma': radians(self._geo.getLongitude())}

    def getSunPosition(self):
        try:
            ra = self._internal['alpha']
        except KeyError:
            self.__computeSunRA()
            ra = self._internal['alpha']
        try:
            dec = self._internal['delta']
        except KeyError:
            self.__computeSunDec()
            dec = self._internal['delta']
        try:
            R = self._internal['R']
        except KeyError:
            self.__computeHelioRadius()
            R = self._internal['R']

        zComp = sin(dec)
        xComp = sqrt((cos(dec) ** 2.0) / (1 + (tan(ra) ** 2.0)))
        xComp = -abs(xComp) if pi / 2 < ra < 3 * pi / 2 else abs(xComp)
        yComp = xComp * tan(ra)
        yComp = abs(yComp) if ra < pi else -abs(yComp)

        return EVector(xComp, yComp, zComp) * R * AU

    def getSunCoordinates(self):
        try:
            ra = self._internal['alpha']
        except KeyError:
            self.__computeSunRA()
            ra = self._internal['alpha']
        try:
            dec = self._internal['delta']
        except KeyError:
            self.__computeSunDec()
            dec = self._internal['delta']
        return ra * 12 / pi, degrees(dec)

    def getAngleTimes(self, target=None):
        """target in degrees"""
        tzOffset = self._JD.getTimeZone() / 24.0
        localVal = self._JD.value() + tzOffset
        num = int(localVal)
        if localVal - int(localVal) < 0.5:
            num -= 1
        time = JulianDate.fromNumber(num, 0.5 - tzOffset, self._JD.getTimeZone())
        spc_time = SunPositionController2(time, self._geo)
        v = spc_time.getApparentSiderealTime()
        if target is None:
            target = -(5.0 / 6.0)
        dt = self._DELTAT / 86400.0
        time_m1 = time.future(dt - 1)
        time_0 = time.future(dt)
        time_p1 = time.future(dt + 1)

        spc_m1 = SunPositionController2(time_m1, self._geo)
        spc_0 = SunPositionController2(time_0, self._geo)
        spc_p1 = SunPositionController2(time_p1, self._geo)
        vals_m1 = spc_m1.getSunCoordinates()
        alpha_m1 = vals_m1[0] * pi / 12.0
        delta_m1 = vals_m1[1] * pi / 180.0
        vals_0 = spc_0.getSunCoordinates()
        alpha_0 = vals_0[0] * pi / 12.0
        delta_0 = vals_0[1] * pi / 180.0
        vals_p1 = spc_p1.getSunCoordinates()
        alpha_p1 = vals_p1[0] * pi / 12.0
        delta_p1 = vals_p1[1] * pi / 180.0

        sigma = self._internal['sigma']
        phi = self._internal['phi']
        m0 = (alpha_0 - sigma - v) / TWOPI
        hp0 = radians(target)
        try:
            H0 = acos((sin(hp0) - sin(phi) * sin(delta_0)) / (cos(phi) * cos(delta_0))) % pi
        except ValueError:
            raise SunRiseSetException('Sun does not rise or set during this time.')
        m1 = (m0 - (H0 / TWOPI)) % 1.0
        m2 = (m0 + (H0 / TWOPI)) % 1.0

        v0 = v + radians(360.985647) * m0
        v1 = v + radians(360.985647) * m1
        v2 = v + radians(360.985647) * m2

        n0 = m0 + dt
        n1 = m1 + dt
        n2 = m2 + dt

        # this is unchecked, if problems exist, check here first
        maxVal = radians(2.0)
        modulator = radians(1.0)
        a = alpha_0 - alpha_m1
        if abs(a) > maxVal:
            a %= modulator
        b = alpha_p1 - alpha_0
        if abs(b) > maxVal:
            b %= modulator
        ap = delta_0 - delta_m1
        if abs(ap) > maxVal:
            ap %= modulator
        bp = delta_p1 - delta_0
        if abs(bp) > maxVal:
            bp %= modulator
        c = b - a
        cp = bp - ap
        alphaP0 = alpha_0 + ((n0 * (a + b + c * n0)) / 2.0)
        alphaP1 = alpha_0 + ((n1 * (a + b + c * n1)) / 2.0)
        alphaP2 = alpha_0 + ((n2 * (a + b + c * n2)) / 2.0)
        deltaP0 = delta_0 + ((n0 * (ap + bp + cp * n0)) / 2.0)
        deltaP1 = delta_0 + ((n1 * (ap + bp + cp * n1)) / 2.0)
        deltaP2 = delta_0 + ((n2 * (ap + bp + cp * n2)) / 2.0)

        Hp0 = (v0 + sigma - alphaP0) % TWOPI
        if Hp0 >= pi:
            Hp0 -= TWOPI
        Hp1 = (v1 + sigma - alphaP1) % TWOPI
        if Hp1 >= pi:
            Hp1 -= TWOPI
        Hp2 = (v2 + sigma - alphaP2) % TWOPI
        if Hp2 >= pi:
            Hp2 -= TWOPI

        h0 = asin(sin(phi) * sin(deltaP0) + cos(phi) * cos(deltaP0) * cos(Hp0))  # sun altitude at transit time
        self._internal['h0'] = h0
        h1 = asin(sin(phi) * sin(deltaP1) + cos(phi) * cos(deltaP1) * cos(Hp1))
        h2 = asin(sin(phi) * sin(deltaP2) + cos(phi) * cos(deltaP2) * cos(Hp2))

        T = m0 - (Hp0 / TWOPI)
        R = m1 + ((h1 - hp0) / (TWOPI * cos(deltaP1) * cos(phi) * sin(Hp1)))
        S = m2 + ((h2 - hp0) / (TWOPI * cos(deltaP2) * cos(phi) * sin(Hp2)))
        return time.future(R), time.future(S), time.future(T)

    def getApparentSiderealTime(self):
        try:
            v = self._internal['v']
        except KeyError:
            self.__computeAppST()
            v = self._internal['v']
        return v

    def getTrueObliquity(self):
        """in radians"""
        try:
            epsilon = self._internal['epsilon']
        except KeyError:
            self.__computeTrueObliquity()
            epsilon = self._internal['epsilon']
        return epsilon

    def getMeanObliquity(self):
        """in radians"""
        try:
            epsilon0 = self._internal['epsilon0']
        except KeyError:
            self.__computeTrueObliquity()
            epsilon0 = self._internal['epsilon0']
        return epsilon0 / 3600.0

    def getHourAngle(self):
        """topocentric local hour angle"""
        try:
            HP = self._internal['HP']
        except KeyError:
            self.__computeTopoLocalHA()
            HP = self._internal['HP']
        return HP

    def getTransitAltitude(self):
        try:
            h0 = self._internal['h0']
        except KeyError:
            self.getAngleTimes()
            h0 = self._internal['h0']
        return degrees(h0)

    def getSunDistance(self):
        try:
            R = self._internal['R']
        except KeyError:
            self.__computeHelioRadius()
            R = self._internal['R']
        return R * AU

    def __expandTableValues(self, table):
        arr = [self.__tableSum(ti) for ti in table]
        sums = 0.0
        for i in range(len(table)):
            sums += arr[i] * (self._JME ** i)
        return sums / 1e8

    def __tableSum(self, table):
        sums = 0.0
        for ti in table:
            sub = ti[0] * cos(ti[1] + ti[2] * self._JME)
            sums += sub
        return sums

    def __getXArray(self):
        X = []
        for i in range(5):
            Xi = 0.0
            for j in range(4):
                Xi += self._XTable[i][j] * (self._JCE ** j)
            X.append(radians(Xi))
        return X

    def __computeHelioLong(self):
        """Guaranteed to be set: L  """
        self._internal['L'] = self.__expandTableValues(self._LTable) % TWOPI

    def __computeHelioLat(self):
        """Guaranteed to be set: B  """
        self._internal['B'] = self.__expandTableValues(self._BTable) % TWOPI

    def __computeHelioRadius(self):
        """Guaranteed to be set: R  """
        self._internal['R'] = self.__expandTableValues(self._RTable)

    def __computeGeoLong(self):
        """Guaranteed to be set: L, theta   """
        L = self.__expandTableValues(self._LTable)
        self._internal['L'] = L % TWOPI
        self._internal['theta'] = (L + pi) % TWOPI

    def __computeGeoLat(self):
        """Guaranteed to be set: B, beta    """
        try:
            B = self._internal['B']
        except KeyError:
            self.__computeHelioLat()
            B = self._internal['B']
        self._internal['beta'] = -B

    def __computeNutationDeltas(self):
        """Guaranteed to be set: dPsi, dEpsilon """
        arrX = self.__getXArray()
        dPsi, dEpsilon = 0.0, 0.0
        for ti in self._NutationTable:
            sumTerm = 0.0
            for i in range(5):
                sumTerm += arrX[i] * ti[i]
            dPsi += (ti[5] + ti[6] * self._JCE) * sin(sumTerm)
            dEpsilon += (ti[7] + ti[8] * self._JCE) * cos(sumTerm)
        self._internal['dPsi'] = radians(dPsi / 3.6e7)
        self._internal['dEpsilon'] = radians(dEpsilon / 3.6e7)

    def __computeTrueObliquity(self):
        """Guaranteed to be set: dPsi, dEpsilon, epsilon0, epsilon"""
        epsilon0 = 0.0
        U = self._JME / 10.0
        for i in range(11):
            epsilon0 += self._UTable[i] * (U ** i)
        try:
            dEpsilon = self._internal['dEpsilon']
        except KeyError:
            self.__computeNutationDeltas()
            dEpsilon = self._internal['dEpsilon']
        self._internal['epsilon0'] = epsilon0
        self._internal['epsilon'] = radians(epsilon0 / 3600.0) + dEpsilon

    def __computeAbrCor(self):
        """Guaranteed to be set: R, dTau    """
        try:
            R = self._internal['R']
        except KeyError:
            self.__computeHelioRadius()
            R = self._internal['R']
        self._internal['dTau'] = -radians(20.4898 / (3600.0 * R))

    def __computeAppSolarLong(self):
        """Guaranteed to be set: theta, L, dPsi, dEpsilon, R, dTau, lambda """
        try:
            theta = self._internal['theta']
        except KeyError:
            self.__computeGeoLong()
            theta = self._internal['theta']
        try:
            dPsi = self._internal['dPsi']
        except KeyError:
            self.__computeNutationDeltas()
            dPsi = self._internal['dPsi']
        try:
            dTau = self._internal['dTau']
        except KeyError:
            self.__computeAbrCor()
            dTau = self._internal['dTau']
        self._internal['lambda'] = theta + dPsi + dTau

    def __computeAppST(self):
        """Guaranteed to be set: dPsi, dEpsilon, epsilon0, epsilon, v   """
        try:
            epsilon = self._internal['epsilon']
            dPsi = self._internal['dPsi']
        except KeyError:
            self.__computeTrueObliquity()
            epsilon = self._internal['epsilon']
            dPsi = self._internal['dPsi']
        v0 = 4.894961212735793 + 6.300388098984957 * (self._JD.value() - 2451545) + 6.770708127139162e-06 \
            * (self._JC ** 2.0) - ((self._JC ** 3.0) / 675616.953447005)
        self._internal['v'] = (v0 % TWOPI) + dPsi * cos(epsilon)

    def __computeSunRA(self):
        """Guaranteed to be set: theta, L, dPsi, dEpsilon, R, dTau, lambda, epsilon0, epsilon, alpha """
        try:
            lmbda = self._internal['lambda']
        except KeyError:
            self.__computeAppSolarLong()
            lmbda = self._internal['lambda']
        try:
            epsilon = self._internal['epsilon']
        except KeyError:
            self.__computeTrueObliquity()
            epsilon = self._internal['epsilon']
        try:
            beta = self._internal['beta']
        except KeyError:
            self.__computeGeoLat()
            beta = self._internal['beta']
        self._internal['alpha'] = atan3(sin(lmbda) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lmbda))

    def __computeSunDec(self):
        """Guaranteed to be set: theta, L, dPsi, dEpsilon, R, dTau, lambda, epsilon0, epsilon, delta """
        try:
            lmbda = self._internal['lambda']
        except KeyError:
            self.__computeAppSolarLong()
            lmbda = self._internal['lambda']
        try:
            epsilon = self._internal['epsilon']
        except KeyError:
            self.__computeTrueObliquity()
            epsilon = self._internal['epsilon']
        try:
            beta = self._internal['beta']
        except KeyError:
            self.__computeGeoLat()
            beta = self._internal['beta']
        self._internal['delta'] = asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lmbda))

    def __computeLocalHA(self):
        """Guaranteed to be set: theta, L, dPsi, dEpsilon, R, dTau, lambda, epsilon0, epsilon, alpha, v, H """
        try:
            v = self._internal['v']
        except KeyError:
            self.__computeAppST()
            v = self._internal['v']
        sigma = self._internal['sigma']
        try:
            alpha = self._internal['alpha']
        except KeyError:
            self.__computeSunRA()
            alpha = self._internal['alpha']
        self._internal['H'] = (v + sigma - alpha) % TWOPI

    def __computeTopoSunRADec(self):
        """Guaranteed to be set: R, theta, L, dPsi, dEpsilon, dTau, lambda, epsilon0, epsilon, alpha, v, H, delta  """
        try:
            R = self._internal['R']
        except KeyError:
            self.__computeHelioRadius()
            R = self._internal['R']
        try:
            H = self._internal['H']
        except KeyError:
            self.__computeLocalHA()
            H = self._internal['H']
        try:
            alpha = self._internal['alpha']
        except KeyError:
            self.__computeSunRA()
            alpha = self._internal['alpha']
        try:
            delta = self._internal['delta']
        except KeyError:
            self.__computeSunDec()
            delta = self._internal['delta']
        zeta = radians(8.794 / (3600.0 * R))
        u = atan(0.99664719 * tan(self._internal['phi']))
        E = self._geo.getElevation()
        x = cos(u) + (E / 6378140) * cos(self._internal['phi'])
        y = 0.99664719 * sin(u) + (E / 6378140) * sin(self._internal['phi'])
        dAlpha = atan2(-x * sin(zeta) * sin(H), cos(delta) - x * sin(zeta) * cos(H))
        self._internal['alphaP'] = alpha + dAlpha
        self._internal['deltaP'] = atan2((sin(delta) - y * sin(zeta)) * cos(dAlpha),
                                         cos(delta) - x * sin(zeta) * cos(H))

    def __computeTopoLocalHA(self):
        """Guaranteed to be set: theta, L, dPsi, dEpsilon, R, dTau, lambda, epsilon0, epsilon, alpha, v, H, delta,
        HP """
        try:
            H = self._internal['H']
        except KeyError:
            self.__computeLocalHA()
            H = self._internal['H']
        try:
            R = self._internal['R']
        except KeyError:
            self.__computeHelioRadius()
            R = self._internal['R']
        try:
            delta = self._internal['delta']
        except KeyError:
            self.__computeSunDec()
            delta = self._internal['delta']
        zeta = radians(8.794 / (3600.0 * R))
        u = atan(0.99664719 * tan(self._internal['phi']))
        E = self._geo.getElevation()
        x = cos(u) + (E / 6378140) * cos(self._internal['phi'])
        dAlpha = atan2(-x * sin(zeta) * sin(H), cos(delta) - x * sin(zeta) * cos(H))
        self._internal['HP'] = H - dAlpha

    def __computeTopoZenith(self):

        try:
            deltaP = self._internal['deltaP']
        except KeyError:
            self.__computeTopoSunRADec()
            deltaP = self._internal['deltaP']
        try:
            HP = self._internal['HP']
        except KeyError:
            self.__computeTopoLocalHA()
            HP = self._internal['HP']
        phi = self._internal['phi']
        e0 = asin(sin(phi) * sin(deltaP) + cos(phi) * cos(deltaP) * cos(HP))
        # todo: include the deltaE term when we can pass pressure and temperature values
        # dE = (P / 1010) * (283 / (273 + T)) * (1.02 / (60 * tan(e0 + (10.3 / (e0 + 5.11)))))
        # self._internal['e'] = e0 + dE
        self._internal['zenith'] = (pi / 2.0) - e0

    def __computeTopoAz(self):

        try:
            HP = self._internal['HP']
        except KeyError:
            self.__computeTopoLocalHA()
            HP = self._internal['HP']
        try:
            deltaP = self._internal['deltaP']
        except KeyError:
            self.__computeTopoSunRADec()
            deltaP = self._internal['deltaP']
        phi = self._internal['phi']
        # this should work, otherwise use atan3, and limit final answer by modding TWOPI
        self._internal['azimuth'] = atan2(sin(HP), cos(HP) * sin(phi) - tan(deltaP) * cos(phi)) + pi

    def __computeIncidence(self):

        #   dont use this, dont have omega or gamma yet
        try:
            GAMMA = self._internal['azimuth'] - pi
        except KeyError:
            self.__computeTopoAz()
            GAMMA = self._internal['azimuth'] - pi
        theta = self._internal['zenith']
        # todo: figure out how to handle this
        omega = 0
        gamma = 0
        self._internal['I'] = acos(cos(theta) * cos(omega) + sin(omega) * sin(theta) * cos(GAMMA - gamma))

# class SunPositionController:
#
#     def __init__(self):
#         self._DELTAT = 72.6 # seconds from https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm
#         self._JD = None
#         self._JDE = None
#         self._JC = None
#         self._JCE = None
#         self._JME = None
#         self._LTable = [
#             [
#                 [175347046, 0, 0],
#                 [3341656, 4.6692568, 6283.07585],
#                 [34894, 4.6261, 12566.1517],
#                 [3497, 2.7441, 5753.3849],
#                 [3418, 2.8289, 3.5231],
#                 [3136, 3.6277, 77713.7715],
#                 [2676, 4.4181, 7860.4194],
#                 [2343, 6.1352, 3930.2097],
#                 [1324, 0.7425, 11506.7698],
#                 [1273, 2.0371, 529.691],
#                 [1199, 1.1096, 1577.3435],
#                 [990, 5.233, 5884.927],
#                 [902, 2.045, 26.298],
#                 [857, 3.508, 398.149],
#                 [780, 1.179, 5223.694],
#                 [753, 2.533, 5507.553],
#                 [505, 4.583, 18849.228],
#                 [492, 4.205, 775.523],
#                 [357, 2.92, 0.067],
#                 [317, 5.849, 11790.629],
#                 [284, 1.899, 796.298],
#                 [271, 0.315, 10977.079],
#                 [243, 0.345, 5486.778],
#                 [206, 4.806, 2544.314],
#                 [205, 1.869, 5573.143],
#                 [202, 2.458, 6069.777],
#                 [156, 0.833, 213.299],
#                 [132, 3.411, 2942.463],
#                 [126, 1.083, 20.775],
#                 [115, 0.645, 0.98],
#                 [103, 0.636, 4694.003],
#                 [102, 0.976, 15720.839],
#                 [102, 4.267, 7.114],
#                 [99, 6.21, 2146.17],
#                 [98, 0.68, 155.42],
#                 [86, 5.98, 161000.69],
#                 [85, 1.3, 6275.96],
#                 [85, 3.67, 71430.7],
#                 [80, 1.81, 17260.15],
#                 [79, 3.04, 12036.46],
#                 [75, 1.76, 5088.63],
#                 [74, 3.5, 3154.69],
#                 [74, 4.68, 801.82],
#                 [70, 0.83, 9437.76],
#                 [62, 3.98, 8827.39],
#                 [61, 1.82, 7084.9],
#                 [57, 2.78, 6286.6],
#                 [56, 4.39, 14143.5],
#                 [56, 3.47, 6279.55],
#                 [52, 0.19, 12139.55],
#                 [52, 1.33, 1748.02],
#                 [51, 0.28, 5856.48],
#                 [49, 0.49, 1194.45],
#                 [41, 5.37, 8429.24],
#                 [41, 2.4, 19651.05],
#                 [39, 6.17, 10447.39],
#                 [37, 6.04, 10213.29],
#                 [37, 2.57, 1059.38],
#                 [36, 1.71, 2352.87],
#                 [36, 1.78, 6812.77],
#                 [33, 0.59, 17789.85],
#                 [30, 0.44, 83996.85],
#                 [30, 2.74, 1349.87],
#                 [25, 3.16, 4690.48]
#             ],
#             [
#                 [628331966747.0, 0, 0],
#                 [206059, 2.678235, 6283.07585],
#                 [4303, 2.6351, 12566.1517],
#                 [425, 1.59, 3.523],
#                 [119, 5.769, 26.298],
#                 [109, 2.966, 1577.344],
#                 [93, 2.59, 18849.23],
#                 [72, 1.14, 529.69],
#                 [68, 1.87, 398.15],
#                 [67, 4.41, 5507.55],
#                 [59, 2.89, 5223.69],
#                 [56, 2.17, 155.42],
#                 [45, 0.4, 796.3],
#                 [36, 0.47, 775.52],
#                 [29, 2.65, 7.11],
#                 [21, 5.34, 0.98],
#                 [19, 1.85, 5486.78],
#                 [19, 4.97, 213.3],
#                 [17, 2.99, 6275.96],
#                 [16, 0.03, 2544.31],
#                 [16, 1.43, 2146.17],
#                 [15, 1.21, 10977.08],
#                 [12, 2.83, 1748.02],
#                 [12, 3.26, 5088.63],
#                 [12, 5.27, 1194.45],
#                 [12, 2.08, 4694],
#                 [11, 0.77, 553.57],
#                 [10, 1.3, 6286.6],
#                 [10, 4.24, 1349.87],
#                 [9, 2.7, 242.73],
#                 [9, 5.64, 951.72],
#                 [8, 5.3, 2352.87],
#                 [6, 2.65, 9437.76],
#                 [6, 4.67, 4690.48]
#             ],
#             [
#                 [52919, 0, 0],
#                 [8720, 1.0721, 6283.0758],
#                 [309, 0.867, 12566.152],
#                 [27, 0.05, 3.52],
#                 [16, 5.19, 26.3],
#                 [16, 3.68, 155.42],
#                 [10, 0.76, 18849.23],
#                 [9, 2.06, 77713.77],
#                 [7, 0.83, 775.52],
#                 [5, 4.66, 1577.34],
#                 [4, 1.03, 7.11],
#                 [4, 3.44, 5573.14],
#                 [3, 5.14, 796.3],
#                 [3, 6.05, 5507.55],
#                 [3, 1.19, 242.73],
#                 [3, 6.12, 529.69],
#                 [3, 0.31, 398.15],
#                 [3, 2.28, 553.57],
#                 [2, 4.38, 5223.69],
#                 [2, 3.75, 0.98]
#             ],
#             [
#                 [289, 5.844, 6283.076],
#                 [35, 0, 0],
#                 [17, 5.49, 12566.15],
#                 [3, 5.2, 155.42],
#                 [1, 4.72, 3.52],
#                 [1, 5.3, 18849.23],
#                 [1, 5.97, 242.73]
#             ],
#             [
#                 [114, 3.142, 0],
#                 [8, 4.13, 6283.08],
#                 [1, 3.84, 12566.15]
#             ],
#             [
#                 [1, 3.14, 0]
#             ]
#         ]
#
#         self._BTable = [
#             [
#                 [280, 3.199, 84334.662],
#                 [102, 5.422, 5507.553],
#                 [80, 3.88, 5223.69],
#                 [44, 3.7, 2352.87],
#                 [32, 4, 1577.34]
#             ],
#             [
#                 [9, 3.9, 5507.55],
#                 [6, 1.73, 5223.69]
#             ]
#         ]
#
#         self._RTable = [
#             [
#                 [100013989, 0, 0],
#                 [1670700, 3.0984635, 6283.07585],
#                 [13956, 3.05525, 12566.1517],
#                 [3084, 5.1985, 77713.7715],
#                 [1628, 1.1739, 5753.3849],
#                 [1576, 2.8469, 7860.4194],
#                 [925, 5.453, 11506.77],
#                 [542, 4.564, 3930.21],
#                 [472, 3.661, 5884.927],
#                 [346, 0.964, 5507.553],
#                 [329, 5.9, 5223.694],
#                 [307, 0.299, 5573.143],
#                 [243, 4.273, 11790.629],
#                 [212, 5.847, 1577.344],
#                 [186, 5.022, 10977.079],
#                 [175, 3.012, 18849.228],
#                 [110, 5.055, 5486.778],
#                 [98, 0.89, 6069.78],
#                 [86, 5.69, 15720.84],
#                 [86, 1.27, 161000.69],
#                 [65, 0.27, 7260.15],
#                 [63, 0.92, 529.69],
#                 [57, 2.01, 83996.85],
#                 [56, 5.24, 71430.7],
#                 [49, 3.25, 2544.31],
#                 [47, 2.58, 775.52],
#                 [45, 5.54, 9437.76],
#                 [43, 6.01, 6275.96],
#                 [39, 5.36, 4694],
#                 [38, 2.39, 8827.39],
#                 [37, 0.83, 19651.05],
#                 [37, 4.9, 12139.55],
#                 [36, 1.67, 12036.46],
#                 [35, 1.84, 2942.46],
#                 [33, 0.24, 7084.9],
#                 [32, 0.18, 5088.63],
#                 [32, 1.78, 398.15],
#                 [28, 1.21, 6286.6],
#                 [28, 1.9, 6279.55],
#                 [26, 4.59, 10447.39]
#             ],
#             [
#                 [103019, 1.10749, 6283.07585],
#                 [1721, 1.0644, 12566.1517],
#                 [702, 3.142, 0],
#                 [32, 1.02, 18849.23],
#                 [31, 2.84, 5507.55],
#                 [25, 1.32, 5223.69],
#                 [18, 1.42, 1577.34],
#                 [10, 5.91, 10977.08],
#                 [9, 1.42, 6275.96],
#                 [9, 0.27, 5486.78]
#             ],
#             [
#                 [4359, 5.7846, 6283.0758],
#                 [124, 5.579, 12566.152],
#                 [12, 3.14, 0],
#                 [9, 3.63, 77713.77],
#                 [6, 1.87, 5573.14],
#                 [3, 5.47, 18849.23]
#             ],
#             [
#                 [145, 4.273, 6283.076],
#                 [7, 3.92, 12566.15]
#             ],
#             [
#                 [4, 2.56, 6283.08]
#             ]
#         ]
#
#         self._XTable = [
#             [297.85036, 445267.111480, -0.0019142, 1.0 / 189474.0],
#             [357.527772, 35999.050340, -0.0001603, 1.0 / -300000.0],
#             [134.96298, 477198.867398, 0.0086972, 1.0 / 56250.0],
#             [93.27191, 483202.017538, -0.0036825, 1.0 / 327270.0],
#             [125.04452, -1934.136261, 0.0020708, 1.0 / 450000.0]
#         ]
#
#         self._NutationTable = [
#             [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
#             [-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1],
#             [0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5],
#             [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
#             [0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1],
#             [0, 0, 1, 0, 0, 712, 0.1, -7, 0],
#             [-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6],
#             [0, 0, 0, 2, 1, -386, -0.4, 200, 0],
#             [0, 0, 1, 2, 2, -301, 0, 129, -0.1],
#             [-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3],
#             [-2, 0, 1, 0, 0, -158, 0, 0, 0],
#             [-2, 0, 0, 2, 1, 129, 0.1, -70, 0],
#             [0, 0, -1, 2, 2, 123, 0, -53, 0],
#             [2, 0, 0, 0, 0, 63, 0, 0, 0],
#             [0, 0, 1, 0, 1, 63, 0.1, -33, 0],
#             [2, 0, -1, 2, 2, -59, 0, 26, 0],
#             [0, 0, -1, 0, 1, -58, -0.1, 32, 0],
#             [0, 0, 1, 2, 1, -51, 0, 27, 0],
#             [-2, 0, 2, 0, 0, 48, 0, 0, 0],
#             [0, 0, -2, 2, 1, 46, 0, -24, 0],
#             [2, 0, 0, 2, 2, -38, 0, 16, 0],
#             [0, 0, 2, 2, 2, -31, 0, 13, 0],
#             [0, 0, 2, 0, 0, 29, 0, 0, 0],
#             [-2, 0, 1, 2, 2, 29, 0, -12, 0],
#             [0, 0, 0, 2, 0, 26, 0, 0, 0],
#             [-2, 0, 0, 2, 0, -22, 0, 0, 0],
#             [0, 0, -1, 2, 1, 21, 0, -10, 0],
#             [0, 2, 0, 0, 0, 17, -0.1, 0, 0],
#             [2, 0, -1, 0, 1, 16, 0, -8, 0],
#             [-2, 2, 0, 2, 2, -16, 0.1, 7, 0],
#             [0, 1, 0, 0, 1, -15, 0, 9, 0],
#             [-2, 0, 1, 0, 1, -13, 0, 7, 0],
#             [0, -1, 0, 0, 1, -12, 0, 6, 0],
#             [0, 0, 2, -2, 0, 11, 0, 0, 0],
#             [2, 0, -1, 2, 1, -10, 0, 5, 0],
#             [2, 0, 1, 2, 2, -8, 0, 3, 0],
#             [0, 1, 0, 2, 2, 7, 0, -3, 0],
#             [-2, 1, 1, 0, 0, -7, 0, 0, 0],
#             [0, -1, 0, 2, 2, -7, 0, 3, 0],
#             [2, 0, 0, 2, 1, -7, 0, 3, 0],
#             [2, 0, 1, 0, 0, 6, 0, 0, 0],
#             [-2, 0, 2, 2, 2, 6, 0, -3, 0],
#             [-2, 0, 1, 2, 1, 6, 0, -3, 0],
#             [2, 0, -2, 0, 1, -6, 0, 3, 0],
#             [2, 0, 0, 0, 1, -6, 0, 3, 0],
#             [0, -1, 1, 0, 0, 5, 0, 0, 0],
#             [-2, -1, 0, 2, 1, -5, 0, 3, 0],
#             [-2, 0, 0, 0, 1, -5, 0, 3, 0],
#             [0, 0, 2, 2, 1, -5, 0, 3, 0],
#             [-2, 0, 2, 0, 1, 4, 0, 0, 0],
#             [-2, 1, 0, 2, 1, 4, 0, 0, 0],
#             [0, 0, 1, -2, 0, 4, 0, 0, 0],
#             [-1, 0, 1, 0, 0, -4, 0, 0, 0],
#             [-2, 1, 0, 0, 0, -4, 0, 0, 0],
#             [1, 0, 0, 0, 0, -4, 0, 0, 0],
#             [0, 0, 1, 2, 0, 3, 0, 0, 0],
#             [0, 0, -2, 2, 2, -3, 0, 0, 0],
#             [-1, -1, 1, 0, 0, -3, 0, 0, 0],
#             [0, 1, 1, 0, 0, -3, 0, 0, 0],
#             [0, -1, 1, 2, 2, -3, 0, 0, 0],
#             [2, -1, -1, 2, 2, -3, 0, 0, 0],
#             [0, 0, 3, 2, 2, -3, 0, 0, 0],
#             [2, -1, 0, 2, 2, -3, 0, 0, 0]
#         ]
#
#         self._UTable = [
#             84381.448, -4680.93, -1.55, 1999.25,
#             51.38, -249.67, -39.05, 7.12,
#             27.87, 5.79, 2.45
#         ]
#
#     def getHeliocentricLongitude(self, time: JulianDate) -> float:
#         """
#         Computes the heliocentric longitude of the Earth at the given time in radians.
#
#         Args:
#             time: Time to find the longitude.
#
#         Returns:
#             The heliocentric longitude of the Earth in radians.
#         """
#
#         self.__setTimes(time)
#         return self.__expandTableValues(LTable) % TWOPI
#
#     def getHeliocentricLatitude(self, time: JulianDate) -> float:
#         """
#         Computes the heliocentric latitude of the Earth at the given time in radians.
#
#         Args:
#             time: Time to find the latitude.
#
#         Returns:
#             The latitude in radians.
#         """
#
#         self.__setTimes(time)
#         return self.__expandTableValues(BTable) % TWOPI
#
#     def getHeliocentricRadius(self, time: JulianDate) -> float:
#         """
#         Computes the heliocentric radius of the Earth at the given time in astronomical units.
#
#         Args:
#             time: Time to find the distance
#
#         Returns:
#             The distance between the Earth and Sun in astronomical units.
#         """
#
#         self.__setTimes(time)
#         return self.__expandTableValues(RTable)
#
#     def getGeocentricLongitude(self, time: JulianDate) -> float:
#         """
#         Computes the geocentric longitude of the Earth at the given time in radians.
#
#         Args:
#             time: Time to find the longitude.
#
#         Returns:
#             The longitude in radians.
#         """
#
#         self.__setTimes(time)
#         helioLongitude = self.__expandTableValues(LTable) % TWOPI
#         return (helioLongitude + pi) % TWOPI
#
#     def getGeocentricLatitude(self, time: JulianDate) -> float:
#         """
#         Computes the geocentric latitude of the Earth at the given time in radians.
#
#         Args:
#             time: Time to find the longitude.
#
#         Returns:
#             The longitude in radians.
#         """
#         self.__setTimes(time)
#         helioLatitude = self.__expandTableValues(BTable) % TWOPI
#         return -helioLatitude
#
#     def getNutationLongitude(self, time: JulianDate) -> float:
#         self.__setTimes(time)
#         dPsi, UNUSED = self.__getNutationDeltas()
#         return dPsi
#
#     def getNutationObliquity(self, time: JulianDate) -> float:
#         self.__setTimes(time)
#         UNUSED, dEpsilon = self.__getNutationDeltas()
#         return dEpsilon
#
#     def getTrueObliquity(self, time: JulianDate) -> float:
#         """true obliquity of the ecliptic"""
#         self.__setTimes(time)
#         epsilon0 = 0.0
#         for i in range(11):
#             epsilon0 += UTable[i] * ((self._JME / 10.0) ** i)
#         UNUSED, dEpsilon = self.__getNutationDeltas()
#         return radians(epsilon0 / 3600.0) + dEpsilon
#
#     def getAberrationCorrection(self, time: JulianDate, R: float = None) -> float:
#         self.__setTimes(time)
#         if R is None:
#             R = self.__expandTableValues(RTable)
#         return -radians(20.4898 / (3600.0 * R))
#
#     def getApparentSolarLongitude(self, time: JulianDate, R: float = None) -> float:
#         self.__setTimes(time)
#         if R is None:
#             R = self.__expandTableValues(RTable)
#         theta = (self.__expandTableValues(LTable) + pi) % TWOPI
#         dPsi, UNUSED = self.__getNutationDeltas()
#         dTau = -radians(20.4898 / (3600.0 * R))
#         return theta + dPsi + dTau
#
#     # todo: put this in sidereal.py ?
#     def getApparentSiderealTime(self, time: JulianDate) -> float:
#         # in radians
#         # todo: can we convert these to radians and hard code them?
#         self.__setTimes(time)
#         v0 = 280.46061837 + 360.98564736629 * (self._JD.value() - 2451545) + 0.000387933 \
#              * (self._JC ** 2.0) - ((self._JC ** 3.0) / 38710000)
#         dPsi, dEpsilon = self.__getNutationDeltas()
#         epsilon0 = 0.0
#         for i in range(11):
#             epsilon0 += UTable[i] * ((self._JME / 10.0) ** i)
#         epsilon = radians(epsilon0 / 3600.0) + dEpsilon
#         return (radians(v0) % TWOPI) + dPsi * cos(epsilon)
#
#     def getSolarRightAscension(self, time: JulianDate) -> float:
#         # in radians
#         self.__setTimes(time)
#         # todo: optimize this later
#         lmbda = self.getApparentSolarLongitude(time)
#         epsilon = self.getTrueObliquity(time)
#         beta = self.getGeocentricLatitude(time)
#         return atan3(sin(lmbda) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lmbda))
#
#     def getSolarDeclination(self, time: JulianDate) -> float:
#         self.__setTimes(time)
#         # todo: optimize this later
#         lmbda = self.getApparentSolarLongitude(time)
#         epsilon = self.getTrueObliquity(time)
#         beta = self.getGeocentricLatitude(time)
#         return asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lmbda))
#
#     def __setTimes(self, time):
#         if self._JD is None or self._JD.value() != time.value():
#             self._JD = time
#             self._JDE = time.future(self._DELTAT / 86400.0)
#             self._JC = time.difference(J2000) / 36525.0
#             self._JCE = self._JDE.difference(J2000) / 36525.0
#             self._JME = self._JCE / 10.0
#
#     def __expandTableValues(self, table):
#         """Returns the result of the table computations in section 3.2 of the SPA in radians."""
#         arr = [self.__tableSum(ti) for ti in table]
#         sums = 0.0
#         for i in range(len(table)):
#             sums += arr[i] * (self._JME ** i)
#         return sums / 1e8
#
#     def __tableSum(self, table):
#         """Returns table sum of the specified table as laid out in section 3.2 of the SPA in radians."""
#         sums = 0.0
#         for ti in table:
#             sub = ti[0] * cos(ti[1] + ti[2] * self._JME)
#             sums += sub
#         return sums
#
#     def __getNutationDeltas(self):
#         """Returns the delta terms computed from their tables for the nutation terms in section 3.4 of the SPA,
#         in radians."""
#         arrX = self.__getXArray()
#
#         dPsi, dEpsilon = 0.0, 0.0
#         for ti in self._NutationTable:
#             sumTerm = 0.0
#             for i in range(5):
#                 sumTerm += arrX[i] * ti[i]
#             dPsi += (ti[5] + ti[6] * self._JCE) * sin(sumTerm)
#             dEpsilon += (ti[7] + ti[8] * self._JCE) * cos(sumTerm)
#         return radians(dPsi / 3.6e7), radians(dEpsilon / 3.6e7)
#
#     def __getXArray(self):
#         """Computes the table values for the computation of the nutation deltas in section 3.4 of the SPA."""
#         X = []
#         for i in range(5):
#             Xi = 0.0
#             for j in range(4):
#                 Xi += XTable[i][j] * (self._JCE ** j)
#             X.append(radians(Xi))
#         return X


# def test(time: JulianDate, geo: GeoPosition, target = -0.833333333):
#     sc = SunPositionController()
#     # figure the details here later, for now use middle of the day, not noon
#     # A.2.1
#     # jd = JulianDate.fromNumber(time.number())
#     jd = time
#     v = sc.getApparentSiderealTime(jd)
#     # A.2.2
#     dt = DELTAT / 86400.0
#     time_m1 = jd.future(dt - 1)
#     time_0 = jd.future(dt)
#     time_p1 = jd.future(dt + 1)
#     alpha_m1 = sc.getSolarRightAscension(time_m1)
#     alpha_0 = sc.getSolarRightAscension(time_0)
#     alpha_p1 = sc.getSolarRightAscension(time_p1)
#     delta_m1 = sc.getSolarDeclination(time_m1)
#     delta_0 = sc.getSolarDeclination(time_0)
#     delta_p1 = sc.getSolarDeclination(time_p1)
#     # A.2.3
#     sigma = radians(geo.getLongitude())
#     m0 = (alpha_0 - sigma - v) / TWOPI
#     # A.2.4
#     hp0 = radians(target)
#     phi = radians(geo.getLatitude())
#     H0 = acos((sin(hp0) - sin(phi)*sin(delta_0)) / (cos(phi) * cos(delta_0))) % pi
#     # A.2.5
#     m1 = (m0 - (H0 / TWOPI)) % 1.0
#     # A.2.6
#     m2 = (m0 + (H0 / TWOPI)) % 1.0
#     # A.2.8
#     v0 = v + radians(360.985647) * m0
#     v1 = v + radians(360.985647) * m1
#     v2 = v + radians(360.985647) * m2
#     # A.2.9
#     n0 = m0 + dt
#     n1 = m1 + dt
#     n2 = m2 + dt
#     # A.2.10
#     a = degrees(alpha_0 - alpha_m1)
#     if abs(a) > 2:
#         a %= 1.0
#     a = radians(a)
#     b = degrees(alpha_p1 - alpha_0)
#     if abs(b) > 2:
#         b %= 1.0
#     b = radians(b)
#     ap = degrees(delta_0 - delta_m1)
#     if abs(ap) > 2:
#         ap %= 1.0
#     ap = radians(ap)
#     bp = degrees(delta_p1 - delta_0)
#     if abs(bp) > 2:
#         bp %= 1.0
#     bp = radians(bp)
#     c = b - a
#     cp = bp - ap
#     alphap0 = alpha_0 + ((n0 * (a + b + c * n0)) / 2.0)
#     alphap1 = alpha_0 + ((n1 * (a + b + c * n1)) / 2.0)
#     alphap2 = alpha_0 + ((n2 * (a + b + c * n2)) / 2.0)
#     deltap0 = delta_0 + ((n0 * (ap + bp + cp * n0)) / 2.0)
#     deltap1 = delta_0 + ((n1 * (ap + bp + cp * n1)) / 2.0)
#     deltap2 = delta_0 + ((n2 * (ap + bp + cp * n2)) / 2.0)
#     # A.2.11
#     Hp0 = (v0 + sigma - alphap0) % TWOPI
#     if Hp0 >= pi:
#         Hp0 += -TWOPI
#     Hp1 = (v1 + sigma - alphap1) % TWOPI
#     if Hp1 >= pi:
#         Hp1 += -TWOPI
#     Hp2 = (v2 + sigma - alphap2) % TWOPI
#     if Hp2 >= pi:
#         Hp2 += -TWOPI
#     # A.2.12
#     h0 = asin(sin(phi)*sin(deltap0) + cos(phi)*cos(deltap0)*cos(Hp0)) # sun altitude at transit time
#     h1 = asin(sin(phi)*sin(deltap1) + cos(phi)*cos(deltap1)*cos(Hp1))
#     h2 = asin(sin(phi)*sin(deltap2) + cos(phi)*cos(deltap2)*cos(Hp2))
#     # A.2.13
#     T = m0 - (Hp0 / TWOPI)
#     # A.2.14
#     R = m1 + ((h1 - hp0) / (TWOPI*cos(deltap1)*cos(phi)*sin(Hp1)))
#     # A.2.15
#     S = m2 + ((h2 - hp0) / (TWOPI*cos(deltap2)*cos(phi)*sin(Hp2)))
#     return R, T, S


# def getSunAngles(time: JulianDate):
#     JDE = time.future(DELTAT / 86400)
#     # JC = time.difference(J2000) / 36525.0
#     JCE = JDE.difference(J2000) / 36525.0
#     JME = JCE / 10.0
#
#     L = expandTableValues(LTable, JME) % TWOPI
#     B = expandTableValues(BTable, JME) % TWOPI
#     R = expandTableValues(RTable, JME)
#
#     theta = (L + pi) % TWOPI
#     beta = -B
#     dPsi, dEpsilon = getNutationDeltas(JCE)
#     epsilon0 = 0.0
#     for i in range(11):
#         epsilon0 += UTable[i] * ((JME / 10.0) ** i)
#     epsilon = radians(epsilon0 / 3600.0) + dEpsilon
#     dTau = -radians(20.4898 / (3600 * R))
#     lmbda = theta + dPsi + dTau
#     alpha = atan3(sin(lmbda) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lmbda))
#     delta = asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lmbda))
#     return (alpha, delta)


# def __totopos(vec, time, geo):
#     geoVector = geoPositionVector(geo, time)
#     mat = getEulerMatrix(
#         ZYX,
#         EulerAngles(
#             radians(geo.getLongitude()) + earthOffsetAngle(time),
#             radians(90 - geo.getLatitude()),
#             0.0
#         )
#     )
#     return rotateToThenOffset(mat, geoVector, vec)


# def getSunRiseSetTimes(time: JulianDate, geo: GeoPosition):
#     desiredAngle = -0.8333
#     epsilon = 5e-7
#     sunPos = getSunPosition2(time)
#     sunPosSez = __totopos(sunPos, time, geo)
#     sunAlt = degrees(asin(sunPosSez[2] / sunPosSez.mag()))
#
#     # if sun alt < -0.83333 then find NEXT sunrise moving forward and
#         # skip ahead based on the angle value (try and hit 12 the next day) and find next sunset moving forward
#     # else find previous sunrise moving backward and the next sunset moving forward
#     sunSetTime = time
#     sunRiseTime = time
#     if sunAlt < desiredAngle:
#         sunAz = degrees(atan3(sunPosSez[1], -sunPosSez[0]))
#         if sunAz > 180:
#             sunSetTime = sunSetTime.future(0.66)
#             sunRiseTime = sunRiseTime.future(0.66)
#         else:
#             sunSetTime = sunSetTime.future(0.33)
#             sunRiseTime = sunRiseTime.future(0.33)
#
#     # to sunset
#     while True:
#         sunPos = getSunPosition2(sunSetTime)
#         sunPosSez = __totopos(sunPos, sunSetTime, geo)
#         sunAlt = degrees(asin(sunPosSez[2] / sunPosSez.mag()))
#
#         dAngle = sunAlt - desiredAngle
#         if abs(dAngle) < epsilon:
#             break
#
#         dt = dAngle / 360
#         sunSetTime = sunSetTime.future(dt)
#
#     # to sunrise
#     while True:
#         sunPos = getSunPosition2(sunRiseTime)
#         sunPosSez = __totopos(sunPos, sunRiseTime, geo)
#         sunAlt = degrees(asin(sunPosSez[2] / sunPosSez.mag()))
#
#         dAngle = sunAlt - desiredAngle
#         if abs(dAngle) < epsilon:
#             break
#
#         dt = dAngle / 360
#         sunRiseTime = sunRiseTime.future(-dt)
#
#     return sunRiseTime, sunSetTime


# def expandTableValues(table, JME):
#     """returns in radians"""
#     arr = [tableSum(ti, JME) for ti in table]
#     sums = 0.0
#     for i in range(len(table)):
#         sums += arr[i] * (JME ** i)
#     return sums / 1e8
#
#
# def tableSum(table, JME):
#     """returns in radians"""
#     sums = 0.0
#     for ti in table:
#         sub = ti[0] * cos(ti[1] + ti[2] * JME)
#         sums += sub
#     return sums
#
#
# def getNutationDeltas(JCE: float):
#     """returns in radians"""
#     arrX = getXArray(JCE)
#     '''arrDPsi, arrDEpsilon = [], []
#     for row in NutationTable:
#         sumTerm = 0.0:
#         for i in range(5):
#             sumTerm += arrX[i] * row[i]
#     '''
#
#     dPsi, dEpsilon = 0.0, 0.0
#     for ti in NutationTable:
#         sumTerm = 0.0
#         for i in range(5):
#             sumTerm += arrX[i] * ti[i]
#         dPsi += (ti[5] + ti[6] * JCE) * sin(sumTerm)
#         dEpsilon += (ti[7] + ti[8] * JCE) * cos(sumTerm)
#     return radians(dPsi / 36000000.0), radians(dEpsilon / 36000000.0)
#
#
# def getXArray(JCE: float):
#     """returns in radians"""
#     X = []
#     for i in range(5):
#         Xi = 0.0
#         for j in range(4):
#             Xi += XTable[i][j] * (JCE ** j)
#         X.append(radians(Xi))
#     return X
#
#
# LTable = [
#     [
#         [175347046, 0,          0],
#         [3341656,   4.6692568,  6283.07585],
#         [34894,     4.6261,     12566.1517],
#         [3497,      2.7441,     5753.3849],
#         [3418,      2.8289,     3.5231],
#         [3136,      3.6277,     77713.7715],
#         [2676,      4.4181,     7860.4194],
#         [2343,      6.1352,     3930.2097],
#         [1324,      0.7425,     11506.7698],
#         [1273,      2.0371,     529.691],
#         [1199,      1.1096,     1577.3435],
#         [990,       5.233,      5884.927],
#         [902,       2.045,      26.298],
#         [857,       3.508,      398.149],
#         [780,       1.179,      5223.694],
#         [753,       2.533,      5507.553],
#         [505,       4.583,      18849.228],
#         [492,       4.205,      775.523],
#         [357,       2.92,       0.067],
#         [317,       5.849,      11790.629],
#         [284,       1.899,      796.298],
#         [271,       0.315,      10977.079],
#         [243,       0.345,      5486.778],
#         [206,       4.806,      2544.314],
#         [205,       1.869,      5573.143],
#         [202,       2.458,      6069.777],
#         [156,       0.833,      213.299],
#         [132,       3.411,      2942.463],
#         [126,       1.083,      20.775],
#         [115,       0.645,      0.98],
#         [103,       0.636,      4694.003],
#         [102,       0.976,      15720.839],
#         [102,       4.267,      7.114],
#         [99,        6.21,       2146.17],
#         [98,        0.68,       155.42],
#         [86,        5.98,       161000.69],
#         [85,        1.3,        6275.96],
#         [85,        3.67,       71430.7],
#         [80,        1.81,       17260.15],
#         [79,        3.04,       12036.46],
#         [75,        1.76,       5088.63],
#         [74,        3.5,        3154.69],
#         [74,        4.68,       801.82],
#         [70,        0.83,       9437.76],
#         [62,        3.98,       8827.39],
#         [61,        1.82,       7084.9],
#         [57,        2.78,       6286.6],
#         [56,        4.39,       14143.5],
#         [56,        3.47,       6279.55],
#         [52,        0.19,       12139.55],
#         [52,        1.33,       1748.02],
#         [51,        0.28,       5856.48],
#         [49,        0.49,       1194.45],
#         [41,        5.37,       8429.24],
#         [41,        2.4,        19651.05],
#         [39,        6.17,       10447.39],
#         [37,        6.04,       10213.29],
#         [37,        2.57,       1059.38],
#         [36,        1.71,       2352.87],
#         [36,        1.78,       6812.77],
#         [33,        0.59,       17789.85],
#         [30,        0.44,       83996.85],
#         [30,        2.74,       1349.87],
#         [25,        3.16,       4690.48]
#     ],
#     [
#         [628331966747.0, 0,         0],
#         [206059,        2.678235,   6283.07585],
#         [4303,          2.6351,     12566.1517],
#         [425,           1.59,       3.523],
#         [119,           5.769,      26.298],
#         [109,           2.966,      1577.344],
#         [93,            2.59,       18849.23],
#         [72,            1.14,       529.69],
#         [68,            1.87,       398.15],
#         [67,            4.41,       5507.55],
#         [59,            2.89,       5223.69],
#         [56,            2.17,       155.42],
#         [45,            0.4,        796.3],
#         [36,            0.47,       775.52],
#         [29,            2.65,       7.11],
#         [21,            5.34,       0.98],
#         [19,            1.85,       5486.78],
#         [19,            4.97,       213.3],
#         [17,            2.99,       6275.96],
#         [16,            0.03,       2544.31],
#         [16,            1.43,       2146.17],
#         [15,            1.21,       10977.08],
#         [12,            2.83,       1748.02],
#         [12,            3.26,       5088.63],
#         [12,            5.27,       1194.45],
#         [12,            2.08,       4694],
#         [11,            0.77,       553.57],
#         [10,            1.3,        6286.6],
#         [10,            4.24,       1349.87],
#         [9,             2.7,        242.73],
#         [9,             5.64,       951.72],
#         [8,             5.3,        2352.87],
#         [6,             2.65,       9437.76],
#         [6,             4.67,       4690.48]
#     ],
#     [
#         [52919, 0,      0],
#         [8720,  1.0721, 6283.0758],
#         [309,   0.867,  12566.152],
#         [27,    0.05,   3.52],
#         [16,    5.19,   26.3],
#         [16,    3.68,   155.42],
#         [10,    0.76,   18849.23],
#         [9,     2.06,   77713.77],
#         [7,     0.83,   775.52],
#         [5,     4.66,   1577.34],
#         [4,     1.03,   7.11],
#         [4,     3.44,   5573.14],
#         [3,     5.14,   796.3],
#         [3,     6.05,   5507.55],
#         [3,     1.19,   242.73],
#         [3,     6.12,   529.69],
#         [3,     0.31,   398.15],
#         [3,     2.28,   553.57],
#         [2,     4.38,   5223.69],
#         [2,     3.75,   0.98]
#     ],
#     [
#         [289,   5.844,  6283.076],
#         [35,    0,      0],
#         [17,    5.49,   12566.15],
#         [3,     5.2,    155.42],
#         [1,     4.72,   3.52],
#         [1,     5.3,    18849.23],
#         [1,     5.97,   242.73]
#     ],
#     [
#         [114,   3.142,  0],
#         [8,     4.13,   6283.08],
#         [1,     3.84,   12566.15]
#     ],
#     [
#         [1, 3.14, 0]
#     ]
# ]
#
# BTable = [
#     [
#         [280,   3.199,  84334.662],
#         [102,   5.422,  5507.553],
#         [80,    3.88,   5223.69],
#         [44,    3.7,    2352.87],
#         [32,    4,      1577.34]
#     ],
#     [
#         [9, 3.9,    5507.55],
#         [6, 1.73,   5223.69]
#     ]
# ]
#
# RTable = [
#     [
#         [100013989, 0,          0],
#         [1670700,   3.0984635,  6283.07585],
#         [13956,     3.05525,    12566.1517],
#         [3084,      5.1985,     77713.7715],
#         [1628,      1.1739,     5753.3849],
#         [1576,      2.8469,     7860.4194],
#         [925,       5.453,      11506.77],
#         [542,       4.564,      3930.21],
#         [472,       3.661,      5884.927],
#         [346,       0.964,      5507.553],
#         [329,       5.9,        5223.694],
#         [307,       0.299,      5573.143],
#         [243,       4.273,      11790.629],
#         [212,       5.847,      1577.344],
#         [186,       5.022,      10977.079],
#         [175,       3.012,      18849.228],
#         [110,       5.055,      5486.778],
#         [98,        0.89,       6069.78],
#         [86,        5.69,       15720.84],
#         [86,        1.27,       161000.69],
#         [65,        0.27,       7260.15],
#         [63,        0.92,       529.69],
#         [57,        2.01,       83996.85],
#         [56,        5.24,       71430.7],
#         [49,        3.25,       2544.31],
#         [47,        2.58,       775.52],
#         [45,        5.54,       9437.76],
#         [43,        6.01,       6275.96],
#         [39,        5.36,       4694],
#         [38,        2.39,       8827.39],
#         [37,        0.83,       19651.05],
#         [37,        4.9,        12139.55],
#         [36,        1.67,       12036.46],
#         [35,        1.84,       2942.46],
#         [33,        0.24,       7084.9],
#         [32,        0.18,       5088.63],
#         [32,        1.78,       398.15],
#         [28,        1.21,       6286.6],
#         [28,        1.9,        6279.55],
#         [26,        4.59,       10447.39]
#     ],
#     [
#         [103019, 1.10749, 6283.07585],
#         [1721, 1.0644, 12566.1517],
#         [702, 3.142, 0],
#         [32, 1.02, 18849.23],
#         [31, 2.84, 5507.55],
#         [25, 1.32, 5223.69],
#         [18, 1.42, 1577.34],
#         [10, 5.91, 10977.08],
#         [9, 1.42, 6275.96],
#         [9, 0.27, 5486.78]
#     ],
#     [
#         [4359,  5.7846, 6283.0758],
#         [124,   5.579,  12566.152],
#         [12,    3.14,   0],
#         [9,     3.63,   77713.77],
#         [6,     1.87,   5573.14],
#         [3,     5.47,   18849.23]
#     ],
#     [
#         [145,   4.273,  6283.076],
#         [7,     3.92,   12566.15]
#     ],
#     [
#         [4, 2.56, 6283.08]
#     ]
# ]
#
# XTable = [
#     [297.85036, 445267.111480, -0.0019142, 1.0/189474.0],
#     [357.527772, 35999.050340, -0.0001603, 1.0/-300000.0],
#     [134.96298, 477198.867398, 0.0086972, 1.0/56250.0],
#     [93.27191, 483202.017538, -0.0036825, 1.0/327270.0],
#     [125.04452, -1934.136261, 0.0020708, 1.0/450000.0]
# ]
#
# NutationTable = [
#     [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
#     [-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1],
#     [0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5],
#     [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
#     [0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1],
#     [0, 0, 1, 0, 0, 712, 0.1, -7, 0],
#     [-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6],
#     [0, 0, 0, 2, 1, -386, -0.4, 200, 0],
#     [0, 0, 1, 2, 2, -301, 0, 129, -0.1],
#     [-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3],
#     [-2, 0, 1, 0, 0, -158, 0, 0, 0],
#     [-2, 0, 0, 2, 1, 129, 0.1, -70, 0],
#     [0, 0, -1, 2, 2, 123, 0, -53, 0],
#     [2, 0, 0, 0, 0, 63, 0, 0, 0],
#     [0, 0, 1, 0, 1, 63, 0.1, -33, 0],
#     [2, 0, -1, 2, 2, -59, 0, 26, 0],
#     [0, 0, -1, 0, 1, -58, -0.1, 32, 0],
#     [0, 0, 1, 2, 1, -51, 0, 27, 0],
#     [-2, 0, 2, 0, 0, 48, 0, 0, 0],
#     [0, 0, -2, 2, 1, 46, 0, -24, 0],
#     [2, 0, 0, 2, 2, -38, 0, 16, 0],
#     [0, 0, 2, 2, 2, -31, 0, 13, 0],
#     [0, 0, 2, 0, 0, 29, 0, 0, 0],
#     [-2, 0, 1, 2, 2, 29, 0, -12, 0],
#     [0, 0, 0, 2, 0, 26, 0, 0, 0],
#     [-2, 0, 0, 2, 0, -22, 0, 0, 0],
#     [0, 0, -1, 2, 1, 21, 0, -10, 0],
#     [0, 2, 0, 0, 0, 17, -0.1, 0, 0],
#     [2, 0, -1, 0, 1, 16, 0, -8, 0],
#     [-2, 2, 0, 2, 2, -16, 0.1, 7, 0],
#     [0, 1, 0, 0, 1, -15, 0, 9, 0],
#     [-2, 0, 1, 0, 1, -13, 0, 7, 0],
#     [0, -1, 0, 0, 1, -12, 0, 6, 0],
#     [0, 0, 2, -2, 0, 11, 0, 0, 0],
#     [2, 0, -1, 2, 1, -10, 0, 5, 0],
#     [2, 0, 1, 2, 2, -8, 0, 3, 0],
#     [0, 1, 0, 2, 2, 7, 0, -3, 0],
#     [-2, 1, 1, 0, 0, -7, 0, 0, 0],
#     [0, -1, 0, 2, 2, -7, 0, 3, 0],
#     [2, 0, 0, 2, 1, -7, 0, 3, 0],
#     [2, 0, 1, 0, 0, 6, 0, 0, 0],
#     [-2, 0, 2, 2, 2, 6, 0, -3, 0],
#     [-2, 0, 1, 2, 1, 6, 0, -3, 0],
#     [2, 0, -2, 0, 1, -6, 0, 3, 0],
#     [2, 0, 0, 0, 1, -6, 0, 3, 0],
#     [0, -1, 1, 0, 0, 5, 0, 0, 0],
#     [-2, -1, 0, 2, 1, -5, 0, 3, 0],
#     [-2, 0, 0, 0, 1, -5, 0, 3, 0],
#     [0, 0, 2, 2, 1, -5, 0, 3, 0],
#     [-2, 0, 2, 0, 1, 4, 0, 0, 0],
#     [-2, 1, 0, 2, 1, 4, 0, 0, 0],
#     [0, 0, 1, -2, 0, 4, 0, 0, 0],
#     [-1, 0, 1, 0, 0, -4, 0, 0, 0],
#     [-2, 1, 0, 0, 0, -4, 0, 0, 0],
#     [1, 0, 0, 0, 0, -4, 0, 0, 0],
#     [0, 0, 1, 2, 0, 3, 0, 0, 0],
#     [0, 0, -2, 2, 2, -3, 0, 0, 0],
#     [-1, -1, 1, 0, 0, -3, 0, 0, 0],
#     [0, 1, 1, 0, 0, -3, 0, 0, 0],
#     [0, -1, 1, 2, 2, -3, 0, 0, 0],
#     [2, -1, -1, 2, 2, -3, 0, 0, 0],
#     [0, 0, 3, 2, 2, -3, 0, 0, 0],
#     [2, -1, 0, 2, 2, -3, 0, 0, 0]
# ]
#
# UTable = [
#     84381.448, -4680.93, -1.55, 1999.25,
#     51.38, -249.67, -39.05, 7.12,
#     27.87, 5.79, 2.45
# ]
