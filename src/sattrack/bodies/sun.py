from enum import Enum
from math import sin, cos, tan, asin, degrees, pi, acos

from pyevspace import Vector

from sattrack.bodies._tables import SUN_L_TABLE, SUN_B_TABLE, SUN_R_TABLE
from sattrack.bodies.body import BodyOrbitController, Body
from sattrack.bodies.exceptions import SunRiseSetException
from sattrack.bodies.position import computeNutationDeltas, \
    computeTrueObliquity, JulianTimes, computeMeanObliquity, computeApparentSiderealTime
from sattrack.bodies.topocentric import toTopocentricOffset
from sattrack.core.coordinates import CelestialCoordinates, AltAz
from sattrack.core.juliandate import JulianDate
from sattrack.util.constants import AU, RAD_TO_HOURS, TWOPI, DELTAT, SUN_MU, SUN_RADIUS
from sattrack.util.helpers import atan3

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.coordinates import GeoPosition


def _computeFromTable(table: list[list[float]], time: JulianDate | JulianTimes) -> list[float]:
    if isinstance(time, JulianTimes):
        JME = time.JME
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525
        JME = JCE / 10

    output = []
    for subTable in table:
        Ln = 0
        for (a, b, c) in subTable:
            Ln += a * cos(b + c * JME)
        output.append(Ln)

    return output


def computeEarthHeliocentricLongitude(time: JulianDate | JulianTimes) -> float:
    # Returns the Earth heliocentric longitude in radians.

    if isinstance(time, JulianTimes):
        JME = time.JME
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525
        JME = JCE / 10

    (L0, L1, L2, L3, L4, L5) = _computeFromTable(SUN_L_TABLE, time)
    rtn = L0 + JME * (L1 + JME * (L2 + JME * (L3 + JME * (L4 + JME * L5))))

    return (rtn / 1e8) % TWOPI


def computeEarthHeliocentricLatitude(time: JulianDate | JulianTimes) -> float:
    # Returns the Earth heliocentric latitude in radians.

    if isinstance(time, JulianTimes):
        JME = time.JME
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525
        JME = JCE / 10

    (B0, B1) = _computeFromTable(SUN_B_TABLE, time)
    rtn = ((B0 + JME * B1) / 1e8) % TWOPI
    if rtn > pi:
        rtn -= TWOPI

    return rtn


def computeSunDistance(time: JulianDate | JulianTimes) -> float:
    # Returns the Earth radius in Astronomical Units.

    if isinstance(time, JulianTimes):
        JME = time.JME
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525
        JME = JCE / 10

    (R0, R1, R2, R3, R4) = _computeFromTable(SUN_R_TABLE, time)
    rtn = R0 + JME * (R1 + JME * (R2 + JME * (R3 + JME * R4)))

    return rtn / 1e8


def computeSunGeocentricLongitude(time: JulianDate) -> float:
    # Returns the Sun's geocentric longitude in radians.

    L = computeEarthHeliocentricLongitude(time)

    return (L + pi) % TWOPI


def computeSunGeocentricLatitude(time: JulianDate) -> float:
    # Returns the Sun's geocentric latitude in radians.

    B = computeEarthHeliocentricLatitude(time)

    # computeEarthHeliocentricLatitude modulates output to [0-TWOPI), so we don't need to.
    return -B


# Don't think this will be needed publicly, but can be changed to.
def _computeAberrationCorrection(arg: JulianDate | float) -> float:
    """Computes the aberration correction deltaTau in radians. If arg is a JulianDate instance, Earth
    radius will be computed, otherwise arg should be the Earth's radius in AU."""

    if isinstance(arg, (JulianDate, JulianTimes)):
        R = computeSunDistance(arg)
    else:
        R = arg

    return -9.93373536319817e-05 / R


def computeSunApparentLongitude(time: JulianDate) -> float:
    # Compute the apparent solar longitude in radians.

    solarLongitude = computeSunGeocentricLongitude(time)
    nutationLongitude, _ = computeNutationDeltas(time)
    aberrationCorrection = _computeAberrationCorrection(time)

    return solarLongitude + nutationLongitude + aberrationCorrection


def computeSunRightAscension(time: JulianDate) -> float:
    # Compute the geocentric sun right ascension in radians.

    apparentLongitude = computeSunApparentLongitude(time)
    trueObliquity = computeTrueObliquity(time)
    # Save overhead of another function call by computing heliocentric latitude from geocentric.
    sunLatitude = -computeEarthHeliocentricLatitude(time)

    numerator = sin(apparentLongitude) * cos(trueObliquity) - tan(sunLatitude) * sin(trueObliquity)
    return atan3(numerator, cos(apparentLongitude))


def computeSunDeclination(time: JulianDate) -> float:
    # Compute the geocentric sun declination in radians.

    # Save overhead of another function call by computing heliocentric latitude from geocentric.
    sunLatitude = -computeEarthHeliocentricLatitude(time)
    trueObliquity = computeTrueObliquity(time)
    apparentLongitude = computeSunApparentLongitude(time)

    term1 = sin(sunLatitude) * cos(trueObliquity)
    term2 = cos(sunLatitude) * sin(trueObliquity) * sin(apparentLongitude)

    return asin(term1 + term2)


def _computeSunCoordinatesFast(time: JulianTimes, R: float = None) -> CelestialCoordinates:
    # Compute solar celestial coordinates directly, saving some function calling overhead.
    # If R is not None, it should be the Earth Sun distance in AU to pass to
    # _computeAberrationCorrection()

    solarLongitude = (computeEarthHeliocentricLongitude(time) + pi) % TWOPI
    nutationLongitude, nutationObliquity = computeNutationDeltas(time)
    if R is None:
        aberrationCorrection = _computeAberrationCorrection(time)
    else:
        aberrationCorrection = _computeAberrationCorrection(R)

    apparentLongitude = solarLongitude + nutationLongitude + aberrationCorrection

    meanObliquity = computeMeanObliquity(time)
    trueObliquity = meanObliquity + nutationObliquity

    sunLatitude = -computeEarthHeliocentricLatitude(time)

    numerator = sin(apparentLongitude) * cos(trueObliquity) - tan(sunLatitude) * sin(trueObliquity)
    # numerator = sin(solarLongitude) * cos(trueObliquity) - tan(sunLatitude) * sin(trueObliquity)
    rightAscension = atan3(numerator, cos(apparentLongitude))
    # rightAscension = atan3(numerator, cos(solarLongitude))

    term1 = sin(sunLatitude) * cos(trueObliquity)
    term2 = cos(sunLatitude) * sin(trueObliquity) * sin(apparentLongitude)
    # term2 = cos(sunLatitude) * sin(trueObliquity) * sin(solarLongitude)

    declination = asin(term1 + term2)

    return rightAscension, declination


def computeSunPosition(time: JulianDate) -> Vector:
    # Compute the Sun's geocentric position vector:

    rightAscension, declination = _computeSunCoordinatesFast(time)
    sunDistance = computeSunDistance(time)

    zTerm = sunDistance * sin(declination)
    topoProjection = sunDistance * cos(declination)
    xTerm = topoProjection * cos(rightAscension)
    yTerm = topoProjection * sin(rightAscension)

    return Vector(xTerm, yTerm, zTerm) * AU


def computeSunCoordinates(time: JulianDate) -> CelestialCoordinates:
    # Compute the Sun's celestial coordinates.

    rightAscension, declination = _computeSunCoordinatesFast(time)

    raHours = rightAscension * RAD_TO_HOURS
    decDegrees = degrees(declination)

    return CelestialCoordinates(raHours, decDegrees)


# todo: this can be generalized to compute rise/set times of any object
def _computeSunAngleData(geo: 'GeoPosition', time: JulianDate, target: float) \
        -> (JulianDate, JulianDate, JulianDate, float):
    # Returns rise time, set time, transit time and transit altitude.

    tzOffset = time.timezone / 24.0
    localVal = time.value + tzOffset
    num = int(localVal)
    if localVal - int(localVal) < 0.5:
        num -= 1
    ut = JulianDate.fromNumber(num + 0.5 - tzOffset, time.timezone)

    apparentSiderealTime = computeApparentSiderealTime(ut)

    dt = DELTAT / 86400
    time_m1 = JulianTimes(ut.future(dt - 1))
    time_0 = JulianTimes(ut.future(dt))
    time_p1 = JulianTimes(ut.future(dt + 1))

    alpha_m1, delta_m1 = _computeSunCoordinatesFast(time_m1)
    alpha_0, delta_0 = _computeSunCoordinatesFast(time_0)
    alpha_p1, delta_p1 = _computeSunCoordinatesFast(time_p1)

    m0 = ((alpha_0 - geo.longitudeRadians - apparentSiderealTime) / TWOPI) % 1.0

    try:
        numerator = sin(target) - sin(geo.latitudeRadians) * sin(delta_0)
        H0 = acos(numerator / (cos(geo.latitudeRadians) * cos(delta_0))) % pi
    except ValueError:
        # todo: check the exception message here
        raise SunRiseSetException('sun does not rise or set')

    m1 = (m0 - H0 / TWOPI) % 1.0
    m2 = (m0 + H0 / TWOPI) % 1.0

    v0 = apparentSiderealTime + 6.300388092591991 * m0
    v1 = apparentSiderealTime + 6.300388092591991 * m1
    v2 = apparentSiderealTime + 6.300388092591991 * m2

    n0 = m0 + dt
    n1 = m1 + dt
    n2 = m2 + dt

    maxValue = 0.03490658503988659  # radians(2.0)
    modulator = 0.017453292519943295  # radians(1.0)
    a = alpha_0 - alpha_m1
    if abs(a) > maxValue:
        a %= modulator
    aPrime = delta_0 - delta_m1
    if abs(aPrime) > maxValue:
        aPrime %= modulator
    b = alpha_p1 - alpha_0
    if abs(b) > maxValue:
        b %= modulator
    bPrime = delta_p1 - delta_0
    if abs(bPrime) > maxValue:
        bPrime %= modulator
    c = b - a
    cPrime = bPrime - aPrime

    alphaPrime0 = alpha_0 + (n0 * (a + b + c * n0) / 2)
    alphaPrime1 = alpha_0 + (n1 * (a + b + c * n1) / 2)
    alphaPrime2 = alpha_0 + (n2 * (a + b + c * n2) / 2)
    deltaPrime0 = delta_0 + (n0 * (aPrime + bPrime + cPrime * n0) / 2)
    deltaPrime1 = delta_0 + (n1 * (aPrime + bPrime + cPrime * n1) / 2)
    deltaPrime2 = delta_0 + (n2 * (aPrime + bPrime + cPrime * n2) / 2)

    HPrime0 = (v0 + geo.longitudeRadians - alphaPrime0) % TWOPI
    if HPrime0 <= -pi:
        HPrime0 += TWOPI
    elif HPrime0 >= pi:
        HPrime0 -= TWOPI
    HPrime1 = (v1 + geo.longitudeRadians - alphaPrime1) % TWOPI
    if HPrime1 <= -pi:
        HPrime1 += TWOPI
    elif HPrime1 >= pi:
        HPrime1 -= TWOPI
    HPrime2 = (v2 + geo.longitudeRadians - alphaPrime2) % TWOPI
    if HPrime2 <= -pi:
        HPrime2 += TWOPI
    elif HPrime2 >= pi:
        HPrime2 -= TWOPI

    h0 = asin(sin(geo.latitudeRadians) * sin(deltaPrime0)
              + cos(geo.latitudeRadians) * cos(deltaPrime0) * cos(HPrime0))
    h1 = asin(sin(geo.latitudeRadians) * sin(deltaPrime1)
              + cos(geo.latitudeRadians) * cos(deltaPrime1) * cos(HPrime1))
    h2 = asin(sin(geo.latitudeRadians) * sin(deltaPrime2)
              + cos(geo.latitudeRadians) * cos(deltaPrime2) * cos(HPrime2))

    T = m0 - HPrime0 / TWOPI

    denominator = TWOPI * cos(deltaPrime1) * cos(geo.latitudeRadians) * sin(HPrime1)
    R = m1 + (h1 - target) / denominator
    denominator = TWOPI * cos(deltaPrime2) * cos(geo.latitudeRadians) * sin(HPrime2)
    S = m2 + (h2 - target) / denominator

    return ut.future(R), ut.future(S), ut.future(T), h0


# fixme: name this better
# Public interface to _computeSunAngleData if wanting all four values.
def computeSunData(geo: 'GeoPosition', time: JulianDate) -> (JulianDate, JulianDate, JulianDate, float):
    target = -0.01454441043328608
    data = _computeSunAngleData(geo, time, target)

    return data


def computeSunRiseSetTimes(geo: 'GeoPosition', time: JulianDate) -> (JulianDate, JulianDate):
    # target is -0.833333333333 degrees in radians.
    target = -0.01454441043328608
    data = _computeSunAngleData(geo, time, target)

    return data[:2]


def computeSunTransitInfo(geo: 'GeoPosition', time: JulianDate) -> (JulianDate, float):
    data = _computeSunAngleData(geo, time, 0.0)

    return data[2:]


def computeSunTwilightTimes(geo: 'GeoPosition', time: JulianDate, twilightType: 'Twilight') -> (JulianDate, JulianDate):
    if twilightType == Twilight.Civil:
        data = _computeSunAngleData(geo, time, -0.10471975511965978)
    elif twilightType == Twilight.Nautical:
        data = _computeSunAngleData(geo, time, -0.20943951023931956)
    elif twilightType == Twilight.Astronomical:
        data = _computeSunAngleData(geo, time, -0.3141592653589793)
    else:
        data = _computeSunAngleData(geo, time, -0.01454441043328608)

    return data[:2]


class Twilight(Enum):
    Day = 0
    Civil = 1
    Nautical = 2
    Astronomical = 3
    Night = 4


def computeTwilightType(geo: 'GeoPosition', time: JulianDate) -> Twilight:
    sunPosition = computeSunPosition(time)
    topoSunPosition = toTopocentricOffset(sunPosition, geo, time)
    sunAngle = asin(topoSunPosition[2] / topoSunPosition.mag())

    if sunAngle < -0.3141592653589793:  # 18 degrees
        return Twilight.Night
    elif sunAngle < -0.20943951023931956:   # 12 degrees
        return Twilight.Astronomical
    elif sunAngle < -0.10471975511965978:   # 6 degrees
        return Twilight.Nautical
    elif sunAngle < -0.01454441043328608:    # 50 arc-minutes
        return Twilight.Civil
    else:
        return Twilight.Day


class SunController(BodyOrbitController):

    def __init__(self):
        pass

    def computePosition(self, time: 'JulianDate') -> Vector:
        timeEphem = JulianTimes(time)
        sunDistance = computeSunDistance(timeEphem)
        rightAscension, declination = _computeSunCoordinatesFast(timeEphem, sunDistance)

        zTerm = sunDistance * sin(declination)
        topoProjection = sunDistance * cos(declination)
        xTerm = topoProjection * cos(rightAscension)
        yTerm = topoProjection * sin(rightAscension)

        return Vector(xTerm, yTerm, zTerm) * AU

    def computeCelestialCoordinates(self, time: 'JulianDate') -> CelestialCoordinates:
        timeEphem = JulianTimes(time)
        rightAscension, declination = _computeSunCoordinatesFast(timeEphem)

        raHours = rightAscension * RAD_TO_HOURS
        decDegrees = degrees(declination)

        return CelestialCoordinates(raHours, decDegrees)

    def computeAltAz(self, geo: 'GeoPosition', time: 'JulianDate') -> 'AltAz':
        # todo: compute this correctly
        position = self.computePosition(time)
        topoPosition = toTopocentricOffset(position, geo, time)

        alt = degrees(asin(topoPosition[2] / topoPosition.mag()))
        az = degrees(atan3(topoPosition[1], -topoPosition[0]))

        return AltAz(alt, az)

    @staticmethod
    def computeRiseSetTimes(geo: 'GeoPosition', time: 'JulianDate') -> (JulianDate, JulianDate):
        return computeSunRiseSetTimes(geo, time)

    @staticmethod
    def computeTransitInfo(geo: 'GeoPosition', time: 'JulianDate') -> (JulianDate, float):
        return computeSunTransitInfo(geo, time)

    @staticmethod
    def computeTwilightTimes(geo: 'GeoPosition', time: 'JulianDate', twilightType: 'Twilight'):
        return computeSunTwilightTimes(geo, time, twilightType)

    @staticmethod
    def computeTwilightType(geo: 'GeoPosition', time: 'JulianDate') -> Twilight:
        return computeTwilightType(geo, time)


Sun = Body('Sun', SUN_MU, SUN_RADIUS, 0.0, SunController)
