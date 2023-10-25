from math import cos

from pyevspace import Vector

from sattrack.bodies.body import BodyOrbitController, Body
from sattrack.bodies.position import JulianTimes, computeNutationDeltas, computeMeanObliquity
from sattrack.util.constants import TWOPI, EARTH_MU, EARTH_EQUITORIAL_RADIUS, SIDEREAL_PER_SOLAR

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.coordinates import GeoPosition
    from sattrack.core.juliandate import JulianDate


def computeMeanSiderealTime(time: 'JulianDate | JulianTimes') -> float:
    if isinstance(time, JulianTimes):
        JD = time.JD
        JC = time.JC
    else:
        JD = time.value
        JC = (JD - 2451545) / 36525

    rtn = 4.894961212735793 + 6.300388098984957 * (JD - 2451545) \
        + JC * JC * (6.770708127139162e-06 - JC * 4.508729661571505e-10)

    return rtn % TWOPI


def computeApparentSiderealTime(time: 'JulianDate | JulianTimes') -> float:
    deltaPsi, deltaEpsilon = computeNutationDeltas(time)
    meanObliquity = computeMeanObliquity(time)
    trueObliquity = meanObliquity + deltaEpsilon
    meanSiderealTime = computeMeanSiderealTime(time)

    return (meanSiderealTime + deltaPsi * cos(trueObliquity)) % TWOPI


class EarthController(BodyOrbitController):
    def __init__(self):
        pass

    def computePosition(self, time: 'JulianDate') -> Vector:
        return Vector(0, 0, 0)

    def computeCelestialCoordinates(self, time: 'JulianDate'):
        pass

    def computeAltAz(self, geo: 'GeoPosition', time: 'JulianDate'):
        pass

    @staticmethod
    def computeMeanSiderealTime(time: 'JulianDate') -> float:
        return computeMeanSiderealTime(time)

    @staticmethod
    def computeApparentSiderealTime(time: 'JulianDate') -> float:
        return computeApparentSiderealTime(time)

    @staticmethod
    def computeLocalSiderealTime(time: 'JulianDate', longitude: float) -> float:
        # longitude in degrees

        siderealTime = computeApparentSiderealTime(time)

        return (siderealTime + longitude) % TWOPI

    @staticmethod
    def offsetAngle(time: 'JulianDate') -> float:
        return computeApparentSiderealTime(time)


Earth = Body('Earth', EARTH_MU, EARTH_EQUITORIAL_RADIUS, SIDEREAL_PER_SOLAR, EarthController)
