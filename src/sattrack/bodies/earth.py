from math import radians

from pyevspace import Vector

from sattrack.bodies.body import BodyOrbitController, Body
from sattrack.bodies.position import computeMeanSiderealTime, computeApparentSiderealTime
from sattrack.bodies.sun import Sun
from sattrack.util.constants import TWOPI, EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS, \
    EARTH_SIDEREAL_PERIOD

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.coordinates import GeoPosition
    from sattrack.core.juliandate import JulianDate


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
        lngRadians = radians(longitude)

        return (siderealTime + lngRadians) % TWOPI


Earth = Body('Earth', EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_SIDEREAL_PERIOD, EarthController,
             Rp=EARTH_POLAR_RADIUS, parent=Sun)
