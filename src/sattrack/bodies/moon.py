from pyevspace import Vector
from sattrack.core.juliandate import JulianDate
from sattrack.core.coordinates import CelestialCoordinates

__all__ = []


def getMoonPosition(jd: JulianDate) -> Vector:
    pass


def getMoonTimes(jd: JulianDate) -> tuple[JulianDate]:
    pass


def getMoonCelestialCoordinates(jd: JulianDate, topocentric=True) -> CelestialCoordinates:
    pass


def getMoonHourAngle(jd: JulianDate, topocentric=True) -> CelestialCoordinates:
    pass