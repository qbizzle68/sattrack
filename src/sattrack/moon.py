from pyevspace import EVector
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.coordinates import CelestialCoordinates


def getMoonPosition(jd: JulianDate) -> EVector:
    pass


def getMoonTimes(jd: JulianDate) -> tuple[JulianDate]:
    pass


def getMoonCelestialCoordinates(jd: JulianDate, topocentric=True) -> CelestialCoordinates:
    pass


def getMoonHourAngle(jd: JulianDate, topocentric=True) -> CelestialCoordinates:
    pass