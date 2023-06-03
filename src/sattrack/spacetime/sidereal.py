from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import TWOPI

__all__ = ['siderealTime', 'localSiderealTime', 'earthOffsetAngle']


def siderealTime(jd: JulianDate) -> float:
    """Converts a solar time (via JulianDate) into sidereal time in hours."""

    dt = jd - J2000
    tmp = (18.697_374_558 + 24.065_709_824_419_08 * dt) % 24.0
    return tmp + 24.0 if tmp < 0 else tmp


def localSiderealTime(jd: JulianDate, lng: float) -> float:
    """Converts a solar time (via JulianDate) to a local sidereal time in hours, defined by a longitude in degrees."""

    if lng < 0:
        lng += 360.0
    lng2RA = lng / 360.0 * 24.0
    lst = siderealTime(jd) + lng2RA
    return lst % 24.0


def earthOffsetAngle(jd: JulianDate) -> float:
    """Computes the earth offset angle from the celestial reference frame in radians."""

    return siderealTime(jd) / 24.0 * TWOPI
