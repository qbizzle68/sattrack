from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import TWOPI

__all__ = ('siderealTime', 'localSiderealTime', 'earthOffsetAngle')


def siderealTime(jd: JulianDate) -> float:
    """
    Converts a solar time into sidereal time, most useful for computing the right-ascension of earth's zero longitude,
    equal to the geographic reference frame offset from the geocentric equitorial celestial reference frame.

    Args:
        jd: Solar time to convert to sidereal time.

    Returns:
        The Earth's sidereal time at the given solar time.
    """

    dt = jd - J2000
    tmp = (18.697_374_558 + 24.065_709_824_419_08 * dt) % 24.0
    return tmp + 24.0 if tmp < 0 else tmp


def localSiderealTime(jd: JulianDate, lng: float) -> float:
    """
    Computes the sidereal time for a line of longitude on Earth, which is equivalent to the offset angle between a line
    of longitude and the first point of Ares, or the celestial x-axis.

    Args:
        jd: Solar time to convert to sidereal time.
        lng: Longitude of the GeoPosition whose local sidereal time is being found.

    Returns:
        The local sidereal time of a longitude at the given solar time.
    """

    if lng < 0:
        lng += 360.0
    lng2RA = lng / 360.0 * 24.0
    lst = siderealTime(jd) + lng2RA
    return lst % 24.0


def earthOffsetAngle(jd: JulianDate) -> float:
    """Computes the earth offset angle from the celestial reference frame in radians."""
    return siderealTime(jd) / 24.0 * TWOPI
