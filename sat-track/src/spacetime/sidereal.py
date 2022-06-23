from spacetime.juliandate import JulianDate, J2000


def siderealTime(jd: JulianDate) -> float:
    """Converts a certain solar time into sidereal time, most useful for computing Earth's
    geographic reference frame offset from the geocentric equitorial celestial reference frame.
    I.e. this is the offset angle of GMT to the first point of Aries measured in hours.
    Parameters:
    jd: The time to convert.
    return: The local sidereal time of Greenwich, England."""
    dt = jd.difference(J2000)
    tmp = (18.697_374_558 + 24.065_709_824_419_08 * dt) % 24.0
    return tmp + 24.0 if tmp < 0 else tmp


def localSiderealTime(jd: JulianDate, lng: float) -> float:
    """Converts a certain solar time into local sidereal time, useful for computing your
    longitudinal offset from the geocentric equitorial celestial reference frame.
    Parameters:
    jd:     The time to convert.
    lng:    The local longitude.
    Returns the local sidereal time in hours."""
    if lng < 0:
        lng += 360.0
    lng2RA = lng / 360.0 * 24.0
    lst = siderealTime(jd) + lng2RA
    return lst % 24.0


def earthOffsetAngle(jd: JulianDate) -> float:
    """Computes the earth offset angle from the celestial reference frame for a
    given solar time.
    Parameters:
    jd: The time to compute.
    Returns the offset angle in degrees."""
    return siderealTime(jd) / 24.0 * 360.0
