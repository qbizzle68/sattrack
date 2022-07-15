from math import sqrt, sin, cos, floor

from sattrack.spacetime.juliandate import JulianDate
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import atan2


def meanToTrue(meanAnomaly: float, ecc: float) -> float:
    """Converts a mean anomaly to its true anomaly.
    Parameters:
    meanAnomaly:    The mean anomaly measured in radians.
    ecc:            Eccentricity of the orbit.
    returns:        The true anomaly measured in radians."""

    eAnom = __m2ENewtonRaphson(meanAnomaly, ecc) % TWOPI
    # eAnom = eTemp if eTemp >= 0 else eTemp + 360.0
    y = sqrt(1 - (ecc * ecc)) * sin(eAnom)
    return atan2(y, cos(eAnom) - ecc)


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly to its eccentric anomaly.
    Parameters:
    meanAnomaly:    The mean anomaly measured in radians.
    ecc:            Eccentricity of the orbit.
    returns:        The eccentric anomaly measured in radians."""

    return __m2ENewtonRaphson(meanAnomaly, eccentricity) % TWOPI
    # e = __m2ENewtonRaphson(meanAnomaly, eccentricity) % 360.0
    # return e if e >= 0 else e + 360.0


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    """Converts an eccentric anomaly to its true anomaly.
    Parameters:
    eccAnom:    The eccentric anomaly measured in radians.
    ecc:        Eccentricity of the orbit.
    returns:    The true anomaly measured in radians."""

    y = sqrt(1 - (ecc * ecc)) * sin(eccAnom)
    return atan2(y, cos(eccAnom) - ecc)


def trueToMean(trueAnom: float, ecc: float) -> float:
    """Converts a true anomaly to its mean anomaly.
    Parameters:
    trueAnom:   The true anomaly measured in radians.
    ecc:        Eccentricity of the orbit.
    returns:    The mean anomaly measured in radians."""

    y = sqrt(1 - (ecc * ecc)) * sin(trueAnom)
    eAnom = atan2(y, cos(trueAnom) + ecc)
    return eAnom - ecc * sin(eAnom) % TWOPI
    # return m if m >= 0 else m + 360.0


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    """Converts a true anomaly to its eccentric anomaly.
    Parameters:
    trueAnom:   The true anomaly measured in radians.
    ecc:        Eccentricity of the orbit.
    returns:    The eccentric anomaly measured in radians."""

    y = sqrt(1 - (ecc * ecc)) * sin(trueAnom)
    return atan2(y, cos(trueAnom) + ecc)


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    """Converts an eccentric anomaly to its mean anomaly.
    Parameters:
    eccAnom:    The eccentric anomaly measured in radians.
    ecc:        Eccentricity of the orbit.
    returns:    The mean anomaly measured in radians."""

    return (eccAnom - ecc * sin(eccAnom)) % TWOPI
    # return m if m >= 0 else m + 360.0


def __m2ENewtonRaphson(M: float, ecc: float) -> float:
    """Computes the eccentric anomaly from a given mean anomaly via Kepler's equation.
    The computation is done using the Newton-Raphson method to a hard-coded epsilon value
    of 1e-7. Inputs and outputs are degrees, the conversion is handed in the implementation.
    Parameters:
    M:      The mean anomaly measured in radians.
    ecc:    Eccentricity of the orbit."""

    # M = radians(M)
    Ej = M
    while True:
        num = Ej - (ecc * sin(Ej)) - M
        den = 1 - ecc * cos(Ej)
        Ej1 = Ej - (num / den)
        if abs(Ej1 - Ej) <= 1e-7:
            return Ej1
        else:
            Ej = Ej1


def timeToNextMeanAnomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    """Computes the next time a satellite passes through the mean anomaly after the given time.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        m0 -- a known mean anomaly at epoch0 in radians
        epoch0 -- the epoch the satellite has a mean anomaly of m0
        m1 -- mean anomaly in radians
        time -- relative time to find the next anomaly
    Returns the next time the satellite's position is m1.
    """

    n = meanMotion * 86400 / TWOPI
    pePass = epoch0.future(-m0 / (TWOPI * n))
    revs = time.difference(pePass) * n
    m0 = (revs - floor(revs)) * TWOPI
    if m1 < m0:
        dm = m1 + TWOPI - m0
    else:
        dm = m1 - m0
    return time.future((dm / n) / TWOPI)


def timeToPrevMeanAnomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    """Computes the previous time the satellite passed through the mean anomaly before the given time.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        m0 -- a known mean anomaly at epoch0 in radians
        epoch0 -- the epoch the satellite has a mean anomaly of m0
        m1 -- mean anomaly in radians
        time -- relative time to find the previous anomaly
    Returns the previous time the satellite's position was m1.s
    """

    n = meanMotion * 86400 / TWOPI
    pePass = epoch0.future(-m0 / (TWOPI * n))
    revs = time.difference(pePass) * n
    m0 = (revs - floor(revs)) * TWOPI
    if m0 < m1:
        dm = m1 - TWOPI - m0
    else:
        dm = m1 - m0
    return time.future((dm / n) / TWOPI)


def timeToNextTrueAnomaly(meanMotion: float, ecc: float, t0: float, epoch0: JulianDate, t1: float,
                          time: JulianDate) -> JulianDate:
    """Computes the next time the satellite passes through the true anomaly after the time given.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        ecc -- eccentricity of the orbit
        t0 -- a known true anomaly at epoch0 in radians
        epoch0 -- the epoch the satellite has a true anomaly of t0
        t1 -- true anomaly in radians
        time -- relative time to find the next anomaly
    Returns the next time the satellites position is t1.
    """

    return timeToNextMeanAnomaly(meanMotion, trueToMean(t0, ecc), epoch0, trueToMean(t1, ecc), time)


def timeToPrevTrueAnomaly(meanMotion: float, ecc: float, t0: float, epoch0: JulianDate, t1: float,
                          time: JulianDate) -> JulianDate:
    """Computes the previous time the satellite passes through the true anomaly before the time given.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        ecc -- eccentricity of the orbit
        t0 -- a known true anomaly at epoch0 in radians
        epoch0 -- the epoch the satellite has a true anomaly of t0
        t1 -- true anomaly in radians
        time -- relative time to find the previous anomaly
    Returns the previous time the satellites position is t1.
    """

    return timeToPrevMeanAnomaly(meanMotion, trueToMean(t0, ecc), epoch0, trueToMean(t1, ecc), time)


def meanAnomalyAt(meanMotion: float, m0: float, epoch0: JulianDate, epoch: JulianDate) -> float:
    """Computes the mean anomaly of a satellite at a given epoch.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        m0 -- a known mean anomaly at epoch0
        epoch0 -- epoch of m0
        epoch -- epoch to find the mean anomaly
    Returns the mean anomaly in radians.
    """

    dM = meanMotion * epoch.difference(epoch0) * 86400.0
    return (m0 + dM) % TWOPI


def trueAnomalyAt(meanMotion: float, ecc: float, t0: float, epoch0: JulianDate, epoch: JulianDate) -> float:
    """Computes the true anomaly of a satellite at a given epoch.
    Parameters
        meanMotion -- mean motion of the satellite in radians / second
        ecc -- eccentricity of the orbit
        t0 -- a known true anomaly at epoch0
        epoch0 -- epoch of t0
        epoch -- epoch to find the true anomaly
    Returns the true anomaly in radians.
    """

    return meanToTrue(meanAnomalyAt(meanMotion, trueToMean(t0, ecc), epoch0, epoch), ecc)
