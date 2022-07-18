from math import sqrt, sin, cos, floor

from sattrack.spacetime.juliandate import JulianDate
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import atan2


def meanToTrue(meanAnomaly: float, ecc: float) -> float:
    """
    Converts a mean anomaly to the same position's true anomaly.

    Args:
        meanAnomaly: The mean anomaly in radians.
        ecc: Eccentricity of the orbit.

    Returns:
        The true anomaly in radians.
    """

    eAnom = __m2ENewtonRaphson(meanAnomaly, ecc) % TWOPI
    y = sqrt(1 - (ecc * ecc)) * sin(eAnom)
    return atan2(y, cos(eAnom) - ecc)


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    """
    Converts a mean anomaly to the same position's eccentric anomaly.

    Args:
        meanAnomaly: The mean anomaly in radians.
        eccentricity: Eccentricity of the orbit.

    Returns:
        The eccentric anomaly measured in radians.
    """

    return __m2ENewtonRaphson(meanAnomaly, eccentricity) % TWOPI


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    """
    Converts an eccentric anomaly to the same position's true anomaly.

    Args:
        eccAnom: The eccentric anomaly in radians.
        ecc: Eccentricity of the orbit.

    Returns:
        The true anomaly in radians.
    """

    y = sqrt(1 - (ecc * ecc)) * sin(eccAnom)
    return atan2(y, cos(eccAnom) - ecc)


def trueToMean(trueAnom: float, ecc: float) -> float:
    """
    Converts a true anomaly to the same position's mean anomaly.

    Args:
        trueAnom: The true anomaly in radians.
        ecc: Eccentricity of the orbit.

    Returns:
        The mean anomaly in radians.
    """

    y = sqrt(1 - (ecc * ecc)) * sin(trueAnom)
    eAnom = atan2(y, cos(trueAnom) + ecc)
    return eAnom - ecc * sin(eAnom) % TWOPI


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    """
    Converts a true anomaly to the same position's eccentric anomaly.

    Args:
        trueAnom: The true anomaly in radians.
        ecc: Eccentricity of the orbit.

    Returns:
        The eccentric anomaly in radians.
    """

    y = sqrt(1 - (ecc * ecc)) * sin(trueAnom)
    return atan2(y, cos(trueAnom) + ecc)


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    """
    Converts an eccentric anomaly to the same position's mean anomaly.

    Args:
        eccAnom: The eccentric anomaly in radians.
        ecc: Eccentric of the orbit.

    Returns:
        The mean anomaly in radians.
    """

    return (eccAnom - ecc * sin(eccAnom)) % TWOPI


def __m2ENewtonRaphson(M: float, ecc: float) -> float:
    """
    Computes the eccentric anomaly from a given mean anomaly via Kepler's equation. The computation is done using the
    Newton-Raphson method to a hard-coded epsilon value of 1e-7.

    Args:
        M: The mean anomaly measured in radians.
        ecc: Eccentricity of the orbit.

    Returns:
        An approximation of the eccentric anomaly.
    """

    # initial guess = M
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
    """
    Computes the next time a satellite passes through the mean anomaly after the given time.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        m0: A known mean anomaly at epoch0 in radians.
        epoch0: The epoch the satellite has a mean anomaly of m0.
        m1: Mean anomaly in radians.
        time: Relative time to find the next anomaly.

    Returns:
        The next time the satellite's position is m1.
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
    """
    Computes the previous time the satellite passed through the mean anomaly before the given time.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        m0: A known mean anomaly at epoch0 in radians.
        epoch0: The epoch the satellite has a mean anomaly of m0.
        m1: Mean anomaly in radians.
        time: Relative time to find the previous anomaly.

    Returns:
        The previous time the satellite's position was m1.
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
    """
    Computes the next time the satellite passes through the true anomaly after the time given.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        ecc: Eccentricity of the orbit.
        t0: A known true anomaly at epoch0 in radians.
        epoch0: The epoch the satellite has a true anomaly of t0.
        t1: True anomaly in radians.
        time: Relative time to find the next anomaly.

    Returns:
        The next time the satellite's position is t1.
    """

    return timeToNextMeanAnomaly(meanMotion, trueToMean(t0, ecc), epoch0, trueToMean(t1, ecc), time)


def timeToPrevTrueAnomaly(meanMotion: float, ecc: float, t0: float, epoch0: JulianDate, t1: float,
                          time: JulianDate) -> JulianDate:
    """
    Computes the previous time the satellite passes through the true anomaly before the time given.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        ecc: Eccentricity of the orbit.
        t0: A known true anomaly at epoch0 in radians.
        epoch0: The epoch the satellite has a true anomaly of t0.
        t1: True anomaly in radians.
        time: Relative time to find the previous anomaly.

    Returns:
        The previous time the satellite position is t1.
    """

    return timeToPrevMeanAnomaly(meanMotion, trueToMean(t0, ecc), epoch0, trueToMean(t1, ecc), time)


def meanAnomalyAt(meanMotion: float, m0: float, epoch0: JulianDate, epoch: JulianDate) -> float:
    """
    Computes the mean anomaly of a satellite at a given epoch.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        m0: A known mean anomaly at epoch0.
        epoch0: Epoch of m0.
        epoch: Epoch to find the mean anomaly.

    Returns:
        The mean anomaly in radians.
    """

    dM = meanMotion * epoch.difference(epoch0) * 86400.0
    return (m0 + dM) % TWOPI


def trueAnomalyAt(meanMotion: float, ecc: float, t0: float, epoch0: JulianDate, epoch: JulianDate) -> float:
    """
    Computes the true anomaly of a satellite at a given epoch.

    Args:
        meanMotion: Mean motion of the satellite in radians per second.
        ecc: Eccentricity of the orbit.
        t0: A known true anomaly at epoch0.
        epoch0: Epoch of t0.
        epoch: Epoch to find the true anomaly.

    Returns:
        The true anomaly in radians.
    """

    return meanToTrue(meanAnomalyAt(meanMotion, trueToMean(t0, ecc), epoch0, epoch), ecc)
