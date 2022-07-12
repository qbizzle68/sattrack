from math import sqrt, sin, cos, radians, degrees, pi, floor

from sattrack.spacetime.juliandate import JulianDate
from sattrack.util.conversions import atan2


def meanToTrue(meanAnomaly: float, ecc: float) -> float:
    """Converts a mean anomaly to its true anomaly.
    Parameters:
    meanAnomaly:    The mean anomaly measured in degrees.
    ecc:            Eccentricity of the orbit.
    returns:        The true anomaly measured in degrees."""

    eTemp = __m2ENewtonRaphson(meanAnomaly, ecc) % 360.0
    eAnom = eTemp if eTemp >= 0 else eTemp + 360.0
    y = sqrt(1 - (ecc * ecc)) * sin(radians(eAnom))
    return degrees(atan2(y, cos(radians(eAnom)) - ecc))


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly to its eccentric anomaly.
    Parameters:
    meanAnomaly:    The mean anomaly measured in degrees.
    ecc:            Eccentricity of the orbit.
    returns:        The eccentric anomaly measured in degrees."""

    e = __m2ENewtonRaphson(meanAnomaly, eccentricity) % 360.0
    return e if e >= 0 else e + 360.0


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    """Converts an eccentric anomaly to its true anomaly.
    Parameters:
    eccAnom:    The eccentric anomaly measured in degrees.
    ecc:        Eccentricity of the orbit.
    returns:    The true anomaly measured in degrees."""

    y = sqrt(1 - (ecc * ecc)) * sin(radians(eccAnom))
    return degrees(atan2(y, cos(radians(eccAnom)) - ecc))


def trueToMean(trueAnom: float, ecc: float) -> float:
    """Converts a true anomaly to its mean anomaly.
    Parameters:
    trueAnom:   The true anomaly measured in degrees.
    ecc:        Eccentricity of the orbit.
    returns:    The mean anomaly measured in degrees."""

    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    eAnom = degrees(atan2(y, cos(radians(trueAnom)) + ecc))
    m = degrees(radians(eAnom) - ecc * sin(radians(eAnom))) % 360.0
    return m if m >= 0 else m + 360.0


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    """Converts a true anomaly to its eccentric anomaly.
    Parameters:
    trueAnom:   The true anomaly measured in degrees.
    ecc:        Eccentricity of the orbit.
    returns:    The eccentric anomaly measured in degrees."""

    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    return degrees(atan2(y, cos(radians(trueAnom)) + ecc))


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    """Converts an eccentric anomaly to its mean anomaly.
    Parameters:
    eccAnom:    The eccentric anomaly measured in degrees.
    ecc:        Eccentricity of the orbit.
    returns:    The mean anomaly measured in degrees."""

    m = (degrees(radians(eccAnom) - ecc * sin(radians(eccAnom)))) % 360.0
    return m if m >= 0 else m + 360.0


def __m2ENewtonRaphson(M: float, ecc: float) -> float:
    """Computes the eccentric anomaly from a given mean anomaly via Kepler's equation.
    The computation is done using the Newton-Raphson method to a hard-coded epsilon value
    of 1e-7. Inputs and outputs are degrees, the conversion is handed in the implementation.
    Parameters:
    M:      The mean anomaly measured in degrees.
    ecc:    Eccentricity of the orbit."""

    M = radians(M)
    Ej = M
    while True:
        num = Ej - (ecc * sin(Ej)) - M
        den = 1 - ecc * cos(Ej)
        Ej1 = Ej - (num / den)
        if abs(Ej1 - Ej) <= 1e-7:
            return degrees(Ej1)
        else:
            Ej = Ej1


def timeToNextMeanAnomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    """meanMotion in rad/s. m0 - epoch0 mean anomaly and time, time to m1 after time"""
    twoPi = 2 * pi
    n = meanMotion * 86400 / twoPi
    pePass = epoch0.future(-radians(m0) / (twoPi * n))
    revs = time.difference(pePass) * n
    m0 = (revs - floor(revs)) * twoPi
    m1 = radians(m1)
    if m1 < m0:
        dm = m1 + twoPi - m0
    else:
        dm = m1 - m0
    return time.future((dm / n) / twoPi)
