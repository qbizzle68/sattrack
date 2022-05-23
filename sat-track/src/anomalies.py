from math import atan2 as at2, sqrt, pi, sin, cos, radians, degrees
from constants import EARTH_MU


def atan2(y: float, x: float) -> float:
    """Computes atan2(y, x) but with the output shifted to (0, 2π).
    Parameters:
    x:  The x coordinate.
    y:  The y coordinate."""
    theta = at2(y, x)
    return (theta + (2 * pi)) if theta < 0 else theta


def meanMotionToSma(meanMotion: float) -> float:
    """Converts mean motion from a two-line element to semi-major axis using
    Kepler's 2nd law.
    Parameters:
    meanMotion: Mean motion measured in revolutions per day.
    Returns the semi-major axis in meters."""
    mMotionRad = meanMotion * 2 * pi / 86400.0
    return (EARTH_MU ** (1.0 / 3.0)) / (mMotionRad ** (2.0 / 3.0))


def smaToMeanMotion(sma: float) -> float:
    """Converts semi-major axis to mean motion.
    Parameters:
    sma: Semi-major axis measured in meters.
    Returns the mean motion in radians per second."""
    return sqrt(EARTH_MU / (sma ** 3))


def meanToTrue(meanAnomaly: float, eccentricity: float) -> float:
    return eccentricToTrue(
        meanToEccentric(meanAnomaly, eccentricity),
        eccentricity
    )


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    e = m2ENewtonRaphson(meanAnomaly, meanAnomaly, eccentricity) % 360.0
    return e if e >= 0 else e + 360.0


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(eccAnom))
    return degrees(atan2(y, cos(radians(eccAnom)) - ecc))


def trueToMean(trueAnom: float, ecc: float) -> float:
    return eccentricToMean(
        trueToEccentric(trueAnom, ecc), ecc
    )


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    return degrees(atan2(y, cos(radians(trueAnom)) + ecc))


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    m = (degrees(radians(eccAnom) - ecc * sin(radians(eccAnom)))) % 360.0
    return m if m >= 0 else m + 360.0


def m2ENewtonRaphson(M: float, Ej: float, ecc: float) -> float:
    Ej1 = Ej - ((Ej - ecc * sin(radians(Ej)) - M) / 1 - ecc * cos(radians(Ej)))
    return Ej1 if abs(Ej1 - Ej) <= 1e-7 else m2ENewtonRaphson(M, Ej1, ecc)
