from math import atan2 as at2, sqrt, pi, sin, cos, radians, degrees
from util.constants import EARTH_MU

''' MOVE THIS TO A MORE APPROPRIATE MODULE'''
def atan2(y: float, x: float) -> float:
    """Computes atan2(y, x) but with the output shifted to (0, 2π).
    Parameters:
    x:  The x coordinate.
    y:  The y coordinate."""
    theta = at2(y, x)
    return (theta + (2 * pi)) if theta < 0 else theta


''' MOVE THIS TO A MORE APPROPRIATE MODULE'''
def meanMotionToSma(meanMotion: float) -> float:
    """Converts mean motion from a two-line element to semi-major axis using
    Kepler's 2nd law.
    Parameters:
    meanMotion: Mean motion measured in revolutions per day.
    Returns the semi-major axis in meters."""
    mMotionRad = meanMotion * 2 * pi / 86400.0
    return (EARTH_MU ** (1.0 / 3.0)) / (mMotionRad ** (2.0 / 3.0))


''' MOVE THIS TO A MORE APPROPRIATE MODULE'''
def smaToMeanMotion(sma: float) -> float:
    """Converts semi-major axis to mean motion.
    Parameters:
    sma: Semi-major axis measured in meters.
    Returns the mean motion in radians per second."""
    return sqrt(EARTH_MU / (sma ** 3))


def meanToTrue(meanAnomaly: float, ecc: float) -> float:
    eTemp = __m2ENewtonRaphson(meanAnomaly, meanAnomaly, ecc) % 360.0
    eAnom = eTemp if eTemp >= 0 else eTemp + 360.0
    y = sqrt(1 - (ecc * ecc)) * sin(radians(eAnom))
    return degrees(atan2(y, cos(radians(eAnom)) - ecc))


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    e = __m2ENewtonRaphson(meanAnomaly, meanAnomaly, eccentricity) % 360.0
    return e if e >= 0 else e + 360.0


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(eccAnom))
    return degrees(atan2(y, cos(radians(eccAnom)) - ecc))


def trueToMean(trueAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    eAnom = degrees(atan2(y, cos(radians(trueAnom)) + ecc))
    m = degrees(radians(eAnom) - ecc * sin(radians(eAnom))) % 360.0
    return m if m >= 0 else m + 360.0


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    return degrees(atan2(y, cos(radians(trueAnom)) + ecc))


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    m = (degrees(radians(eccAnom) - ecc * sin(radians(eccAnom)))) % 360.0
    return m if m >= 0 else m + 360.0


def __m2ENewtonRaphson(M: float, Ej: float, ecc: float) -> float:
    M = radians(M)
    Ej = radians(Ej)
    while True:
        num = Ej - (ecc * sin(Ej)) - M
        den = 1 - ecc * cos(Ej)
        Ej1 = Ej - (num / den)
        if abs(Ej1 - Ej) <= 1e-7:
            return degrees(Ej1)
        else:
            Ej = Ej1

