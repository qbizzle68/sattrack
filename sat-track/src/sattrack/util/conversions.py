from math import atan2 as at2, pi, sqrt

from sattrack.util.constants import EARTH_MU


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
