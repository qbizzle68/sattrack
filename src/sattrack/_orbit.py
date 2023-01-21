from math import sqrt, sin, cos, pi

from pyevspace import Vector, dot, vang, norm, cross
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import atan3
from sattrack.sgp4 import computeEccentricVector


def _vectorAlmostEqual(lhs: Vector, rhs: Vector, epsilon=1e-5):
    """Compares Vectors that are very close to equal but have enough rounding error
    to evaluate to False. This should be addressed in pyevspace but is a patch for now."""
    if lhs == rhs:
        return True

    if abs(lhs.mag2() - rhs.mag2()) < epsilon:
        relativeProj = dot(lhs, rhs) / (lhs.mag() * rhs.mag())
        return relativeProj > (1 - epsilon)

    return False


def _trueAnomalyFromState(position: Vector, velocity: Vector, mu: float) -> float:
    """Computes the true anomaly of a position based on its state vectors."""
    eccentricVector = computeEccentricVector(position, velocity, mu)
    trueAnomaly = vang(position, eccentricVector)
    if _vectorAlmostEqual(norm(cross(position, eccentricVector)), norm(cross(position, velocity))):
        trueAnomaly = TWOPI - trueAnomaly
    return trueAnomaly


def _trueToMeanAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly in radians to an eccentric anomaly in radians."""
    eccentricAnomaly = _trueToEccentricAnomaly(trueAnomaly, eccentricity)
    return _eccentricToMeanAnomaly(eccentricAnomaly, eccentricity)


def _trueToEccentricAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly to an eccentric anomaly in radians."""
    y = sqrt(1 - (eccentricity * eccentricity)) * sin(trueAnomaly)
    return atan3(y, cos(trueAnomaly) + eccentricity)


def _eccentricToMeanAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly to a mean anomaly in radians."""
    return eccentricAnomaly - eccentricity * sin(eccentricAnomaly)


def _smaToMeanMotion(semiMajorAxis: float, mu: float) -> float:
    """Converts a semi-major axis in kilometers to mean motion in radians / second."""
    return sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis))


def _nearestTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0, t1: float):
    """Finds the time of the nearest true anomaly (radians), either forward or backward in time. Mean motion
    must be in revolutions / day."""
    m0 = _trueToMeanAnomaly(t0, eccentricity)
    m1 = _trueToMeanAnomaly(t1, eccentricity)
    return _nearestMeanAnomaly(meanMotion, m0, epoch0, m1)


def _nearestMeanAnomaly(meanMotion: float, m0: float, epoch0, m1: float):
    """Finds the time of the nearest mean anomaly (radians), either forward or backward in time. Mean motion
    must be in revolutions / day."""
    if m1 < m0:
        if m0 - m1 < pi:
            dma = m1 - m0
        else:
            dma = m1 + TWOPI - m0
    else:
        if m1 - m0 > pi:
            dma = m1 - TWOPI - m0
        else:
            dma = m1 - m0
    return epoch0.future(dma / (meanMotion * TWOPI))
