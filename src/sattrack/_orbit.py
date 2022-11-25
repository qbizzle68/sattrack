from math import sqrt, sin, cos, acos, pi

from pyevspace import EVector, dot, vang, norm, cross
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import atan3
from sattrack._sgp4 import _compute_eccentric_vector, _elements_from_state


# def _compute_eccentric_vector(position: EVector, velocity: EVector, MU: float) -> EVector:
#     lhs = position * ((velocity.mag2() / MU) - (1 / position.mag()))
#     rhs = velocity * (dot(position, velocity) / MU)
#     return lhs - rhs


def _true_anomaly_from_state(position: EVector, velocity: EVector, mu: float) -> float:
    # todo: we can make a mean anomaly from state if we can guarantee eccentricVector.mag() is accurate
    eccentricVector = _compute_eccentric_vector(position, velocity, mu)
    trueAnomaly = vang(position, eccentricVector)
    if norm(cross(position, eccentricVector)) == norm(cross(position, velocity)):
        trueAnomaly = TWOPI - trueAnomaly
    return trueAnomaly


def _true_to_mean_anomaly(trueAnomaly: float, eccentricity: float) -> float:
    eccentricAnomaly = _true_to_eccentric_anomaly(trueAnomaly, eccentricity)
    return _eccentric_to_mean_anomaly(eccentricAnomaly, eccentricity)


def _true_to_eccentric_anomaly(trueAnomaly: float, eccentricity: float) -> float:
    y = sqrt(1 - (eccentricity * eccentricity)) * sin(trueAnomaly)
    return atan3(y, cos(trueAnomaly) + eccentricity)


def _eccentric_to_mean_anomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    return eccentricAnomaly - eccentricity * sin(eccentricAnomaly)


# def _elements_from_state(position: EVector, velocity: EVector, MU: float) -> (float, float, float, float, float, float):
#     angularMomentum = cross(position, velocity)
#     lineOfNodes = norm(cross(EVector.e3, angularMomentum))
#     eccentricityVector = _compute_eccentric_vector(position, velocity, MU)
#
#     ecc = eccentricityVector.mag()
#     inc = acos(angularMomentum[2] / angularMomentum.mag())
#     raan = acos(lineOfNodes[0])
#     if lineOfNodes[1] < 0:
#         raan = TWOPI - raan
#     aop = acos(dot(lineOfNodes, eccentricityVector) / eccentricityVector.mag())
#     if eccentricityVector[2] < 0:
#         aop = TWOPI - aop
#     tAnom = acos(dot(eccentricityVector, position) / (eccentricityVector.mag() * position.mag()))
#     if dot(position, velocity) < 0:
#         tAnom = 2 * pi - tAnom
#     mAnom = _true_to_mean_anomaly(tAnom, ecc)
#     sma = (angularMomentum.mag() ** 2) / ((1 - (ecc * ecc)) * MU)
#
#     return raan, inc, aop, ecc, sma, mAnom


def _sma_to_mean_motion(semiMajorAxis: float, mu: float) -> float:
    # semiMajorAxis - km, return radian / s
    return sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis))


def _nearest_true_anomaly(meanMotion: float, eccentricity: float, t0: float, epoch0, t1: float):
    m0 = _true_to_mean_anomaly(t0, eccentricity)
    m1 = _true_to_mean_anomaly(t1, eccentricity)
    return _nearest_mean_anomaly(meanMotion, m0, epoch0, m1)


def _nearest_mean_anomaly(meanMotion: float, m0: float, epoch0, m1: float):
    # meanMotion - rev / day; m0, m1 - radians
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
