from math import radians, cos, sin, atan, sqrt

from pyevspace import EVector, vang, norm, cross

from sattrack.structures.elements import computeEccentricVector, OrbitalElements


# todo: why does this return NaN for 0 and 180
from sattrack.util.anomalies import meanToTrue
from sattrack.util.constants import EARTH_MU


def computeTrueAnomaly(position: EVector, velocity: EVector) -> float:
    eccVec = computeEccentricVector(position, velocity)
    ang = vang(position, eccVec)
    return 360 - ang if norm(cross(position, eccVec)) == norm(cross(position, velocity)) else ang


def computeRadius(elements: OrbitalElements) -> float:
    tAnom = meanToTrue(elements.getMeanAnomaly(), elements.getEcc())
    ecc = elements.getEcc()
    return elements.getSma() * (1 - ecc * ecc) / (1 + ecc * cos(radians(tAnom)))


def computeFlightAngle(elements: OrbitalElements) -> float:
    tAnomRad = radians(meanToTrue(elements.getMeanAnomaly(), elements.getEcc()))
    return atan((elements.getEcc() * sin(tAnomRad)) / (1 + elements.getEcc() * cos(tAnomRad)))


def computeVelocity(elements: OrbitalElements) -> float:
    radius = computeRadius(elements)
    return sqrt(EARTH_MU * ((2 / radius) - (1 / elements.getSma())))
