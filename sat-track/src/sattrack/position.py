from math import radians, cos, sin, atan, sqrt, degrees, atan2 as matan2

from pyevspace import EVector, vang, norm, cross

from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.coordinates import GeoPosition, geocentricToGeodetic
from sattrack.structures.elements import computeEccentricVector, OrbitalElements


from sattrack.util.anomalies import meanToTrue
from sattrack.util.constants import EARTH_MU
from sattrack.util.conversions import smaToMeanMotion, atan2 as catan2


def timeToMeanAnomaly(elements: OrbitalElements, meanAnom: float) -> float:
    n = smaToMeanMotion(elements.getSma())
    dM = (meanAnom - elements.getMeanAnomaly()) % 360.0
    return radians(dM) / n


def meanAnomalyAt(elements: OrbitalElements, jd: JulianDate) -> float:
    if elements.getEpoch() == 0:
        raise ValueError('Epoch was not set for this instance.')
    n = smaToMeanMotion(elements.getSma())
    dt = jd.difference(elements.getEpoch())
    mAnomRad = n * dt + radians(elements.getMeanAnomaly())
    return degrees(mAnomRad) % 360.0


# todo: why does this return NaN for 0 and 180
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


def getSubPoint(position: EVector, velocity: EVector, jd: JulianDate) -> GeoPosition:
    xyMag = sqrt(position[0]*position[0] + position[1]*position[1])
    dec = degrees(matan2(position[2], xyMag))
    lng = (degrees(catan2(position[1], position[0])) - earthOffsetAngle(jd)) % 360.0
    if lng > 180.0:
        lng = lng - 360.0
    return GeoPosition(geocentricToGeodetic(dec), lng)
