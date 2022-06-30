from math import sqrt, degrees, atan2 as matan2

from pyevspace import EVector, vang, norm, cross

from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.coordinates import GeoPosition, geocentricToGeodetic
from sattrack.structures.elements import computeEccentricVector

from sattrack.util.conversions import atan2 as catan2


def computeTrueAnomaly(position: EVector, velocity: EVector) -> float:
    eccVec = computeEccentricVector(position, velocity)
    # todo: fix this in pyevspace module
    if norm(position) == eccVec:
        return 0
    elif norm(-position) == eccVec:
        return 180
    ang = vang(position, eccVec)
    return 360 - ang if norm(cross(position, eccVec)) == norm(cross(position, velocity)) else ang


def getSubPoint(position: EVector, jd: JulianDate) -> GeoPosition:
    xyMag = sqrt(position[0]*position[0] + position[1]*position[1])
    dec = degrees(matan2(position[2], xyMag))
    lng = (degrees(catan2(position[1], position[0])) - earthOffsetAngle(jd)) % 360.0
    if lng > 180.0:
        lng = lng - 360.0
    return GeoPosition(geocentricToGeodetic(dec), lng)
