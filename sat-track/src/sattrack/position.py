from math import sqrt, degrees, atan2 as matan2, radians, pi

from pyevspace import EVector, vang, norm, cross

from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.coordinates import GeoPosition, geocentricToGeodetic
from sattrack.structures.elements import computeEccentricVector
from sattrack.structures.satellite import Satellite
from sattrack.util.anomalies import trueToMean

from sattrack.util.conversions import atan2 as catan2

# todo: put this in anomalies
def computeTrueAnomaly(position: EVector, velocity: EVector) -> float:
    eccVec = computeEccentricVector(position, velocity)
    # todo: fix this in pyevspace module
    if norm(position) == eccVec:
        return 0
    elif norm(-position) == eccVec:
        return pi
    ang = radians(vang(position, eccVec))
    return (2*pi) - ang if norm(cross(position, eccVec)) == norm(cross(position, velocity)) else ang

# todo: put this in anomalies
def nearestTrueAnomaly(sat: Satellite, time: JulianDate, trueAnom: float) -> JulianDate:
    state = sat.getState(time)
    meanAnom0 = trueToMean(computeTrueAnomaly(state[0], state[1]), sat.tle().getEcc())
    meanAnom1 = trueToMean(trueAnom, sat.tle().getEcc())
    if meanAnom1 < meanAnom0:
        if meanAnom0 - meanAnom1 < pi:
            dma = meanAnom1 - meanAnom0
        else:
            dma = meanAnom1 + (2*pi) - meanAnom0
    else:
        if meanAnom1 - meanAnom0 > pi:
            dma = meanAnom1 - (2*pi) - meanAnom0
        else:
            dma = meanAnom1 - meanAnom0
    n = sat.tle().getMeanMotion() * 2 * pi
    return time.future(dma / n)

# todo: what do we do with this? do we need it?
def getSubPoint(position: EVector, jd: JulianDate) -> GeoPosition:
    xyMag = sqrt(position[0]*position[0] + position[1]*position[1])
    dec = degrees(matan2(position[2], xyMag))
    lng = (degrees(catan2(position[1], position[0]) - earthOffsetAngle(jd))) % 360.0
    if lng > 180.0:
        lng = lng - 360.0
    return GeoPosition(geocentricToGeodetic(dec), lng)
