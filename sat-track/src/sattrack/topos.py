from math import degrees, asin

from pyevspace import EVector, cross, dot, norm

from sattrack.structures.coordinates import GeoPosition, geoPositionVector, zenithVector
from sattrack.rotation.order import Order
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateToThenOffset
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.satellite import Satellite
from sattrack.util.conversions import atan2


def toTopocentric(vec: EVector, jd: JulianDate, geo: GeoPosition) -> EVector:
    geoVector = geoPositionVector(geo, jd)
    mat = getEulerMatrix(
        Order.ZYX,
        EulerAngles(
            geo.getLongitude() + earthOffsetAngle(jd),
            90 - geo.getLatitude(),
            0.0
        )
    )
    return rotateToThenOffset(mat, geoVector, vec)


# todo: think of a better name for this
def getPVector(geo: GeoPosition, state: tuple[EVector], jd: JulianDate) -> EVector:
    zeta = norm(zenithVector(geo, jd))
    gamma = geoPositionVector(geo, jd)
    lamb = norm(cross(state[0], state[1]))
    x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
    y = -lamb[0] * x / lamb[1]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    return r + v * t


#   todo: moved these, not sure if i want to keep them like this
def getToposPosition(satellite: Satellite, jd: JulianDate, geo: GeoPosition) -> EVector:
    state = satellite.getState(jd)
    return toTopocentric(state[0], jd, geo)


def getAltitude(satellite: Satellite, jd: JulianDate, geo: GeoPosition) -> float:
    state = satellite.getState(jd)
    sez = toTopocentric(state[0], jd, geo)
    return degrees(asin(sez[2] / sez.mag()))


def getAzimuth(satellite: Satellite, jd: JulianDate, geo: GeoPosition) -> float:
    state = satellite.getState(jd)
    sez = toTopocentric(state[0], jd, geo)
    return degrees(atan2(sez[1], -sez[0]))
