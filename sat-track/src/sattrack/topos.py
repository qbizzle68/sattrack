from pyevspace import EVector, cross, dot, norm

from structures.coordinates import GeoPosition, geoPositionVector, zenithVector
from rotation.order import Order
from rotation import getEulerMatrix, EulerAngles, rotateToThenOffset
from spacetime import JulianDate, earthOffsetAngle


def ijkToSEZ(vec: EVector, jd: JulianDate, geo: GeoPosition) -> EVector:
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


def getPVector(geo: GeoPosition, state: tuple[EVector], jd: JulianDate) -> EVector:
    zeta = zenithVector(geo, jd)
    gamma = geoPositionVector(geo, jd)
    #state = sat.getState(jd)
    lamb = norm(cross(state[0], state[1]))
    x = dot(zeta, gamma) * (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
    y = -lamb[0] * x / lamb[1]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)
    t = dot(r, (gamma - r)) / dot(r, r)
    return r + v * t
