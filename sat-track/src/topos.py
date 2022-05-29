from pyevspace import EVector

from coordinates import GeoPosition, geoPositionVector
from order import Order
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
