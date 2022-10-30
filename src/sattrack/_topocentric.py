from math import radians, pi

from sattrack.rotation.order import ZYX
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateToThenOffset
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures._coordinates import _compute_position_vector


def _to_topocentric(position, geo, time):
    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    geoVector = _compute_position_vector(latitude, longitude, geo.elevation, time)
    matrix = getEulerMatrix(
        ZYX,
        EulerAngles(
            longitude + earthOffsetAngle(time),
            pi / 2 - latitude,
            0.0
        )
    )
    return rotateToThenOffset(matrix, geoVector, position)
