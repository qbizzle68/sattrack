from math import radians, pi

from pyevspace import getMatrixEuler, ZYX, Angles, rotateOffsetTo, rotateEulerTo
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack._coordinates import _computePositionVector


def __getTopocentricAngles(geo, time):
    """Computes the Euler angles required for a rotation to a topocentric reference frame."""

    lng = radians(geo.longitude) + earthOffsetAngle(time)
    lat = pi / 2 - radians(geo.latitude)
    return Angles(lng, lat, 0.0)


def _toTopocentricOffset(position, geo, time):
    """Converts a position vector to a topocentric (SEZ) reference frame by offsetting by the
    geo-position's vector."""

    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)

    geoVector = _computePositionVector(latitude, longitude, geo.elevation, time)
    angs = __getTopocentricAngles(geo, time)
    matrix = getMatrixEuler(ZYX, angs)
    return rotateOffsetTo(matrix, geoVector, position)


def _toTopocentric(position, geo, time):
    """Converts a vector to a topocentric reference frame that doesn't require offsetting, such as a
    velocity vector."""

    angs = __getTopocentricAngles(geo, time)
    return rotateEulerTo(ZYX, angs, position)
