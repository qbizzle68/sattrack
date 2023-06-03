from math import radians, pi

from pyevspace import getMatrixEuler, ZYX, Angles, rotateEulerTo, cross, Vector, rotateOffsetTo

from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack._coordinates import _computePositionVector
from sattrack.util.constants import EARTH_SIDEREAL_PERIOD, TWOPI

__all__ = ['__getTopocentricAngles', '_toTopocentricOffset', '_toTopocentric', '_toTopocentricState']


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


def _toTopocentric(vector, geo, time):
    """Converts a vector to a topocentric reference frame that doesn't require offsetting, such as a
    velocity vector."""

    angs = __getTopocentricAngles(geo, time)
    return rotateEulerTo(ZYX, angs, vector)


def _toTopocentricState(position, velocity, geo, time):
    """Converts state vectors to the topocentric reference frame defined by geo and time."""

    w = TWOPI / EARTH_SIDEREAL_PERIOD
    angularMomentum = Vector(0, 0, w)
    coriolisEffect = cross(angularMomentum, position)
    topocentricPosition = _toTopocentricOffset(position, geo, time)
    topocentricVelocity = _toTopocentric(velocity - coriolisEffect, geo, time)

    return topocentricPosition, topocentricVelocity
