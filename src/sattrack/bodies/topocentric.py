from math import asin, degrees, pi

from pyevspace import Vector, ZYX, getMatrixEuler, Angles, rotateEulerTo, rotateOffsetTo, rotateOffsetFrom, cross

from sattrack.bodies.position import computeEarthOffsetAngle
from sattrack.util.constants import TWOPI, EARTH_SIDEREAL_PERIOD
from sattrack.util.helpers import atan3

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.coordinates import GeoPosition
    from sattrack.orbit.satellite import Orbitable
    from sattrack.core.juliandate import JulianDate


def _getTopocentricAngles(geo: 'GeoPosition', time: 'JulianDate') -> Angles:
    """Computes the Euler angles required for a rotation to a topocentric reference frame."""

    # lng = geo.longitudeRadians + Earth.offsetAngle(time)
    lng = geo.longitudeRadians + computeEarthOffsetAngle(time)
    lat = pi / 2 - geo.latitudeRadians

    return Angles(lng, lat, 0.0)


def toTopocentric(vector: Vector, geo: 'GeoPosition', time: 'JulianDate') -> Vector:
    """Converts a vector to a topocentric reference frame that doesn't require offsetting, such as a
    velocity vector."""

    angs = _getTopocentricAngles(geo, time)
    return rotateEulerTo(ZYX, angs, vector)


# todo: add an offset keyword parameter to toTopocentric to replace this method
def toTopocentricOffset(position: Vector, geo: 'GeoPosition', time: 'JulianDate') -> Vector:
    """Converts a position vector to a topocentric (SEZ) reference frame by offsetting by the
    geo-position's vector."""

    geoVector = geo.getPositionVector(time)
    angs = _getTopocentricAngles(geo, time)
    matrix = getMatrixEuler(ZYX, angs)

    return rotateOffsetTo(matrix, geoVector, position)


def toTopocentricState(position: Vector, velocity: Vector, geo: 'GeoPosition', time: 'JulianDate') -> (Vector, Vector):
    """Converts state vectors to the topocentric reference frame defined by geo and time."""

    w = TWOPI / EARTH_SIDEREAL_PERIOD
    angularMomentum = Vector(0, 0, w)
    coriolisEffect = cross(angularMomentum, position)
    topocentricPosition = toTopocentricOffset(position, geo, time)
    topocentricVelocity = toTopocentric(velocity - coriolisEffect, geo, time)

    return topocentricPosition, topocentricVelocity


# todo: add an offset argument to this and add appropriate logic for it
def fromTopocentric(vector: Vector, geo: 'GeoPosition', time: 'JulianDate') -> Vector:
    """Converts a vector from a topocentric reference frame defined by a geo-position at a specified time."""

    geoVector = geo.getPositionVector(time)
    # angs = Angles(geo.longitudeRadians + Earth.offsetAngle(time), pi / 2 - geo.latitudeRadians, 0.0)
    angs = Angles(geo.longitudeRadians + computeEarthOffsetAngle(time), pi / 2 - geo.latitudeRadians, 0.0)
    matrix = getMatrixEuler(ZYX, angs)

    return rotateOffsetFrom(matrix, geoVector, vector)


def getAltitude(satellite: 'Orbitable', geo: 'GeoPosition', time: 'JulianDate') -> float:
    """Returns the altitude of an orbitable from a geo-position at a specified time."""

    position = satellite.getState(time)[0]
    topocentricPosition = toTopocentricOffset(position, geo, time)
    arg = topocentricPosition[2] / topocentricPosition.mag()
    return degrees(asin(arg))


def getAzimuth(satellite: 'Orbitable', geo: 'GeoPosition', time: 'JulianDate') -> float:
    """Returns the azimuth of an orbitable from a geo-position at a specified time."""

    position = satellite.getState(time)[0]
    topocentricPosition = toTopocentricOffset(position, geo, time)
    return degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))
