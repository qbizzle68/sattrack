from math import asin, degrees, pi

from pyevspace import Vector, ZYX, getMatrixEuler, Angles, rotateEulerTo, rotateOffsetTo, rotateOffsetFrom, cross

from sattrack.core.sidereal import earthOffsetAngle
from sattrack.util.constants import TWOPI, EARTH_SIDEREAL_PERIOD
from sattrack.util.helpers import atan3

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.coordinates import GeoPosition
    from sattrack.orbit.satellite import Orbitable
    from sattrack.core.juliandate import JulianDate


def _getTopocentricAngles(geo: 'GeoPosition', time: 'JulianDate') -> Angles:
    """Computes the Euler angles required for a rotation to a topocentric reference frame."""

    lng = geo.longitudeRadians + earthOffsetAngle(time)
    lat = pi / 2 - geo.latitudeRadians

    return Angles(lng, lat, 0.0)


def toTopocentric(vector: Vector, geo: 'GeoPosition', time: 'JulianDate') -> Vector:
    """Converts a vector to a topocentric reference frame that doesn't require offsetting, such as a
    velocity vector."""

    angs = _getTopocentricAngles(geo, time)
    return rotateEulerTo(ZYX, angs, vector)


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


def fromTopocentric(vector: Vector, geo: 'GeoPosition', time: 'JulianDate') -> Vector:
    """Converts a vector from a topocentric reference frame defined by a geo-position at a specified time."""

    geoVector = geo.getPositionVector(time)
    angs = Angles(geo.longitudeRadians + earthOffsetAngle(time), pi / 2 - geo.latitudeRadians, 0.0)
    # angs = Angles(radians(geo.longitude) + earthOffsetAngle(time), radians(90 - geo.latitude), 0.0)
    matrix = getMatrixEuler(ZYX, angs)

    return rotateOffsetFrom(matrix, geoVector, vector)

# todo: make a fromTopocentricOffset() function


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


class AltAz:
    __slots__ = '_altitude', '_azimuth', '_direction'

    def __init__(self, altitude: float, azimuth: float):
        # Args in degrees, will be modded to valid ranges.

        if not -90 <= altitude <= 90:
            raise ValueError(f'altitude must be in (-90, 90), was {altitude}')

        self._altitude = altitude
        self._azimuth = azimuth % 360.0
        self._direction = self.azimuthAngleString(self._azimuth)

    @staticmethod
    def azimuthAngleString(azimuth: float):
        """Converts an azimuth angle in degrees to a compass direction string."""
        if azimuth > 348.25 or azimuth <= 11.25:
            return 'N'
        elif azimuth <= 33.75:
            return 'NNE'
        elif azimuth <= 56.25:
            return 'NE'
        elif azimuth <= 78.75:
            return 'ENE'
        elif azimuth <= 101.25:
            return 'E'
        elif azimuth <= 123.75:
            return 'ESE'
        elif azimuth <= 146.25:
            return 'SE'
        elif azimuth <= 168.75:
            return 'SSE'
        elif azimuth <= 191.25:
            return 'S'
        elif azimuth <= 213.75:
            return 'SSW'
        elif azimuth <= 236.25:
            return 'SW'
        elif azimuth <= 258.75:
            return 'WSW'
        elif azimuth <= 281.25:
            return 'W'
        elif azimuth <= 303.75:
            return 'WNW'
        elif azimuth <= 326.25:
            return 'NW'
        else:
            return 'NNW'

    def toDict(self):
        return {'az': self._azimuth, 'alt': self._altitude, 'direction': self._direction}

    def __str__(self):
        return str(self.toDict())

    def __repr__(self):
        return f'AltAz({self._altitude}, {self._azimuth})'

    @property
    def altitude(self):
        return self._altitude

    @property
    def azimuth(self):
        return self._azimuth

    @property
    def direction(self):
        return self._direction
