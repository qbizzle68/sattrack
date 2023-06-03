import json
import abc as _abc
import math as _math

from sattrack._coordinates import _computeZenithVector, _computePositionVector, _radiusAtLat, \
    _geodeticToGeocentric
from sattrack.util.constants import EARTH_FLATTENING, TWOPI
from pyevspace import Vector
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.conversions import atan3

if __debug__ is True:
    debug__all__ = ['_RAD_TO_HOURS', '_DEG_TO_HOURS', '_geocentricToGeodetic']
else:
    debug__all__ = []

__all__ = ['Coordinates', 'GeoPosition', 'CelestialCoordinates', 'radiusAtLatitude', 'geocentricToGeodetic',
           'geodeticToGeocentric', 'getSubPoint'] + debug__all__

_RAD_TO_HOURS = 12 / _math.pi
_DEG_TO_HOURS = 15


class Coordinates(_abc.ABC):
    """Base class for all coordinate type classes."""

    __slots__ = '_lat', '_lng'

    def __init__(self, latitude: float, longitude: float):
        """Initialize coordinates to latitude and longitudes in degrees."""
        self._lat = self._checkLatitude(latitude)
        self._lng = self._checkLongitude(longitude)

    # def __repr__(self) -> str:
    #     return _json.dumps(self, default=lambda o: dict(o))

    @property
    def latitude(self):
        return self._fmtLatitude()

    @latitude.setter
    def latitude(self, value):
        self._lat = self._checkLatitude(value)

    @property
    def longitude(self):
        return self._fmtLongitude()

    @longitude.setter
    def longitude(self, value):
        self._lng = self._checkLongitude(value)

    @_abc.abstractmethod
    def _fmod(self, val: float) -> float:
        """Modulates a longitude between valid values."""
        pass

    @_abc.abstractmethod
    def _checkLatitude(self, lat):
        pass

    @_abc.abstractmethod
    def _checkLongitude(self, lng):
        """Should call self._fmod() before returning a value."""
        pass

    @_abc.abstractmethod
    def _fmtLatitude(self):
        """Formats the latitude before returning it for display purposes."""
        pass

    @_abc.abstractmethod
    def _fmtLongitude(self):
        """Formats the longitude before returning it for display purposes."""
        pass


class GeoPosition(Coordinates):
    """
    A container class for an Earth geo-position. The class contains a latitude, longitude and an optional elevation. The
    elevation is used for positions relative to a geo-position, as well as rotating to a topocentric reference frame.
    Since the Earth's shape is approximated as an oblate spheroid, there are two types of latitudes, geocentric and
    geodetic. Geocentric latitude is the angle between the equator and a point on Earth's surface through the Earth's
    center. Geodetic latitude is the angle between the equator and a line normal to Earth's surface at a given point.
    The geodetic latitude is also known as geographic latitude, and is the one used in a geo-position coordinate.
    Latitudes are measured zero at the equator, positive to the north, and negative to the south. Longitude is measured
    zero on the longitude passing through Greenwich, England, and is positive to the east, and negative to the west.
    """

    __slots__ = '_elevation',

    def __init__(self, latitude: float, longitude: float, elevation: float = 0):
        """Initializes the GeoPosition with the position angles in degrees and optional elevation in kilometers."""
        super().__init__(latitude, longitude)
        if not isinstance(elevation, (int, float)):
            raise TypeError('elevation parameter must be an int or float type', type(elevation))
        self._elevation = elevation

    def __str__(self):
        """Returns a string representation of the geo-position."""
        return f'latitude: {self._fmtLatitude()}, longitude: {self._fmtLongitude()}, elevation: {self._elevation}'

    def __repr__(self):
        """Returns a string representation of the geo-position."""
        return f'GeoPosition({self._fmtLatitude()}, {self._fmtLongitude()}, {self._elevation})'

    def __reduce__(self):
        """Allows the GeoPosition to be pickled."""
        return self.__class__, (self._fmtLatitude(), self._fmtLongitude(), self._elevation)

    def toJson(self):
        """Returns a string of the GeoPosition in json format."""
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        """Returns a dictionary of the GeoPosition to create json formats of other types containing a GeoPosition."""
        return {"latitude": self._lat, "longitude": self._lng, "elevation": self._elevation}

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError('value must be an int or float type')
        self._elevation = value

    def getRadius(self) -> float:
        """Returns the radius of the position including elevation in kilometers."""
        return _radiusAtLat(self._lat) + self._elevation

    def getVelocityVector(self):
        """Computes the velocity vector in kilometers / second of the GeoPosition in the topocentric reference frame."""
        vel = self.getRadius() * TWOPI / 86164.090531
        return Vector(0, vel, 0)

    def getGeocentricLatitude(self) -> float:
        """Computes the geocentric latitude of the position."""
        return _geodeticToGeocentric(self._lat)

    def getPositionVector(self, jd: JulianDate = None) -> Vector:
        """Gets the spatial position vector of the GeoPosition. If jd is not None the offset angle of the Earth is taken
        into account."""
        if jd is not None and not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type', type(jd))
        return _computePositionVector(self._lat, self._lng, self._elevation, jd)

    def getZenithVector(self, jd: JulianDate = None) -> Vector:
        """Gets the normal vector to the Earth's surface at this GeoPosition. The zenith vector is guaranteed to have
        a magnitude of 1."""
        if jd is not None and not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type', type(jd))
        return _computeZenithVector(self._lat, self._lng, self._elevation, jd)

    def _fmod(self, value: float) -> float:
        """Modulates the latitude value from -180 to 180."""
        if value > 180.0:
            tmp = value - int(value / 360.0) * 360.0
            return tmp - 360.0 if tmp > 180.0 else tmp
        elif value < -180.0:
            tmp = value - int(value / 360.0) * 360.0
            return tmp + 360.0 if tmp < -180.0 else tmp
        else:
            return value

    def _checkLatitude(self, latitude):
        if not isinstance(latitude, (int, float)):
            raise TypeError('lat parameter must be an int or float type', type(latitude))
        # a latitude outside the bounds [-90, 90] doesn't make any sense, so raise error if it occurs
        if not -90 <= latitude <= 90:
            raise ValueError('lat argument must be between -90 and 90 degrees', latitude)
        return _math.radians(latitude)

    def _checkLongitude(self, longitude):
        if not isinstance(longitude, (int, float)):
            raise TypeError('lng parameter must be an int or float type', type(longitude))
        return _math.radians(self._fmod(longitude))

    def _fmtLatitude(self):
        """Re-convert latitude back to degrees for printing."""
        return _math.degrees(self._lat)

    def _fmtLongitude(self):
        """Re-convert longitude back to degrees for printing."""
        return _math.degrees(self._lng)


def _geocentricToGeodetic(geocentricLatitude):
    """Converts a geocentric latitude in radians to a geodetic latitude in radians."""
    if not -_math.pi / 2 <= geocentricLatitude <= _math.pi / 2:
        raise ValueError('lat argument must be between -90 and 90 degrees', geocentricLatitude)
    return _math.atan(_math.tan(geocentricLatitude) / ((1 - EARTH_FLATTENING) ** 2))


def geocentricToGeodetic(geocentricLatitude: float) -> float:
    """Converts a geocentric latitude in degrees to a geodetic latitude in degrees."""
    return _math.degrees(_geocentricToGeodetic(_math.radians(geocentricLatitude)))


def geodeticToGeocentric(geodeticLatitude: float) -> float:
    """Converts a geodetic latitude in degrees to a geocentric latitude in degrees."""
    # parameter and return value are in degrees
    return _math.degrees(_geodeticToGeocentric(_math.radians(geodeticLatitude)))


def radiusAtLatitude(latitude: float) -> float:
    """Computes the radius of the Earth in kilometers from a given geodetic latitude in degrees."""
    return _radiusAtLat(_math.radians(latitude))


class CelestialCoordinates(Coordinates):
    """
    A container clas for a set of celestial coordinates. The coordinate names are declination and right-ascension.
    Declination is the equivalent of latitude, and is the angle between the celestial equator and a point in space.
    Right-ascension is the equivalent of longitude, and is the angle from the first point of Ares, or the x-axis. It is
    measured in time units, from 0 to 24 hours, and is positive to the east.
    """

    def __init__(self, rightAscension: float, declination: float):
        """Initializes the coordinate object with the right-ascension in hours and declination in degrees. Note the
        order of the parameters follows the idiom 'right-ascension and declination', not 'latitude and longitude'."""
        super().__init__(declination, rightAscension)

    def __str__(self):
        """Returns a string representation of the coordinates."""
        return f'right-ascension: {self._fmtLongitude()}, declination: {self._fmtLatitude()}'

    def __repr__(self):
        """Returns a string representation of the coordinates."""
        return f'CelestialCoordinates({self._fmtLongitude()}, {self._fmtLatitude()})'

    def __reduce__(self):
        """Allows the coordinate object to be pickled."""
        return self.__class__, (self._fmtLatitude(), self._fmtLongitude())

    def toJson(self):
        """Returns a string of the coordinate in json format."""
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        """Returns a dictionary of the coordinate to create json formats of other types containing a GeoPosition."""
        return {"right-ascension": self._fmtLongitude(), "declination": self._fmtLatitude()}

    def _fmod(self, value: float) -> float:
        """Modulates the right-ascension between 0 and 24 hours."""
        return value % 24.0

    def _checkLatitude(self, declination):
        if not isinstance(declination, (int, float)):
            raise TypeError('declination parameter must be an int or float type', type(declination))
        if not -90 <= declination <= 90:
            raise ValueError('declination argument must be between -90 and 90 degrees', declination)
        return _math.radians(declination)

    def _checkLongitude(self, rightAscension):
        if not isinstance(rightAscension, (int, float)):
            raise TypeError('rightAscension parameter must be an int or float type', type(rightAscension))
        return self._fmod(rightAscension) / _RAD_TO_HOURS

    def _fmtLatitude(self):
        """Formats the declination to degrees for printing."""
        return _math.degrees(self._lat)

    def _fmtLongitude(self):
        """Formats the right-ascension to degrees for printing."""
        return self._lng * _RAD_TO_HOURS


def getSubPoint(position: Vector, jd: JulianDate) -> GeoPosition:
    """Computes the geo-position directly below a satellite at a given time.
    """
    declination = _math.degrees(_math.asin(position[2] / position.mag()))
    # longitude is equal to right-ascension minus earth's offset angle at the time
    longitude = (_math.degrees(atan3(position[1], position[0]) - earthOffsetAngle(jd))) % 360.0
    if longitude > 180.0:
        longitude = longitude - 360.0
    # declination is the same angle as geocentric latitude
    return GeoPosition(geocentricToGeodetic(declination), longitude)
