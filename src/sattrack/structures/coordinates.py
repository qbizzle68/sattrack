import json as _json
import abc as _abc
import math as _math
from sattrack.util.constants import EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS
from pyevspace import EVector
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.conversions import atan3

__all__ = ['Coordinates', 'GeoPosition', 'CelestialCoordinates', 'radiusAtLatitude', 'geocentricToGeodetic',
           'geodeticToGeocentric', 'getSubPoint']

_RAD_TO_HOURS = 12 / _math.pi
_DEG_TO_HOURS = 15


class Coordinates(_abc.ABC):
    """Base class for all coordinate type classes."""

    __slots__ = '_lat', '_lng'

    def __init__(self, latitude: float, longitude: float):
        # constructor
        self._lat = self._check_latitude(latitude)
        self._lng = self._check_longitude(longitude)

    def __repr__(self) -> str:
        return _json.dumps(self, default=lambda o: dict(o))

    @property
    def latitude(self):
        return self._fmt_latitude()

    @latitude.setter
    def latitude(self, value):
        self._lat = self._check_latitude(value)

    @property
    def longitude(self):
        return self._fmt_longitude()

    @longitude.setter
    def longitude(self, value):
        self._lng = self._check_longitude(value)

    # abstract methods to be overridden
    @_abc.abstractmethod
    def _fmod(self, val: float) -> float:
        pass

    @_abc.abstractmethod
    def _check_latitude(self, lat):
        pass

    # should call self._fmod() before returning a value
    @_abc.abstractmethod
    def _check_longitude(self, lng):
        pass

    @_abc.abstractmethod
    def _fmt_latitude(self):
        pass

    @_abc.abstractmethod
    def _fmt_longitude(self):
        pass


class GeoPosition(Coordinates):
    """
    A container class for an Earth geo-position. The class contains a latitude, longitude and an optional elevation. The
    elevation is used for positions relative to a geo-position, as well as rotating to a topocentric reference frame.
    Since the Earth's shape is approximated as an oblate sphere, there are two types of latitudes, geocentric and
    geodetic. Geocentric latitude is the angle between the equator and a point on Earth's surface through the Earth's
    center. Geodetic latitude is the angle between the equator and a line normal to Earth's surface at a given point.
    The geodetic latitude is also known as geographic latitude, and is the one used in a geo-position coordinate.
    Latitudes are measured zero at the equator, positive to the north, and negative to the south. Longitude is measured
    zero on the longitude passing through Greenwich, England, and is positive to the east, and negative to the west.
    """

    __slots__ = '_elevation',

    def __init__(self, latitude: float, longitude: float, elevation: float = 0):
        """
        Initializes the object with the given values.

        Args:
            latitude: The latitude in degrees.
            longitude: The longitude in degrees.
            elevation: The elevation of the GeoPosition in km (Default = 0).
        """
        super().__init__(latitude, longitude)
        if not isinstance(elevation, (int, float)):
            raise TypeError('elevation parameter must be an int or float type')
        self._elevation = elevation

    def __iter__(self):
        yield from {
            'latitude': self._lat,
            'longitude': self._lng,
            'elevation': self._elevation
        }.items()

    def __str__(self):
        return f'latitude: {self._fmt_latitude()}, longitude: {self._fmt_longitude()}, elevation: {self._elevation}'

    def __reduce__(self):
        return self.__class__, (self._fmt_latitude(), self._fmt_longitude(), self._elevation)

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        self._elevation = value

    def getRadius(self) -> float:
        return _radius_at_lat(self._lat) + self._elevation

    def getGeocentricLatitude(self) -> float:
        return _geodetic_to_geocentric(self._lat)

    def getPositionVector(self, jd: JulianDate = None) -> EVector:
        # if jd is not None and not isinstance(jd, JulianDate):
        #     raise TypeError('jd parameter must be a JulianDate type')
        # position vector needs geocentric latitude
        return _compute_normal_vector(
            _geodetic_to_geocentric(self._lat),
            self._lng,
            _radius_at_lat(self._lat) + self._elevation,
            jd
        )

    def getZenithVector(self, jd: JulianDate = None) -> EVector:
        if jd is not None and not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type')
        return _compute_normal_vector(
            self._lat,
            self._lng,
            _radius_at_lat(self._lat) + self._elevation,
            jd
        )

    def _fmod(self, value: float) -> float:
        # find the equivalent angle from -180 to 180
        if value > 180.0:
            tmp = value - int(value / 360.0) * 360.0
            return tmp - 360.0 if tmp > 180.0 else tmp
        elif value < -180.0:
            tmp = value - int(value / 360.0) * 360.0
            return tmp + 360.0 if tmp < -180.0 else tmp
        else:
            return value

    def _check_latitude(self, latitude):
        if not isinstance(latitude, (int, float)):
            raise TypeError('lat parameter must be an int or float type')
        if not -90 <= latitude <= 90:
            raise ValueError('lat argument must be between -90 and 90 degrees')
        return _math.radians(latitude)

    def _check_longitude(self, longitude):
        if not isinstance(longitude, (int, float)):
            raise TypeError('lng parameter must be an int or float type')
        return _math.radians(self._fmod(longitude))

    def _fmt_latitude(self):
        return _math.degrees(self._lat)

    def _fmt_longitude(self):
        return _math.degrees(self._lng)


def _geocentric_to_geodetic(geocentricLatitude):
    # parameter and return value are in radians
    if not -_math.pi / 2 <= geocentricLatitude <= _math.pi / 2:
        raise ValueError('lat argument must be between -90 and 90 degrees')
    return _math.atan(_math.tan(geocentricLatitude) / ((1 - EARTH_FLATTENING) ** 2))


def _geodetic_to_geocentric(geodeticLatitude):
    # parameter and return value are in radians
    if not -_math.pi / 2 <= geodeticLatitude <= _math.pi / 2:
        raise ValueError('lat argument must be between -90 and 90 degrees')
    return _math.atan(((1 - EARTH_FLATTENING) ** 2) * _math.tan(geodeticLatitude))


def geocentricToGeodetic(geocentricLatitude: float) -> float:
    # parameter and return value are in degrees
    return _math.degrees(_geocentric_to_geodetic(_math.radians(geocentricLatitude)))


def geodeticToGeocentric(geodeticLatitude: float) -> float:
    # parameter and return value are in degrees
    return _math.degrees(_geodetic_to_geocentric(_math.radians(geodeticLatitude)))


def _radius_at_lat(latitude):
    # latitude refers to geodetic latitude, returns radius in kilometers
    if latitude == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif latitude == 90 or latitude == -90:
        return EARTH_POLAR_RADIUS

    # _geodetic_to_geocentric() checks latitude range for us
    geocentric = _geodetic_to_geocentric(latitude)
    # radius = (Rp * Re) / sqrt((Rp*cos(geocentric))^2 + (Re*sin(geocentric))^2)
    bcos2 = (EARTH_POLAR_RADIUS * _math.cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * _math.sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / _math.sqrt(bcos2 + asin2)


def _compute_normal_vector(latitude, longitude, radius, jd):
    # only account for earth offset when jd is not None
    if jd is not None:
        if isinstance(jd, JulianDate):
            longitude += earthOffsetAngle(jd)
        else:
            raise TypeError('jd parameter must be a JulianDate type')
    return EVector(
        radius * _math.cos(latitude) * _math.cos(longitude),
        radius * _math.cos(latitude) * _math.sin(longitude),
        radius * _math.sin(latitude)
    )


def radiusAtLatitude(latitude: float) -> float:
    # latitude in degrees, return value in kilometers
    return _radius_at_lat(_math.radians(latitude))


class CelestialCoordinates(Coordinates):
    """
    A container clas for a set of celestial coordinates. The coordinate names are declination and right-ascension.
    Declination is the equivalent of latitude, and is the angle between the celestial equator and a point in space.
    Right-ascension is the equivalent of longitude, and is the angle from the first point of Ares, or the x-axis. It is
    measured in time units, from 0 to 24 hours, and is positive to the east.
    """

    def __init__(self, rightAscension: float, declination: float):
        # order of parameters goes with the idiom 'ra, dec'
        super().__init__(declination, rightAscension)

    def __iter__(self):
        yield from {
            'right-ascension': self._lng,
            'declination': self._lat
        }.items()

    def __str__(self):
        return f'right-ascension: {self._fmt_longitude()}, declination: {self._fmt_latitude()}'

    def __reduce__(self):
        # return self.__class__, (_math.degrees(self._lat), _math.degrees(self._lng))
        return self.__class__, (self._fmt_latitude(), self._fmt_longitude())

    def _fmod(self, value: float) -> float:
        # keeps longitude between 0 and 24 hours
        return value % 24.0

    def _check_latitude(self, declination):
        if not isinstance(declination, (int, float)):
            raise TypeError('declination parameter must be an int or float type')
        if not -90 <= declination <= 90:
            raise ValueError('declination argument must be between -90 and 90 degrees')
        return _math.radians(declination)

    def _check_longitude(self, rightAscension):
        if not isinstance(rightAscension, (int, float)):
            raise TypeError('rightAscension parameter must be an int or float type')
        return self._fmod(rightAscension) / _RAD_TO_HOURS

    def _fmt_latitude(self):
        return _math.degrees(self._lat)

    def _fmt_longitude(self):
        return self._lng * _RAD_TO_HOURS


def getSubPoint(position: EVector, jd: JulianDate) -> GeoPosition:
    declination = _math.degrees(_math.asin(position[2] / position.mag()))
    # longitude is equal to right-ascension minus earth's offset angle at the time
    longitude = (_math.degrees(atan3(position[1], position[0]) - earthOffsetAngle(jd))) % 360.0
    if longitude > 180.0:
        longitude = longitude - 360.0
    # declination is the same angle as geocentric latitude
    return GeoPosition(geocentricToGeodetic(declination), longitude)
