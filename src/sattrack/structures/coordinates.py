import json as _json
import abc as _abc
import math as _math
from sattrack.util.constants import EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS
from pyevspace import EVector
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.conversions import atan3

__all__ = ['Coordinates', 'GeoPosition', 'CelestialCoordinates', 'radiusAtLatitude', 'geocentricToGeodetic',
           'geodeticToGeocentric', 'geoPositionVector', 'zenithVector', 'getSubPoint']

_RAD_TO_HOURS = 12 / _math.pi
_DEG_TO_HOURS = 15


class Coordinates(_abc.ABC):
    """Base class for all coordinate type classes."""

    __slots__ = '_lat', '_lng'

    def __init__(self, lat: float, lng: float):
        """Initializes the object to the latitude and longitude given in degrees."""
        self._lat = self._check_lat(lat)
        self._lng = self._check_lng(lng)

    def __str__(self) -> str:
        """Returns a string representation of the object."""
        return f'Latitude: {self._fmt_lat()}, Longitude: {self._fmt_lng()}'

    def __repr__(self) -> str:
        return _json.dumps(self, default=lambda o: dict(o))

    @property
    def latitude(self):
        return self._fmt_lat()

    @latitude.setter
    def latitude(self, value):
        self._lat = self._check_lat(value)

    @property
    def longitude(self):
        return self._fmt_lng()

    @longitude.setter
    def longitude(self, value):
        self._lng = self._check_lng(value)

    @_abc.abstractmethod
    def _fmod(self, val: float) -> float:
        """Abstract method to be overridden in derived classes, which modulates the longitude value when set."""
        pass

    @_abc.abstractmethod
    def _check_lat(self, lat):
        pass

    @_abc.abstractmethod
    def _check_lng(self, lng):
        pass

    @_abc.abstractmethod
    def _fmt_lat(self):
        pass

    @_abc.abstractmethod
    def _fmt_lng(self):
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

    def __init__(self, lat: float, lng: float, elevation: float = 0):
        """
        Initializes the object with the given values.

        Args:
            lat: The latitude in degrees.
            lng: The longitude in degrees.
            elevation: The elevation of the GeoPosition in km (Default = 0).
        """
        super().__init__(lat, lng)
        self._elevation = elevation

    def __iter__(self):
        yield from {
            'latitude': self._lat,
            'longitude': self._lng,
            'elevation': self._elevation
        }.items()

    def __reduce__(self):
        return self.__class__, (_math.degrees(self._lat), _math.degrees(self._lng), self._elevation)

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        self._elevation = value

    def getRadius(self) -> float:
        """
        Returns the radius of the GeoPosition from the Earth's center in kilometers. This accounts for Earth's
        oblateness as well as the elevation of the position (if set).
        """
        return _radius_at_lat(self._lat) + self._elevation

    def getGeocentricLatitude(self) -> float:
        return _geodetic_to_geocentric(self._lat)

    def getPositionVector(self, jd: JulianDate = None) -> EVector:
        if jd is not None and not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type')
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

    def _fmod(self, val: float) -> float:
        """Modulates the longitude to a value of +/- 180.0 degrees."""
        if val > 180.0:
            tmp = val - int(val / 360.0) * 360.0
            return tmp - 360.0 if tmp > 180.0 else tmp
        elif val < -180.0:
            tmp = val - int(val / 360.0) * 360.0
            return tmp + 360.0 if tmp < -180.0 else tmp
        else:
            return val

    def _check_lat(self, lat):
        if not isinstance(lat, (int, float)):
            raise TypeError('lat parameter must be an int or float type')
        if not -90 <= lat <= 90:
            raise ValueError('lat argument must be between -90 and 90 degrees')
        return _math.radians(lat)

    def _check_lng(self, lng):
        if not isinstance(lng, (int, float)):
            raise TypeError('lng parameter must be an int or float type')
        return _math.radians(self._fmod(lng))

    def _fmt_lat(self):
        return _math.degrees(self._lat)

    def _fmt_lng(self):
        return _math.degrees(self._lng)


def _geocentric_to_geodetic(lat):
    """radians to radians"""
    if not -_math.pi / 2 <= lat <= _math.pi / 2:
        raise ValueError('lat argument must be between -90 and 90 degrees')
    return _math.atan(_math.tan(lat) / ((1 - EARTH_FLATTENING) ** 2))


def _geodetic_to_geocentric(lat):
    """radians to radians"""
    if not -_math.pi / 2 <= lat <= _math.pi / 2:
        raise ValueError('lat argument must be between -90 and 90 degrees')
    return _math.atan(((1 - EARTH_FLATTENING) ** 2) * _math.tan(lat))


def geocentricToGeodetic(lat: float) -> float:
    """Converts a geocentric latitude to a geodetic latitude.
    (lat, return) in degrees"""
    return _math.degrees(_geocentric_to_geodetic(_math.radians(lat)))


def geodeticToGeocentric(lat: float) -> float:
    """Converts a geodetic latitude to a geocentric latitude.
    (lat, return) in degrees"""
    return _math.degrees(_geodetic_to_geocentric(_math.radians(lat)))


def _radius_at_lat(lat):
    """radians to km"""
    if lat == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif lat == 90 or lat == -90:
        return EARTH_POLAR_RADIUS
    geocentric = _geodetic_to_geocentric(lat)  # checks range of lat for us
    bcos2 = (EARTH_POLAR_RADIUS * _math.cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * _math.sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / _math.sqrt(bcos2 + asin2)


def _compute_normal_vector(lat, lng, radius, jd):
    if jd is not None:
        lng += earthOffsetAngle(jd)
    return EVector(
        radius * _math.cos(lat) * _math.cos(lng),
        radius * _math.cos(lat) * _math.sin(lng),
        radius * _math.sin(lat)
    )


def radiusAtLatitude(lat: float) -> float:
    """Computes the radius at a given geodetic latitude.
    lat in degrees, return in kilometers"""
    return _radius_at_lat(_math.radians(lat))


def geoPositionVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    """Computes the position vector of a geo-position at a given time.
    jd is None to ignore earth offset (local sidereal time = 0), return in kilometers"""
    lat = _math.radians(geo.latitude)
    return _compute_normal_vector(
        _geodetic_to_geocentric(lat),
        _math.radians(geo.longitude),
        _radius_at_lat(lat),
        jd
    )


def zenithVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    """Computes the zenith vector of a geo-position at a given time.
    jd is None to ignore earth offset angle (local sidereal time = 0), return in kilometers"""
    lat = _math.radians(geo.latitude)
    return _compute_normal_vector(
        lat,
        _math.radians(geo.longitude),
        _radius_at_lat(lat),
        jd
    )


class CelestialCoordinates(Coordinates):
    """
    A container clas for a set of celestial coordinates. The coordinate names are declination and right-ascension.
    Declination is the equivalent of latitude, and is the angle between the celestial equator and a point in space.
    Right-ascension is the equivalent of longitude, and is the angle from the first point of Ares, or the x-axis. It is
    measured in time units, from 0 to 24 hours, and is positive to the east.
    """

    def __init__(self, ra: float, dec: float):
        """
        Initializes the values to the arguments given.

        Args:
            ra: Right-ascension of the celestial position measured in time units.
            dec: Declination of the celestial position measured in degrees.
        """

        super().__init__(dec, ra)

    def __iter__(self):
        yield from {
            'right-ascension': self._lng,
            'declination': self._lat
        }.items()

    def __reduce__(self):
        return self.__class__, (_math.degrees(self._lat), _math.degrees(self._lng))

    def _fmod(self, val: float) -> float:
        """Modulates the right-ascension value to between (0, 24}."""
        return val % 24.0

    def _check_lat(self, dec):
        if not isinstance(dec, (int, float)):
            raise TypeError('lat parameter must be an int or float type')
        if not -90 <= dec <= 90:
            raise ValueError('dec argument must be between -90 and 90 degrees')
        return _math.radians(dec)

    def _check_lng(self, ra):
        if not isinstance(ra, (int, float)):
            raise TypeError('ra parameter must be an int or float type')
        return self._fmod(ra) / _RAD_TO_HOURS

    def _fmt_lat(self):
        return _math.degrees(self._lat)

    def _fmt_lng(self):
        return self._lng * _RAD_TO_HOURS


def getSubPoint(position: EVector, time: JulianDate) -> GeoPosition:
    """
    Computes the geo-position directly under a satellite.

    Args:
        position: Position vector of the satellite.
        time: Time the satellite is at the position vector.

    Returns:
        The GeoPosition directly below the satellite.
    """

    dec = _math.degrees(_math.asin(position[2] / position.mag()))
    lng = (_math.degrees(atan3(position[1], position[0]) - earthOffsetAngle(time))) % 360.0
    if lng > 180.0:
        lng = lng - 360.0
    return GeoPosition(geocentricToGeodetic(dec), lng)
