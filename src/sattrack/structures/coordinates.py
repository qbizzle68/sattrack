import json as _json
from abc import ABC, abstractmethod
from math import sin, cos, tan, atan, radians, degrees, sqrt, asin
from sattrack.util.constants import EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS
from pyevspace import EVector
from sattrack.rotation.rotation import rotateOrderFrom, EulerAngles
from sattrack.rotation.order import ZYX
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.conversions import atan3


class Coordinates(ABC):
    """Base class for all coordinate type classes."""

    def __init__(self, lat: float, lng: float):
        """Initializes the object to the latitude and longitude given in degrees."""
        self._check_angles(lat, lng)
        self._lat = lat
        self._lng = self._fmod(lng)

    def __str__(self) -> str:
        """Returns a string representation of the object."""
        return f'Latitude: {self._lat}, Longitude: {self._lng}'

    def __repr__(self) -> str:
        return _json.dumps(self, default=lambda o: dict(o))

    # def toJSON(self):
    #     return _json.dumps(self, indent=4, default=lambda o: o.__dict__)

    def getLatitude(self) -> float:
        """Returns the latitude in degrees."""
        return self._lat

    def setLatitude(self, lat: float) -> None:
        """Sets the latitude to the given value in degrees."""
        self._lat = lat

    def getLongitude(self) -> float:
        """Returns the longitude in degrees."""
        return self._lng

    def setLongitude(self, lng: float) -> None:
        """Sets the longitude to the given value in degrees."""
        self._lng = self._fmod(lng)

    @abstractmethod
    def _fmod(self, val: float) -> float:
        """Abstract method to be overridden in derived classes, which modulates the longitude value when set."""
        pass

    def _check_angles(self, lat, lng):
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

    def __init__(self, lat: float, lng: float, elevation: float = 0):
        """
        Initializes the object with the given values.

        Args:
            lat: The latitude in degrees.
            lng: The longitude in degrees.
            elevation: The elevation of the GeoPosition in meters (Default = 0).
        """
        super().__init__(lat, lng)
        self._elevation = elevation

    @classmethod
    def fromJSON(cls, json):
        return cls(json['_lat'], json['_lng'], json['_elevation'])

    def getElevation(self) -> float:
        """Returns the elevation in meters."""
        return self._elevation

    def setElevation(self, elevation: float) -> None:
        """Sets the elevation to the given value in meters."""
        self._elevation = elevation

    def getRadius(self) -> float:
        """
        Returns the radius of the GeoPosition from the Earth's center in kilometers. This accounts for Earth's
        oblateness as well as the elevation of the position (if set).
        """

        return radiusAtLatitude(self._lat) + self._elevation

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


def geocentricToGeodetic(lat: float) -> float:
    """
    Converts a geocentric latitude to a geodetic latitude.

    Args:
        lat: Geocentric latitude in degrees.

    Returns:
        Geodetic latitude of the equivalent position in degrees.
    """

    lhs = tan(radians(lat)) / (1 - EARTH_FLATTENING) ** 2
    return degrees(atan(lhs))


def geodeticToGeocentric(lat: float) -> float:
    """
    Converts a geodetic latitude to a geocentric latitude.

    Args:
        lat: Geodetic latitude in degrees.

    Returns:
        Geocentric latitude of the equivalent position in degrees.
    """
    rhs = ((1 - EARTH_FLATTENING) ** 2) * tan(radians(lat))
    return degrees(atan(rhs))


def radiusAtLatitude(lat: float) -> float:
    """
    Computes the radius at a given latitude due to the oblate shape of the Earth.

    Args:
        lat: The geodetic latitude in degrees.

    Returns:
        The Earth radius at the given geodetic latitude in kilometers.
    """

    if lat == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif lat == 90:
        return EARTH_POLAR_RADIUS
    geocentric = radians(geodeticToGeocentric(lat))
    bcos2 = (EARTH_POLAR_RADIUS * cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / sqrt(bcos2 + asin2)


def geoPositionVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    """
    Computes the position vector in the geocentric equitorial reference frame, of a GeoPosition at a given time.

    Args:
        geo: GeoPosition of the position to find.
        jd: Time of the desired position (used to incorporate the Earth offset). (Default = None, won't include the
        Earth's offset from the inertial reference frame).

    Returns:
        The geocentric equitorial reference frame position vector of the GeoPosition.
    """

    if not jd:
        return rotateOrderFrom(
            ZYX,
            EulerAngles(
                radians(geo.getLatitude()),
                radians(-geodeticToGeocentric(geo.getLatitude())),
                0.0),
            EVector(
                # radians(radiusAtLatitude(geo.getLatitude())),
                radiusAtLatitude(geo.getLatitude()) + geo.getElevation(),
                0.0,
                0.0)
        )
    else:
        radiusAtLat = radiusAtLatitude(geo.getLatitude()) + geo.getElevation()
        geocentricLat = radians(geodeticToGeocentric(geo.getLatitude()))
        lst = radians(geo.getLongitude()) + earthOffsetAngle(jd)
        return EVector(
            radiusAtLat * cos(geocentricLat) * cos(lst),
            radiusAtLat * cos(geocentricLat) * sin(lst),
            radiusAtLat * sin(geocentricLat)
        )


def zenithVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    """
    Computes the zenith vector in the geocentric equitorial reference frame, of a GeoPosition at a given time.

    Args:
        geo: GeoPosition of the zenith to find.
        jd: Time of the desired zenith (used to incorporate the Earth offset). (Default = None, won't include the
        Earth's offset from the inertial reference frame).

    Returns:
        The geocentric equitorial reference frame vector of the zenith for the GeoPosition.
    """

    if not jd:
        return rotateOrderFrom(
            ZYX,
            EulerAngles(
                radians(geo.getLatitude()),
                radians(geo.getLatitude()),
                0.0),
            EVector(
                radiusAtLatitude(geo.getLatitude()) + geo.getElevation(),
                0.0,
                0.0)
        )
    else:
        radiusAtLat = radiusAtLatitude(geo.getLatitude()) + geo.getElevation()
        lat = radians(geo.getLatitude())
        lst = radians(geo.getLongitude()) + earthOffsetAngle(jd)
        return EVector(
            radiusAtLat * cos(lat) * cos(lst),
            radiusAtLat * cos(lat) * sin(lst),
            radiusAtLat * sin(lat)
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
        self.setRightAscension = super().setLongitude
        self.getRightAscension = super().getLongitude
        self.setDeclination = super().setLatitude
        self.getDeclination = super().getLatitude

    def __str__(self) -> str:
        """Returns a string representation of the object."""
        return f'Right-Ascension: {self._lng}, Declination: {self._lat}'

    def getLongitude(self) -> float:
        """Returns the right-ascension as a longitude in degrees."""
        return self._lng * 15.0

    def _fmod(self, val: float) -> float:
        """Modulates the right-ascension value to between (0, 24}."""
        return val % 24.0


def getSubPoint(position: EVector, time: JulianDate) -> GeoPosition:
    """
    Computes the geo-position directly under a satellite.

    Args:
        position: Position vector of the satellite.
        time: Time the satellite is at the position vector.

    Returns:
        The GeoPosition directly below the satellite.
    """

    dec = degrees(asin(position[2] / position.mag()))
    lng = (degrees(atan3(position[1], position[0]) - earthOffsetAngle(time))) % 360.0
    if lng > 180.0:
        lng = lng - 360.0
    return GeoPosition(geocentricToGeodetic(dec), lng)
