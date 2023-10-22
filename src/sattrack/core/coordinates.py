import json
from math import sqrt, cos, sin
from abc import ABC, abstractmethod
from math import tan, pi, asin, radians, atan, degrees

from pyevspace import Vector, norm

from sattrack.util.constants import EARTH_FLATTENING, TWOPI, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS,\
    HOURS_TO_RAD, EARTH_SIDEREAL_PERIOD
from sattrack.core.juliandate import JulianDate
from sattrack.core.sidereal import earthOffsetAngle
from sattrack.util.helpers import atan3


class Coordinates(ABC):
    __slots__ = '_lat', '_lng', '_latOther', '_lngOther'

    def __init__(self, latitude: float, longitude: float):
        """Angles in degrees, elevation in kilometers."""

        if -90 > latitude > 90:
            raise ValueError(f'latitude must be between -90 and 90 degrees, not {latitude}')

        # latKey and lngKey map the latitude and longitude arguments to radians.
        self._lat = self._latKey(latitude)
        self._lng = self._lngKey(longitude)
        self._latOther = latitude
        self._lngOther = longitude

    @abstractmethod
    def _latKey(self, value: float) -> float:
        pass

    @abstractmethod
    def _lngKey(self, value: float) -> float:
        pass

    def __str__(self):
        return f'latitude: {self._latOther}, longitude: {self._lngOther}'

    def __repr__(self):
        return f'{self.__class__.__name__}({self._latOther}, {self._lngOther})'

    def toDict(self):
        return {"latitude": self._lat, "longitude": self._lng}

    def toJson(self):
        return json.dumps(self, default=lambda o: o.toDict())

    @property
    def latitude(self):
        return self._latOther

    @property
    def longitude(self):
        return self._lngOther

    @property
    def latitudeRadians(self):
        return self._lat

    @property
    def longitudeRadians(self):
        return self._lng


# todo: make a base class for body like coordinates that implement the vector computations
class GeoPosition(Coordinates):
    __slots__ = '_latGeocentric', '_radius', '_elv'

    def __init__(self, latitude: float, longitude: float, elevation: float = 0.0):

        super().__init__(latitude, longitude)

        self._latGeocentric = _geodeticToGeocentric(self._lat)
        self._radius = self._computeRadius()
        self._elv = elevation

    def _latKey(self, value: float) -> float:
        if -90 > value > 90:
            raise ValueError(f'latitude must be between -90 and 90, not {value}')

        return radians(value)

    def _lngKey(self, value: float) -> float:
        return radians(value)

    def __str__(self) -> str:
        return super().__str__() + f', elevation: {self._elv}'

    def __repr__(self) -> str:
        # Strip the closing parenthesis.
        return super().__repr__()[:-1] + f', {self._elv})'

    def toDict(self) -> dict:
        tmp = super().toDict()
        tmp["elevation"] = self._elv

        return tmp

    def toJson(self) -> str:
        return json.dumps(self, default=lambda o: o.toDict())

    @property
    def latitudeGeocentric(self) -> float:
        return self._latGeocentric

    @property
    def elevation(self) -> float:
        return self._elv

    @property
    def radius(self):
        return self._radius

    def _computeRadius(self) -> float:
        piOverTwo = pi / 2
        if self._lat == 0:
            return EARTH_EQUITORIAL_RADIUS
        elif self._lat == piOverTwo or self._lat == -piOverTwo:
            return EARTH_POLAR_RADIUS

        bCos = (EARTH_POLAR_RADIUS * cos(self._latGeocentric))
        aSin = (EARTH_EQUITORIAL_RADIUS * sin(self._latGeocentric))

        numerator = EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS
        denominator = sqrt(bCos * bCos + aSin * aSin)

        return numerator / denominator

    def computeVelocityVector(self) -> Vector:
        """Topocentric units in kilometers per second."""

        velocity = self._radius * TWOPI / EARTH_SIDEREAL_PERIOD
        return Vector(0, velocity, 0)

    def getPositionVector(self, time: JulianDate) -> Vector:
        return _computeNormalVector(self._latGeocentric, self._lng, self._radius, time)

    def getZenithVector(self, time: JulianDate) -> Vector:
        normalVector = _computeNormalVector(self._lat, self._lng, self._radius, time)
        return norm(normalVector)


class CelestialCoordinates(Coordinates):
    def __init__(self, rightAscension: float, declination: float):
        super().__init__(rightAscension, declination)

    def __str__(self) -> str:
        return f'right-ascension: {self._lngOther}, declination: {self._latOther}'

    def toDict(self) -> dict:
        return {"right-ascension": self._lngOther, "declination": self._latOther}

    def toJson(self) -> str:
        return json.dumps(self, default=lambda o: o.toDict())

    def _latKey(self, value: float) -> float:
        if -90 > value > 90:
            raise ValueError(f'declination must be between -90 and 90, not {value}')

        return radians(value)

    def _lngKey(self, value: float) -> float:
        return value * HOURS_TO_RAD

    @property
    def rightAscension(self):
        return self._lngOther

    @property
    def declination(self):
        return self._latOther

    @property
    def rightAscensionRadians(self):
        return self._lng

    @property
    def declinationRadians(self):
        return self._lat


def geocentricToGeodetic(geocentricLatitude: float) -> float:
    """Converts a geocentric latitude in degrees to a geodetic latitude in degrees."""

    latRadians = radians(geocentricLatitude)
    if -90 > geocentricLatitude > 90:
        raise ValueError(f'lat argument must be between -90 and 90 degrees, not {geocentricLatitude}')

    numerator = tan(latRadians)
    tmp = 1 - EARTH_FLATTENING
    denominator = tmp * tmp

    return degrees(atan(numerator / denominator))


def geodeticToGeocentric(geodeticLatitude: float) -> float:
    """Converts a geodetic latitude in degrees to a geocentric latitude in degrees."""

    return degrees(_geodeticToGeocentric(radians(geodeticLatitude)))


def getSubPoint(position: Vector, jd: JulianDate) -> GeoPosition:
    """Computes the geo-position directly below a satellite at a given time."""

    declination = degrees(asin(position[2] / position.mag()))

    # longitude is equal to right-ascension minus earth's offset angle at the time
    longitude = (degrees(atan3(position[1], position[0]) - earthOffsetAngle(jd))) % 360.0
    if longitude > 180.0:
        longitude = longitude - 360.0

    # declination is the same angle as geocentric latitude
    return GeoPosition(geocentricToGeodetic(declination), longitude)


def _computeNormalVector(latitude: float, longitude: float, radius: float, time: JulianDate) -> Vector:
    """Logic for computing surface vectors determined by latitude type. Angles in radians and radius in kilometers."""
    longitude += earthOffsetAngle(time)
    return Vector(
        radius * cos(latitude) * cos(longitude),
        radius * cos(latitude) * sin(longitude),
        radius * sin(latitude)
    )


def _geodeticToGeocentric(geodeticLatitude):
    """Convert a geodetic latitude in radians to geocentric latitude in radians."""
    term_1 = 1 - EARTH_FLATTENING
    arg = (term_1 * term_1) * tan(geodeticLatitude)
    return atan(arg)
