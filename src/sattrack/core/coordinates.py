import json
from math import sqrt, cos, sin
from abc import ABC, abstractmethod
from math import tan, pi, asin, radians, atan, degrees

from pyevspace import Vector, norm, ReferenceFrame, Angles, ZYX

from sattrack.bodies.position import getEarthOffsetAngle
from sattrack.util.constants import EARTH_FLATTENING, TWOPI, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS, \
    HOURS_TO_RAD, EARTH_SIDEREAL_PERIOD
from sattrack.util.helpers import atan3

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate


class Coordinates(ABC):
    __slots__ = '_lat', '_lng', '_latOther', '_lngOther'

    def __init__(self, latitude: float, longitude: float):
        """Angles in degrees, elevation in kilometers."""

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
        return {"latitude": self._latOther, "longitude": self._lngOther}

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

    @property
    def parts(self):
        latParts = computeAngleParts(self._latOther)
        lngParts = computeAngleParts(self._lngOther)

        return latParts, lngParts


# todo: make a base class for body like coordinates that implement the vector computations
class GeoPosition(Coordinates):
    __slots__ = '_latGeocentric', '_radius', '_elv'

    def __init__(self, latitude: float, longitude: float, elevation: float = 0.0):

        super().__init__(latitude, longitude)

        self._latGeocentric = _geodeticToGeocentric(self._lat)
        self._radius = self._computeRadius()
        self._elv = elevation

    def _latKey(self, value: float) -> float:
        if value < -90 or value > 90:
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

    def getPositionVector(self, time: 'JulianDate') -> Vector:
        return _computeNormalVector(self._latGeocentric, self._lng, self._radius, time)

    def getZenithVector(self, time: 'JulianDate') -> Vector:
        normalVector = _computeNormalVector(self._lat, self._lng, self._radius, time)
        return norm(normalVector)

    def getReferenceFrame(self, time: 'JulianDate') -> ReferenceFrame:
        lng = self._lng + getEarthOffsetAngle(time)
        lat = pi / 2 - self._lat
        angles = Angles(lng, lat, 0.0)

        return ReferenceFrame(ZYX, angles)


class CelestialCoordinates(Coordinates):
    # rightAscension (longitude) in hours, declination (latitude) degrees
    def __init__(self, rightAscension: float, declination: float):
        super().__init__(declination, 0.0)

        self._lng = rightAscension * HOURS_TO_RAD
        self._lngOther = rightAscension

    def __str__(self) -> str:
        return f'right-ascension: {self._lngOther}, declination: {self._latOther}'

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._lngOther}, {self._latOther})'

    def toDict(self) -> dict:
        return {"right-ascension": self._lngOther, "declination": self._latOther}

    def _latKey(self, value: float) -> float:
        if value < -90 or value > 90:
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

    @property
    def parts(self):
        raParts = computeAngleParts(self._lngOther)
        decParts = computeAngleParts(self._latOther)

        return raParts, decParts


def geocentricToGeodetic(geocentricLatitude: float) -> float:
    """Converts a geocentric latitude in degrees to a geodetic latitude in degrees."""

    if geocentricLatitude < -90 or geocentricLatitude > 90:
        raise ValueError(f'lat argument must be between -90 and 90 degrees, not {geocentricLatitude}')
    latRadians = radians(geocentricLatitude)

    numerator = tan(latRadians)
    tmp = 1 - EARTH_FLATTENING
    denominator = tmp * tmp

    return degrees(atan(numerator / denominator))


def geodeticToGeocentric(geodeticLatitude: float) -> float:
    """Converts a geodetic latitude in degrees to a geocentric latitude in degrees."""

    if geodeticLatitude < -90 or geodeticLatitude > 90:
        raise ValueError(f'lat argument must be between -90 and 90 degrees, not {geodeticLatitude}')
    latRadians = radians(geodeticLatitude)

    return degrees(_geodeticToGeocentric(latRadians))


def computeSubPoint(position: Vector, jd: 'JulianDate') -> GeoPosition:
    """Computes the geo-position directly below a satellite at a given time."""

    declination = degrees(asin(position[2] / position.mag()))

    # longitude is equal to right-ascension minus earth's offset angle at the time
    longitude = (degrees(atan3(position[1], position[0]) - getEarthOffsetAngle(jd))) % 360.0
    if longitude > 180.0:
        longitude = longitude - 360.0

    # declination is the same angle as geocentric latitude
    return GeoPosition(geocentricToGeodetic(declination), longitude)


def _computeNormalVector(latitude: float, longitude: float, radius: float, time: 'JulianDate') -> Vector:
    """Logic for computing surface vectors determined by latitude type. Angles in radians and radius in kilometers."""

    longitude += getEarthOffsetAngle(time)
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


def computeAngleParts(value: float) -> (float, float, float):
    wholePart = int(value)
    frac = value - wholePart
    if frac < 0:
        frac *= -1

    tmp = frac * 60
    minutesWhole = int(tmp)
    frac = tmp - minutesWhole
    seconds = frac * 60

    return wholePart, minutesWhole, seconds


class AltAz:
    __slots__ = '_altitude', '_azimuth', '_direction'

    def __init__(self, altitude: float, azimuth: float):
        # Args in degrees, will be modded to valid ranges.

        if not -90 <= altitude <= 90:
            raise ValueError(f'altitude must be in (-90, 90), was {altitude}')

        self._altitude = altitude
        self._azimuth = azimuth % 360
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
        return {"altitude": self._altitude, "azimuth": self._azimuth, "direction": self._direction}

    def toJson(self):
        return json.dumps(self, default=lambda o: o.toDict())

    def __str__(self):
        return f'altitude: {self._altitude}, azimuth: {self._azimuth}, direction: {self._direction}'

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
