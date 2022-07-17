from abc import ABC, abstractmethod
from math import sin, cos, tan, atan, radians, degrees, sqrt
from sattrack.util.constants import EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS
from pyevspace import EVector
from sattrack.rotation.rotation import rotateOrderFrom, EulerAngles
from sattrack.rotation.order import ZYX
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle


class Coordinates(ABC):

    def __init__(self, lat: float, lng: float):
        self._lat = lat
        self._lng = self._fmod(lng)

    def __str__(self) -> str:
        return f'Latitude: {self._lat}, Longitude: {self._lng}'

    def getLatitude(self) -> float:
        return self._lat

    def setLatitude(self, lat: float) -> None:
        self._lat = lat

    def getLongitude(self) -> float:
        return self._lng

    def setLongitude(self, lng: float) -> None:
        self._lng = self._fmod(lng)

    @abstractmethod
    def _fmod(self, val: float) -> float:
        pass


class GeoPosition(Coordinates):

    def __init__(self, lat: float, lng: float, elevation: float = 0):
        super().__init__(lat, lng)
        self._elevation = elevation

    def getElevation(self) -> float:
        return self._elevation

    def setElevation(self, el: float) -> None:
        self._elevation = el

    def getRadius(self) -> float:
        return radiusAtLatitude(self._lat) + self._elevation

    def _fmod(self, val: float) -> float:
        if val > 180.0:
            tmp = val - int(val / 360.0) * 360.0
            return tmp - 360.0 if tmp > 180.0 else tmp
        elif val < -180.0:
            tmp = val - int(val / 360.0) * 360.0
            return tmp + 360.0 if tmp < -180.0 else tmp
        else:
            return val


def geocentricToGeodetic(lat: float) -> float:
    lhs = tan(radians(lat)) / (1 - EARTH_FLATTENING) ** 2
    return degrees(atan(lhs))


def geodeticToGeocentric(lat: float) -> float:
    rhs = ((1 - EARTH_FLATTENING) ** 2) * tan(radians(lat))
    return degrees(atan(rhs))


def radiusAtLatitude(lat: float) -> float:
    if lat == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif lat == 90:
        return EARTH_POLAR_RADIUS
    geocentric = radians(geodeticToGeocentric(lat))
    bcos2 = (EARTH_POLAR_RADIUS * cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / sqrt(bcos2 + asin2)


def geoPositionVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    if not jd:
        return rotateOrderFrom(
            ZYX,
            EulerAngles(
                radians(geo.getLatitude()),
                radians(-geodeticToGeocentric(geo.getLatitude())),
                0.0),
            EVector(
                radians(radiusAtLatitude(geo.getLatitude())),
                0.0,
                0.0)
        )
    else:
        radiusAtLat = radiusAtLatitude(geo.getLatitude())
        geocentricLat = radians(geodeticToGeocentric(geo.getLatitude()))
        lst = radians(geo.getLongitude()) + earthOffsetAngle(jd)
        return EVector(
            radiusAtLat * cos(geocentricLat) * cos(lst),
            radiusAtLat * cos(geocentricLat) * sin(lst),
            radiusAtLat * sin(geocentricLat)
        )


def zenithVector(geo: GeoPosition, jd: JulianDate = None) -> EVector:
    if not jd:
        return rotateOrderFrom(
            ZYX,
            EulerAngles(
                radians(geo.getLatitude()),
                radians(geo.getLatitude()),
                0.0),
            EVector(
                radiusAtLatitude(geo.getLatitude()),
                0.0,
                0.0)
        )
    else:
        radiusAtLat = radiusAtLatitude(geo.getLatitude())
        lat = radians(geo.getLatitude())
        lst = radians(geo.getLongitude()) + earthOffsetAngle(jd)
        return EVector(
            radiusAtLat * cos(lat) * cos(lst),
            radiusAtLat * cos(lat) * sin(lst),
            radiusAtLat * sin(lat)
        )


class CelestialCoordinates(Coordinates):

    def __init__(self, ra: float, dec: float):
        super().__init__(dec, ra)
        self.setRightAscension = super().setLongitude
        self.getRightAscension = super().getLongitude
        self.setDeclination = super().setLatitude
        self.getDeclination = super().getLatitude

    def __str__(self) -> str:
        return f'Right-Ascension: {self._lng}, Declination: {self._lat}'

    def getLatitude(self) -> float:
        return self._lat * 15.0

    def getLongitude(self) -> float:
        return self._lng * 15.0

    def _fmod(self, val: float) -> float:
        return val % 24.0
