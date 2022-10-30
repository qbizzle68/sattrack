import math as _math

from pyevspace import EVector
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.constants import EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS, EARTH_FLATTENING


def _compute_zenith_vector(latitude: float, longitude: float, elevation: float, time) -> EVector:
    # geodetic latitude
    radius = _radius_at_lat(latitude)
    return _compute_normal_vector(latitude, longitude, radius + elevation, time)


def _compute_position_vector(latitude: float, longitude: float, elevation: float, time) -> EVector:
    # geodetic latitude
    radius = _radius_at_lat(latitude)
    geocentricLatitude = _geodetic_to_geocentric(latitude)
    return _compute_normal_vector(geocentricLatitude, longitude, radius + elevation, time)


def _radius_at_lat(latitude):
    # latitude refers to geodetic latitude, returns radius in kilometers
    if latitude == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif latitude == 90 or latitude == -90:
        return EARTH_POLAR_RADIUS

    geocentric = _geodetic_to_geocentric(latitude)
    # radius = (Rp * Re) / sqrt((Rp*cos(geocentric))^2 + (Re*sin(geocentric))^2)
    bcos2 = (EARTH_POLAR_RADIUS * _math.cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * _math.sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / _math.sqrt(bcos2 + asin2)


def _compute_normal_vector(latitude, longitude, radius, time):
    longitude += earthOffsetAngle(time)
    return EVector(
        radius * _math.cos(latitude) * _math.cos(longitude),
        radius * _math.cos(latitude) * _math.sin(longitude),
        radius * _math.sin(latitude)
    )


def _geodetic_to_geocentric(geodeticLatitude):
    # parameter and return value are in radians
    term_1 = 1 - EARTH_FLATTENING
    arg = (term_1 * term_1) * _math.tan(geodeticLatitude)
    return _math.atan(arg)
