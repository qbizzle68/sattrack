import math as _math

import pyevspace as evs
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.constants import EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS, EARTH_FLATTENING


def _computeZenithVector(latitude: float, longitude: float, elevation: float, time) -> evs.Vector:
    """Compute a normal vector to earth surface (latitude is geodetic in radians). The magnitude of the zenith vector
    is guaranteed to have a magnitude of 1."""
    radius = _radiusAtLat(latitude)
    return evs.norm(_computeNormalVector(latitude, longitude, radius + elevation, time))


def _computePositionVector(latitude: float, longitude: float, elevation: float, time) -> evs.Vector:
    """Compute the position vector of a geo-position (latitude is geodetic in radians)."""
    radius = _radiusAtLat(latitude)
    geocentricLatitude = _geodeticToGeocentric(latitude)
    return _computeNormalVector(geocentricLatitude, longitude, radius + elevation, time)


def _radiusAtLat(latitude: float):
    """Returns the Earth radius in kilometers from a geodetic latitude in radians."""
    piOver2 = _math.pi / 2
    if latitude == 0:
        return EARTH_EQUITORIAL_RADIUS
    elif latitude == piOver2 or latitude == -piOver2:
        return EARTH_POLAR_RADIUS

    geocentric = _geodeticToGeocentric(latitude)
    # radius = (Rp * Re) / sqrt((Rp*cos(geocentric))^2 + (Re*sin(geocentric))^2)
    bcos2 = (EARTH_POLAR_RADIUS * _math.cos(geocentric)) ** 2
    asin2 = (EARTH_EQUITORIAL_RADIUS * _math.sin(geocentric)) ** 2
    return (EARTH_POLAR_RADIUS * EARTH_EQUITORIAL_RADIUS) / _math.sqrt(bcos2 + asin2)


def _computeNormalVector(latitude: float, longitude: float, radius: float, time):
    """Logic for computing surface vectors determined by latitude type. Angles in radians and radius in kilometers."""
    longitude += earthOffsetAngle(time)
    return evs.Vector((
        radius * _math.cos(latitude) * _math.cos(longitude),
        radius * _math.cos(latitude) * _math.sin(longitude),
        radius * _math.sin(latitude)
    ))


def _geodeticToGeocentric(geodeticLatitude):
    """Convert a geodetic latitude in radians to geocentric latitude in radians."""
    term_1 = 1 - EARTH_FLATTENING
    arg = (term_1 * term_1) * _math.tan(geodeticLatitude)
    return _math.atan(arg)
