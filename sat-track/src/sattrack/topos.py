from math import degrees, asin, radians, sin, cos, tan, sqrt

from pyevspace import EVector, cross, dot, norm, transpose

from sattrack.rotation.order import ZYX
from sattrack.structures.coordinates import GeoPosition, geoPositionVector, zenithVector
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateToThenOffset, undoRotateToThenOffset
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.satellite import Satellite
from sattrack.sun import TwilightType, getSunPosition
from sattrack.util.conversions import atan2


def toTopocentric(vec: EVector, time: JulianDate, geo: GeoPosition) -> EVector:
    """
    Transform a vector from a geocentric equitorial reference frame to a topocentric horizontal reference frame.

    Args:
        vec: Vector to transform to the topocentric reference frame.
        time: Time of the position, to accommodate for sidereal time.
        geo: GeoPosition which will be the origin of the topocentric reference frame.

    Returns:
        The transformed vector.
    """

    geoVector = geoPositionVector(geo, time)
    mat = getEulerMatrix(
        ZYX,
        EulerAngles(
            radians(geo.getLongitude()) + earthOffsetAngle(time),
            radians(90 - geo.getLatitude()),
            0.0
        )
    )
    return rotateToThenOffset(mat, geoVector, vec)


def fromTopocentric(vec: EVector, time: JulianDate, geo: GeoPosition) -> EVector:
    """
    Transform a vector from a topocentric horizontal reference frame to a geocentric equitorial reference frame.

    Args:
        vec: Vector to transform to the geocentric reference frame.
        time: Time of the position, to accommodate for sidereal time.
        geo: GeoPosition which is the origin of the topocentric reference frame.

    Returns:
        The transformed vector.
    """

    geoVector = geoPositionVector(geo, time)
    mat = getEulerMatrix(
        ZYX,
        EulerAngles(
            radians(geo.getLongitude()) + earthOffsetAngle(time),
            radians(90 - geo.getLatitude()),
            0.0
        )
    )
    return undoRotateToThenOffset(mat, geoVector, vec)


# todo: think of a better name for this
def getPVector(geo: GeoPosition, position: EVector, velocity: EVector, jd: JulianDate) -> EVector:
    """
    Computes the vector that lies in the orbital plane P, where for v on the horizontal plane of a geo-position whose
    endpoint lies on the intersection of the orbital and horizontal planes, P - v is minimized. This vector approximates
    the position of the highest altitude the satellite could achieve above (or below) the horizon, if the satellite was
    at that position in its orbit.

    Args:
        geo: GeoPosition of the viewer.
        position: Position vector of a state of the satellite.
        velocity: Velocity vector of a state of the satellite.
        jd: Time of the computation.

    Returns:
        The so-called 'P' vector of the orbit and horizontal planes.
    """

    zeta = norm(zenithVector(geo, jd))
    gamma = geoPositionVector(geo, jd)
    lamb = norm(cross(position, velocity))
    x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
    y = -lamb[0] * x / lamb[1]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    return r + v * t


def getToposPosition(satellite: Satellite, time: JulianDate, geo: GeoPosition) -> EVector:
    """
    Returns the position vector of a satellite in the topocentric reference frame.
    Args:
        satellite: Satellite to find the position of.
        time: Time to find the position.
        geo: GeoPosition of the topocentric reference frame.

    Returns:
        The satellite position vector in SEZ coordinates.
    """

    return toTopocentric(satellite.getState(time)[0], time, geo)


def getAltitude(satellite: Satellite, time: JulianDate, geo: GeoPosition) -> float:
    """
    Computes the altitude of a satellite from a given geo-position.

    Args:
        satellite: Satellite to find the position of.
        time: Time to find the position.
        geo: GeoPosition to find the relative altitude of.

    Returns:
        The satellite altitude above (or below) the horizon in degrees.
    """

    sez = toTopocentric(satellite.getState(time)[0], time, geo)
    return degrees(asin(sez[2] / sez.mag()))


def getAzimuth(satellite: Satellite, time: JulianDate, geo: GeoPosition) -> float:
    """
    Computes the azimuth of a satellite from a given geo-position. The azimuth of zero degrees points to the north, and
    the angle increases to the east, or in a clockwise direction.

    Args:
        satellite: Satellite to find the position of.
        time: Time to find the position.
        geo: GeoPosition to find the relative azimuth of.

    Returns:
        The satellite azimuth measured clockwise from north in degrees.
    """

    sez = toTopocentric(satellite.getState(time)[0], time, geo)
    return degrees(atan2(sez[1], -sez[0]))


def azimuthAngleString(azimuth: float) -> str:
    """
    Converts an azimuth angle to a string of the compass abbreviation.

    Args:
        azimuth: Azimuth angle in degrees.

    Returns:
        String of the compass abbreviation.
    """

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


def getTwilightType(time: JulianDate, geo: GeoPosition) -> TwilightType:
    """
    Computes the twilight-type from the angle of the sun.

    Args:
        time: Time to find the twilight-type.
        geo: GeoPosition whose twilight-type is to be computed.

    Returns:
        One of the TwilightType enumerations.
    """

    sunSEZPos = toTopocentric(getSunPosition(time), time, geo)
    sunAngle = degrees(asin(sunSEZPos[2] / sunSEZPos.mag()))
    if sunAngle < -18:
        return TwilightType.Night
    elif sunAngle < -12:
        return TwilightType.Astronomical
    elif sunAngle < -6:
        return TwilightType.Nautical
    elif sunAngle < -(5.0 / 6.0):
        return TwilightType.Civil
    else:
        return TwilightType.Day


def altAzToPosition(altitude: float, azimuth: float, distance: float) -> EVector:
    """
    Converts a pair of altitude and azimuth angles into a position vector in the topocentric reference frame.

    Args:
        altitude: Altitude of the position in degrees.
        azimuth: Azimuth of the object position in degrees.
        distance: Distance of the object to the viewing position.

    Returns:
        Returns the topocentric position vector of the object in the same units that the distance parameter is in.
    """

    alt = radians(altitude)
    az = radians(azimuth)
    zComp = distance * sin(alt)

    cosAlt = distance * distance * (cos(alt) ** 2)
    xComp = sqrt(cosAlt / (1 + (tan(az) ** 2)))
    if azimuth < 90 or azimuth > 270:
        xComp = -abs(xComp)
    else:
        xComp = abs(xComp)

    yComp = -xComp * tan(az)
    if azimuth < 180:
        yComp = abs(yComp)
    else:
        yComp = -abs(yComp)

    return EVector(xComp, yComp, zComp)
