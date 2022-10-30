from math import atan2, sqrt

from sattrack.structures.body import Body, EARTH_BODY
from sattrack.util.constants import TWOPI


def atan3(y: float, x: float) -> float:
    # atan2 but range is 0 to 2π
    angle = atan2(y, x)

    if y < 0:
        return angle + TWOPI
    return angle


# def meanMotionToSma(meanMotion: float, body: Body = EARTH_BODY) -> float:
#     """
#     Converts a satellites mean motion from a two-line element set to the orbit's semi-major axis using Kepler's 2nd law.
#
#     Args:
#         meanMotion: Mean motion of the satellite measured in revolutions per day.
#         body: Body of the orbiting satellite (Default = EARTH_BODY).
#
#
#     Returns:
#         The semi-major axis in kilometers.
#     """
#
#     mMotionRad = meanMotion * TWOPI / 86400.0
#     return (body.getMu() ** (1.0 / 3.0)) / (mMotionRad ** (2.0 / 3.0))
#
#
# def smaToMeanMotion(sma: float, body: Body = EARTH_BODY) -> float:
#     """
#     Converts a satellite's semi-major axis to mean motion.
#
#     Args:
#         sma: Semi-major axis measured in kilometers.
#         body: Body of the orbiting satellite (Default = EARTH_BODY).
#
#     Returns:
#         The mean motion in radians per second.
#     """
#
#     return sqrt(body.getMu() / (sma ** 3))
