from math import atan2, pi

from sattrack.util.constants import TWOPI


def atan3(y: float, x: float) -> float:
    """A makeshift version of the atan2 method, where the return value is between 0 and 2π."""

    angle = atan2(y, x)

    if y < 0:
        return angle + TWOPI
    return angle


def signOf(value: float) -> int:
    # Zero will return True.

    if value >= 0:
        return 1
    return -1


def computeAngleDifference(angle: float) -> float:
    validAngle = angle % TWOPI
    if validAngle > pi:
        return validAngle - TWOPI
    return validAngle


def hasSameSign(a: (float, int), b: (float, int)) -> bool:
    # Will return NotImplemented if either a or b is zero.

    if a < 0:
        if b < 0:
            return True
        elif b > 0:
            return False
    elif a > 0:
        if b > 0:
            return True
        elif b < 0:
            return False
    return NotImplemented
