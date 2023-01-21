from math import atan2

from sattrack.util.constants import TWOPI

__all__ = ['atan3']


def atan3(y: float, x: float) -> float:
    """A makeshift version of the atan2 method, where the return value is between 0 and 2π."""

    angle = atan2(y, x)

    if y < 0:
        return angle + TWOPI
    return angle
