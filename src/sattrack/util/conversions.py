from math import atan2

from sattrack.util.constants import TWOPI


def atan3(y: float, x: float) -> float:
    # atan2 but range is 0 to 2π
    angle = atan2(y, x)

    if y < 0:
        return angle + TWOPI
    return angle
