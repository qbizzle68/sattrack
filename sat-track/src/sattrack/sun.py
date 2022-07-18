from enum import Enum
from functools import total_ordering
from math import sin, cos

from pyevspace import EVector

from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import AU


@total_ordering
class TwilightType(Enum):
    """
    Enumerated value to distinguish between the various sun angles relative to a horizon.
    The classifications for sun angle to twilight type is:
    Sun < -18 degrees --> Night time
    Sun < -12 degrees --> Astronomical twilight
    Sun < -6 degrees --> Nautical twilight
    Sun < -50 arc-minutes --> Civil twilight
    Otherwise --> Day time
    Most object can't be seen until their position is in civil twilight.
    """

    Day = 0
    Civil = 1
    Nautical = 2
    Astronomical = 3
    Night = 4

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented


def getSunPosition(time: JulianDate) -> EVector:
    """
    Very simple algorithm for computing the position of the Sun in a geocentric equitorial reference frame.

    Args:
        time: Time to find the Sun's position.

    Returns:
        The Sun's position vector in kilometers.
    """

    #   time since noon TT on Jan 1, 2000
    n = time.difference(J2000)
    #   mean longitude of the Sun
    L = 4.89495042 + 0.0172027924 * n
    #   mean anomaly of the Sun
    g = 6.240040768 + 0.0172019703 * n
    #   ecliptic longitude of the Sun
    lmbda = L + 0.0334230552 * sin(g) + 0.00034906585 * sin(2 * g)
    #   distance of the Sun from the Earth
    R = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2 * g)
    #   obliquity of the ecliptic (axial tilt)
    e = 0.4090877234 - 6.981317008e-9 * n

    return EVector(
        R * AU * cos(lmbda),
        R * AU * cos(e) * sin(lmbda),
        R * AU * sin(e) * sin(lmbda)
    )
