from math import sin, cos

from pyevspace import EVector

from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import AU


def getSunPosition(time: JulianDate) -> EVector:
    """Computes the position of the Sun in an earth-centered reference frame,
    from the algorithm found on Wikipedia.
    Parameters:
    time:   Time for finding the Sun's position.
    returns: An earth centered position vector of the Sun in km."""

    #   time since noon TT on Jan 1, 2000
    n = time.difference() - J2000
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
