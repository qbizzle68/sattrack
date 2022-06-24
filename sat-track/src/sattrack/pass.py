from math import cos, radians

from pyevspace import EVector, cross, dot, norm, vang

from structures.coordinates import GeoPosition, zenithVector, geoPositionVector
from structures.elements import OrbitalElements, computeEccentricVector
from spacetime.juliandate import JulianDate


def orbitAltitude(geo: GeoPosition, jd: JulianDate, state: tuple[EVector], ecc: float, sma: float) -> float:
    """Computes the angle above or below the horizon, of the nearest point along the orbit to
    a GeoPosition at a given time. This in essence tells you if the path of the orbit can be seen
    above a given horizon. A negative value indicates the nearest point on the orbital path
    is below the horizon, where a positive path indicates both that an overhead pass is
    possible, and tells you the maximum height a pass can achieve.
    Parameters:
    geo:    GeoPosition of the horizon.
    jd:     Time to of the horizon.
    state:  State vectors of the satellite.
    ecc:    Eccentricity of the orbit.
    sma:    Semi-major axis of the orbit in meters."""

    # zenith vector for the GeoPosition
    zeta = zenithVector(geo, jd)
    # GeoPosition vector in geocentric reference frame
    gamma = geoPositionVector(geo, jd)
    # normalized angular momentum, vector equation for orbital plane
    lamb = norm(cross(state[0], state[1]))

    # compute intermediate values to find solution to parameterized vector intersection
    x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
    y = -lamb[0] * x / lamb[1]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)

    # compute exact solution for parametrized vector which yields the nearest point to intersection
    t = dot(v, gamma - r) / dot(v, v)
    p = v * t + r

    # find the true anomaly of this vector if it were a position vector
    trueAnom = vang(computeEccentricVector(state[0], state[0]), p)
    pSat = norm(p) * ((sma * (1 - ecc * ecc)) / (1 + ecc * cos(radians(trueAnom))))

    ang = vang(p - gamma, pSat - gamma)
    return ang if pSat.mag2() < p.mag2() else -ang
