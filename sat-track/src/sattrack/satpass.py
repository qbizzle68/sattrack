from math import cos, radians, pi, sqrt, acos, sin, degrees, asin

from pyevspace import EVector, cross, dot, norm, vang

from sattrack.position import computeTrueAnomaly, nearestTrueAnomaly
from sattrack.rotation.order import Axis
from sattrack.rotation.rotation import getMatrix
from sattrack.structures.satellite import Satellite
from sattrack.topos import getPVector, getToposPosition
from sattrack.util.anomalies import trueToMean
from sattrack.util.conversions import atan2
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.coordinates import GeoPosition, zenithVector, geoPositionVector
from sattrack.structures.elements import computeEccentricVector, raanProcession


class PositionInfo:
    """Helper class containing the information about a given position during a satellite pass. This
    structure holds the altitude and azimuth data for an instance of a pass, for example at the rise
    point or at the point of maximum altitude.
    Parameters:
    altitude:   Altitude of the satellite above the horizon.
    azimuth:    Azimuth of the satellite measured clockwise from north.
    time:       Time of this instance in the pass.
    visible:    Visibility of the satellite, defaults to false."""

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, visible: bool = False):
        self._altitude = altitude
        self._azimuth = azimuth
        self._time = time
        self._visible = visible

    def __str__(self):
        """Generates a string of the data in this class."""
        return f'Altitude: {self._altitude}\nAzimuth: {self._azimuth}\nTime: {self._time}\nVisibility: {self._visible}'

    def altitude(self) -> float:
        """Returns the altitude measured in degrees.."""
        return self._altitude

    def azimuth(self) -> float:
        """Returns the azimuth, measured clockwise from north in degrees."""
        return self._azimuth

    def time(self) -> JulianDate:
        """Returns a JulianDate representing the time of this instance."""
        return self._time

    def visibility(self) -> bool:
        """Returns the visibility of the satellite."""
        return self._visible


class Pass:
    """set first if rise isn't visible, this determines pass visibility."""

    def __init__(self, riseInfo: PositionInfo, setInfo: PositionInfo, maxInfo: PositionInfo,
                 firstInfo: PositionInfo = None, lastInfo: PositionInfo = None):
        self._riseInfo = riseInfo
        self._setInfo = setInfo
        self._maxInfo = maxInfo
        self._firstVisible = firstInfo
        self._lastVisible = lastInfo
        if firstInfo is not None or riseInfo.visibility():
            self._visible = True
        else:
            self._visible = False

    def __str__(self):
        rtn = f'Rise:\n{self._riseInfo}\nMax:\n{self._maxInfo}\nSet:\n{self._setInfo}'
        if self._firstVisible is not None:
            rtn += f'First Visible:\n{self._firstVisible}'
        if self._lastVisible is not None:
            rtn += f'Last Visible:\n{self._lastVisible}'
        return rtn

    def riseInfo(self) -> PositionInfo:
        return self._riseInfo

    def setInfo(self) -> PositionInfo:
        return self._setInfo

    def maxInfo(self) -> PositionInfo:
        return self._maxInfo

    def firstVisibleInfo(self) -> PositionInfo:
        return self._firstVisible

    def lastVisibleInfo(self) -> PositionInfo:
        return self._lastVisible


def nextPass(sat: Satellite, geo: GeoPosition, time: JulianDate) -> Pass:
    nextPassTime = nextPassMax(sat, geo, time)
    riseTime, setTime = riseSetTimes(sat, geo, nextPassTime)

    risePos = getToposPosition(sat, riseTime, geo)
    riseInfo = PositionInfo(degrees(asin(risePos[2] / risePos.mag())),
                            degrees(atan2(risePos[1], -risePos[0])),
                            riseTime)  # todo: check visibility
    setPos = getToposPosition(sat, setTime, geo)
    setInfo = PositionInfo(degrees(asin(setPos[2] / setPos.mag())),
                           degrees(atan2(setPos[1], -setPos[0])),
                           setTime)  # todo: check visibility
    maxPos = getToposPosition(sat, nextPassTime, geo)
    maxInfo = PositionInfo(degrees(asin(maxPos[2] / maxPos.mag())),
                           degrees(atan2(maxPos[1], -maxPos[0])),
                           nextPassTime)  # todo: check visibility
    return Pass(riseInfo, setInfo, maxInfo)


def nextPassMax(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """Computes the time to the next maximum height that the next satellite pass achieves.
    Parameters:
    sat:    The satellite.
    geo:    GeoPosition the pass is observed from.
    time:   A time between passes. The pass computed will be the soonest pass after time."""
    if orbitAltitude(sat, geo, time) < 0:
        t0 = timeToPlane(sat, geo, time)
    else:
        t0 = time

    # rough estimate to next maximum height moving forward
    state = sat.getState(t0)
    pVec = getPVector(geo, state, t0)
    ma0 = trueToMean(computeTrueAnomaly(state[0], state[1]), sat.tle().eccentricity())
    eccVec = computeEccentricVector(state[0], state[1])
    ta1 = vang(eccVec, pVec)
    if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
        ta1 = 360 - ta1
    ma1 = trueToMean(ta1, sat.tle().eccentricity())
    if ma1 < ma0:
        dma0 = radians(ma1 + 360 - ma0)
    else:
        dma0 = radians(ma1 - ma0)
    tn = t0.future(dma0 / (sat.tle().meanMotion() * 2 * pi))

    # iterate towards answer moving forward or backward
    state = sat.getState(tn)
    pVec = getPVector(geo, state, tn)
    while vang(state[0], pVec) > 2.78e-4:
        pVec = getPVector(geo, state, tn)
        eccVec = computeEccentricVector(state[0], state[1])
        tan = computeTrueAnomaly(state[0], state[1])
        man = trueToMean(tan, sat.tle().eccentricity())
        tan1 = vang(eccVec, pVec)
        if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
            tan1 = 360 - tan1
        man1 = trueToMean(tan1, sat.tle().eccentricity())
        if tan1 <= tan:
            if (tan - tan1) < 180:
                dma = radians(man1 - man)
            else:
                dma = radians(man1 + 360 - man)
        else:
            if (tan1 - tan) > 180:
                dma = radians(man1 - 360 - man)
            else:
                dma = radians(man1 - man)
        tn = tn.future(dma / (sat.tle().meanMotion() * 2 * pi))
        state = sat.getState(tn)
        pVec = getPVector(geo, state, tn)
    if orbitAltitude(sat, geo, tn) < 0:
        return nextPassMax(sat, geo, tn.future(0.001))
    return tn


def riseSetTimes(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
    # time needs to be a during the pass
    a = sat.tle().sma()
    c = a * sat.tle().eccentricity()
    b = sqrt(a * a - c * c)

    state = sat.getState(time)
    hNorm = norm(cross(state[0], state[1]))
    u = norm(computeEccentricVector(state[0], state[1])) * a
    v = norm(cross(hNorm, u)) * b
    ce = -norm(u) * c

    zeta = norm(zenithVector(geo, time))
    gamma = geoPositionVector(geo, time)

    R = sqrt((dot(zeta, u) ** 2) + (dot(zeta, v) ** 2))
    try:
        beta = acos(dot(zeta, gamma - ce) / R)
    except ValueError:
        # create our own exception type here
        raise Exception("No pass during this time.")
    # todo: is this atan2 good enough? (range)
    alpha = atan2(dot(zeta, v), dot(zeta, u))
    w1 = alpha + beta
    w2 = alpha - beta

    rho1 = (u * cos(w1) + v * sin(w1)).mag()
    ta1 = degrees(atan2(rho1 * sin(w1), rho1 * cos(w1) - a * sat.tle().eccentricity()))
    rho2 = (u * cos(w2) + v * sin(w2)).mag()
    ta2 = degrees(atan2(rho2 * sin(w2), rho2 * cos(w2) - a * sat.tle().eccentricity()))
    jd1 = nearestTrueAnomaly(sat, time, ta1)
    jd2 = nearestTrueAnomaly(sat, time, ta2)
    if jd1.value() < jd2.value():
        return jd1, jd2
    else:
        return jd2, jd1


def timeToPlane(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    jd = time
    state = sat.getState(jd)
    alt = orbitAltitude(sat, geo, jd)
    if alt > 0:
        return jd
    dRaan = raanProcession(sat.tle())
    while abs(alt) > 2.7e-4:  # one arc-second on either side
        # todo: improve this guess of dt
        dt = -alt / 360.0
        jd = jd.future(dt)
        mat = getMatrix(Axis.Z_AXIS, dRaan * dt)
        state = (mat @ state[0], mat @ state[1])
        alt = orbitAltitude(sat, geo, jd)
    return jd


def orbitAltitude(sat: Satellite, geo: GeoPosition, jd: JulianDate) -> float:
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

    state = sat.getState(jd)
    ecc = sat.tle().eccentricity()
    sma = sat.tle().sma()

    # zenith vector for the GeoPosition
    zeta = norm(zenithVector(geo, jd))
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
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    p = v * t + r

    # find the true anomaly of this vector if it were a position vector
    trueAnom = vang(computeEccentricVector(state[0], state[1]), p)
    pSat = norm(p) * ((sma * (1 - ecc * ecc)) / (1 + ecc * cos(radians(trueAnom))))

    ang = vang(p - gamma, pSat - gamma)
    return ang if pSat.mag2() > p.mag2() else -ang
