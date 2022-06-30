from math import cos, radians, pi

from pyevspace import EVector, cross, dot, norm, vang

from sattrack.position import computeTrueAnomaly
from sattrack.rotation.order import Axis
from sattrack.rotation.rotation import getMatrix
from sattrack.structures.satellite import Satellite
from sattrack.topos import getPVector
from sattrack.util.anomalies import trueToMean
from sattrack.util.conversions import meanMotionToSma
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.coordinates import GeoPosition, zenithVector, geoPositionVector
from sattrack.structures.elements import computeEccentricVector, raanProcession


class PositionInfo:

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, visible: bool = False):
        self._altitude = altitude
        self._azimuth = azimuth
        self._time = time
        self._visible = visible

    def __str__(self):
        return f'Altitude: {self._altitude}\nAzimuth: {self._azimuth}\nTime: {self._time}\nVisibility: {self._visible}'

    def altitude(self) -> float:
        return self._altitude

    def azimuth(self) -> float:
        return self._azimuth

    def time(self) -> JulianDate:
        return self._time

    def visibility(self) -> bool:
        return self._visible


class Pass:

    """set first if rise isn't visible, this determines pass visibility."""
    def __init__(self, rise: PositionInfo, set: PositionInfo, max: PositionInfo, first: PositionInfo = None, last: PositionInfo = None):
        self._riseInfo = rise
        self._setInfo = set
        self._maxInfo = max
        self._firstVisible = first
        self._lastVisible = last
        if first is not None or rise.visibility():
            self._visible = True
        else:
            self._visible = False

    def __str__(self):
        rtn = f'Rise:\n{self._riseInfo}\nMax:\n{self._maxInfo}\nSet:\n{self._setInfo}'
        if self._firstVisible is not None:
            rtn += f'First Visible:\n{self._firstVisible}'
        if self._lastVisible is not None:
            rtn += f'Last Visible:\n{self._lastVisible}'

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


def nextPass(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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
        return nextPass(sat, geo, tn.future(0.001))
    return tn


def timeToPlane(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    jd = time
    state = sat.getState(jd)
    ecc = sat.tle().eccentricity()
    sma = meanMotionToSma(sat.tle().meanMotion())
    #alt = orbitAltitude(geo, jd, state, ecc, sma)
    alt = orbitAltitude(sat, geo, jd)
    if alt > 0:
        return jd
    dRaan = raanProcession(sat.tle())
    while abs(alt) > 2.7e-4: # one arc-second on either side
        # todo: improve this guess of dt
        dt = -alt / 360.0
        jd = jd.future(dt)
        mat = getMatrix(Axis.Z_AXIS, dRaan * dt)
        state = (mat @ state[0], mat @ state[1])
        #alt = orbitAltitude(geo, jd, state, ecc, sma)
        alt = orbitAltitude(sat, geo, jd)
    return jd


#def orbitAltitude(geo: GeoPosition, jd: JulianDate, state: tuple[EVector], ecc: float, sma: float) -> float:
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
    # t = dot(v, gamma - r) / dot(v, v)
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    p = v * t + r

    # find the true anomaly of this vector if it were a position vector
    trueAnom = vang(computeEccentricVector(state[0], state[1]), p)
    pSat = norm(p) * ((sma * (1 - ecc * ecc)) / (1 + ecc * cos(radians(trueAnom))))

    ang = vang(p - gamma, pSat - gamma)
    return ang if pSat.mag2() > p.mag2() else -ang
