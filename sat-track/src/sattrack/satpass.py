from math import cos, radians, pi, sqrt, acos, sin, degrees, asin

from pyevspace import EVector, cross, dot, norm, vang

from sattrack.position import computeTrueAnomaly, nearestTrueAnomaly
from sattrack.rotation.order import Axis, Order
from sattrack.rotation.rotation import getMatrix, rotateOrderTo, EulerAngles
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.satellite import Satellite
from sattrack.sun import getSunPosition, TwilightType
from sattrack.topos import getPVector, toTopocentric, getTwilightType, getAltitude
from sattrack.util.anomalies import trueToMean
from sattrack.util.constants import EARTH_EQUITORIAL_RADIUS, SUN_RADIUS
from sattrack.util.conversions import atan2
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.coordinates import GeoPosition, zenithVector, geoPositionVector
from sattrack.structures.elements import computeEccentricVector, raanProcessionRate


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
        return f'Altitude: {"%.2f" % self._altitude}\nAzimuth: {"%.2f" % self._azimuth}\nTime: {self._time}' \
               f'\nVisibility: {self._visible}'

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


class PassConstraints:
    """Container class to hold values to constrain satellite passes.

    Attributes:
        minAltitude -- minimum altitude the satellite must achieve in degrees
        minDuration -- minimum duration the satellite must be above the horizon in minutes
        illuminated -- whether the satellite is illuminated at any point of the pass
    """

    def __init__(self, *, minAltitude: float = None, minDuration: float = None, illuminated: bool = None):
        self.minAltitude = minAltitude
        self.minDuration = minDuration
        self.illuminated = illuminated


def nextPass(sat: Satellite, geo: GeoPosition, time: JulianDate,
             constraints: PassConstraints = None) -> Pass:
    nextPassTime = nextPassMax(sat, geo, time)
    maxPos = sat.getState(nextPassTime)[0]
    maxPosSez = toTopocentric(maxPos, nextPassTime, geo)
    maxAlt = degrees(asin(maxPosSez[2] / maxPosSez.mag()))

    if constraints is not None and constraints.minAltitude is not None:
        if constraints.minAltitude > maxAlt:
            return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    riseTime, setTime = riseSetTimes(sat, geo, nextPassTime)
    risePos = sat.getState(riseTime)[0]
    risePosSez = toTopocentric(risePos, riseTime, geo)
    riseIlluminated = not isEclipsed(risePos, getSunPosition(riseTime))
    setPos = sat.getState(setTime)[0]
    setPosSez = toTopocentric(setPos, setTime, geo)
    setIlluminated = not isEclipsed(setPos, getSunPosition(setTime))

    if constraints is not None:
        if constraints.illuminated is not None:
            if constraints.illuminated and not riseIlluminated and not setIlluminated:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
            if not constraints.illuminated and riseIlluminated and setIlluminated:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
        if constraints.minDuration is not None:
            if constraints.minDuration > setTime.difference(riseTime) * 1440:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    riseInfo = PositionInfo(degrees(asin(risePosSez[2] / risePosSez.mag())),
                            degrees(atan2(risePosSez[1], -risePosSez[0])),
                            riseTime,
                            riseIlluminated and getTwilightType(riseTime, geo) > TwilightType.Day)

    setInfo = PositionInfo(degrees(asin(setPosSez[2] / setPosSez.mag())),
                           degrees(atan2(setPosSez[1], -setPosSez[0])),
                           setTime,
                           setIlluminated and getTwilightType(setTime, geo) > TwilightType.Day)

    maxIlluminated = not isEclipsed(maxPos, getSunPosition(nextPassTime))
    maxInfo = PositionInfo(degrees(asin(maxPosSez[2] / maxPosSez.mag())),
                           degrees(atan2(maxPosSez[1], -maxPosSez[0])),
                           nextPassTime,
                           maxIlluminated and getTwilightType(nextPassTime, geo) > TwilightType.Day)

    return Pass(riseInfo, setInfo, maxInfo)


def getPassList(sat: Satellite, geo: GeoPosition, start: JulianDate, duration: float,
                constraints: PassConstraints = None) -> tuple[Pass]:
    passList = []
    nTime = start.future(-0.001)
    while nTime.difference(start) < duration:
        nPass = nextPass(sat, geo, nTime.future(0.001), constraints)
        nTime = nPass.maxInfo().time()
        if nPass.maxInfo().time().difference(start) < duration:
            passList.append(nPass)
    return tuple(passList)


def nextPassMax(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    nextMax = nextPassMaxGuess(sat, geo, time)
    return maxPassRefine(sat, geo, nextMax)


def nextPassMaxGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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
    ma0 = trueToMean(computeTrueAnomaly(state[0], state[1]), sat.tle().getEcc())
    eccVec = computeEccentricVector(state[0], state[1])
    ta1 = radians(vang(eccVec, pVec))
    if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
        ta1 = (2*pi) - ta1
    ma1 = trueToMean(ta1, sat.tle().getEcc())
    if ma1 < ma0:
        dma0 = ma1 + (2*pi) - ma0
        # dma0 = radians(ma1 + 360 - ma0)
    else:
        dma0 = ma1 - ma0
        # dma0 = radians(ma1 - ma0)
    tn = t0.future(dma0 / (sat.tle().getMeanMotion() * 2 * pi))

    # iterate towards answer moving forward or backward
    state = sat.getState(tn)
    pVec = getPVector(geo, state, tn)
    while vang(state[0], pVec) > 2.78e-4:
        pVec = getPVector(geo, state, tn)
        eccVec = computeEccentricVector(state[0], state[1])
        tan = computeTrueAnomaly(state[0], state[1])
        man = trueToMean(tan, sat.tle().getEcc())
        tan1 = radians(vang(eccVec, pVec))
        if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
            tan1 = (2*pi) - tan1
            # tan1 = 360 - tan1
        man1 = trueToMean(tan1, sat.tle().getEcc())
        if tan1 <= tan:
            if (tan - tan1) < pi:
                dma = man1 - man
                # dma = radians(man1 - man)
            else:
                dma = man1 + (2*pi) - man
                # dma = radians(man1 + 360 - man)
        else:
            if (tan1 - tan) > pi:
                dma = man1 - (2*pi) - man
                # dma = radians(man1 - 360 - man)
            else:
                dma = man1 - man
                # dma = radians(man1 - man)
        tn = tn.future(dma / (sat.tle().getMeanMotion() * 2 * pi))
        state = sat.getState(tn)
        pVec = getPVector(geo, state, tn)
    if orbitAltitude(sat, geo, tn) < 0:
        return nextPassMax(sat, geo, tn.future(0.001))
    return tn


def maxPassRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    # todo: utilize the dadt values to future time increase
    '''this works, but feels like overkill (many iterations with small dt) to achieve accurate answer.
    can we use either the dadt or break it down into dzdt and dxydt terms, and create a meaningful empirical guess?'''

    alt = getAltitude(sat, time, geo)
    futureAlt = getAltitude(sat, time.future(0.1 / 86400), geo)
    pastAlt = getAltitude(sat, time.future(-0.1 / 86400), geo)

    parity = 0  # to please the editor
    if alt >= futureAlt and alt >= pastAlt:
        return time
    elif alt < futureAlt:
        parity = 1
    elif alt < pastAlt:
        parity = -1

    jd = time
    nextAlt = getAltitude(sat, time.future(parity * 0.1 / 86400), geo)
    while alt < nextAlt:
        jd = jd.future(parity * 0.1 / 86400)
        alt = getAltitude(sat, jd, geo)
        nextAlt = getAltitude(sat, jd.future(parity * 0.1 / 86400), geo)

    return jd


def riseSetTimes(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
    # time needs to be a during the pass
    riseTime, setTime = riseSetGuess(sat, geo, time)
    riseTime = horizonTimeRefine(sat, geo, riseTime)
    setTime = horizonTimeRefine(sat, geo, setTime)
    return riseTime, setTime


def riseSetGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
    # time needs to be a during the pass
    a = sat.tle().getSma()
    c = a * sat.tle().getEcc()
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
        # todo: create our own exception type here
        raise Exception("No pass during this time.")
    alpha = atan2(dot(zeta, v), dot(zeta, u))
    w1 = alpha + beta
    w2 = alpha - beta

    rho1 = (u * cos(w1) + v * sin(w1)).mag()
    ta1 = atan2(rho1 * sin(w1), rho1 * cos(w1) - a * sat.tle().getEcc())
    rho2 = (u * cos(w2) + v * sin(w2)).mag()
    ta2 = atan2(rho2 * sin(w2), rho2 * cos(w2) - a * sat.tle().getEcc())
    jd1 = nearestTrueAnomaly(sat, time, ta1)
    jd2 = nearestTrueAnomaly(sat, time, ta2)
    if jd1.value() < jd2.value():
        return jd1, jd2
    else:
        return jd2, jd1


def horizonTimeRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    state = sat.getState(time)
    sezPos = toTopocentric(state[0], time, geo)
    alt = asin(sezPos[2] / sezPos.mag())
    while abs(alt) > radians(1 / 3600):
        sezVel = rotateOrderTo(
            Order.ZYX,
            EulerAngles(
                radians(geo.getLongitude() + earthOffsetAngle(time)),
                radians(90 - geo.getLatitude()),
                0.0
            ),
            state[1]
        )

        dz = sezPos.mag() * alt
        dt = (-dz / sezVel[2]) / 86400.0
        time = time.future(dt)

        state = sat.getState(time)
        sezPos = toTopocentric(state[0], time, geo)
        alt = asin(sezPos[2] / sezPos.mag())
    return time


def timeToPlane(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    jd = time
    state = sat.getState(jd)
    alt = orbitAltitude(sat, geo, jd)
    if alt > 0:
        return jd
    dRaan = raanProcessionRate(sat.tle())
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
    ecc = sat.tle().getEcc()
    sma = sat.tle().getSma()

    # zenith vector for the GeoPosition
    zeta = norm(zenithVector(geo, jd))
    # GeoPosition vector in geocentric reference frame
    gamma = geoPositionVector(geo, jd)
    # normalized angular momentum, vector equation for orbital plane
    lamb = norm(cross(state[0], state[1]))

    # todo: ensure lamb[1] != 0 and use other solution if it is
    # compute intermediate values to find solution to parameterized vector intersection
    x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
    y = -lamb[0] * x / lamb[1]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)

    # compute exact solution for parametrized vector which yields the nearest point to intersection
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    p = v * t + r

    # find the true anomaly of this vector if it were a position vector
    trueAnom = radians(vang(computeEccentricVector(state[0], state[1]), p))
    pSat = norm(p) * ((sma * (1 - ecc * ecc)) / (1 + ecc * cos(trueAnom)))

    ang = vang(p - gamma, pSat - gamma)
    return ang if pSat.mag2() > p.mag2() else -ang


def isEclipsed(satPos: EVector, sunPos: EVector) -> bool:
    """Determines whether a position is eclipsed by the earth.
    Parameters:
    satPos: Position vector of object.
    sunPos: Position of the Sun at the same time the object is at satPos.
    returns True if eclipsed by any means, false otherwise."""

    #   vectors are relative to the satellite
    earthPos = -satPos
    sunPos = -satPos + sunPos
    #   semi-diameters of earth and sun
    # todo: calculate real Earth radius from perspective here
    thetaE = asin(6371 / earthPos.mag()) # average earth radius
    thetaS = asin(SUN_RADIUS / sunPos.mag())
    #   angle between earth and sun centers relative to the satellite
    theta = radians(vang(earthPos, sunPos))

    #   umbral eclipse
    if thetaE > thetaS and theta < (thetaE - thetaS):
        return True
    #   penumbral eclipse
    if abs(thetaE - thetaS) < theta < (thetaE + thetaS):
        return True
    return thetaS > thetaE and theta < (thetaS - thetaE)
