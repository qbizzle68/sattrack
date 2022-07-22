from math import cos, radians, pi, sqrt, acos, sin, degrees, asin, atan

from pyevspace import EVector, cross, dot, norm, vang

from sattrack.exceptions import NoPassException
from sattrack.rotation.order import Axis, ZYX, ZXZ
from sattrack.rotation.rotation import getMatrix, rotateOrderTo, EulerAngles
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures.coordinates import GeoPosition, zenithVector, geoPositionVector
from sattrack.structures.elements import computeEccentricVector, raanProcessionRate, OrbitalElements
from sattrack.structures.satellite import Satellite
from sattrack.structures.tle import TwoLineElement
from sattrack.sun import getSunPosition, TwilightType, getSunPosition2
from sattrack.topos import getPVector, toTopocentric, getTwilightType, getAltitude, azimuthAngleString
from sattrack.util.anomalies import trueToMean, timeToNearestTrueAnomaly, computeTrueAnomaly
from sattrack.util.constants import SUN_RADIUS, TWOPI, EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS
from sattrack.util.conversions import atan3


class PositionInfo:
    """
    A helper class used for storing the information about a given position during a satellite pass. This structure holds
    the altitude and azimuth data for an instance of a pass, for example at the rise point or at the point of maximum
    altitude.

    Attributes
        altitude: Altitude of the satellite above the horizon.
        azimuth: Azimuth of the satellite measured clockwise from north.
        direction: Compass direction of the azimuth direction.
        time: Time of this instance in the pass.
        visible: Visibility of the satellite, defaults to false.
    """

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, visible: bool = False):
        """
        Initializes the values to the parameters.
        Args:
            altitude: Altitude of the satellite in degrees.
            azimuth: Azimuth of the satellite in degrees.
            time: Time of this instance in the pass.
            visible: Visibility of the satellite (Default = False).
        """

        self._altitude = altitude
        self._azimuth = azimuth
        self._direction = azimuthAngleString(azimuth)
        self._time = time
        self._visible = visible

    def __str__(self):
        """Generates a string of the data in this class."""
        return f'Altitude: {"%.2f" % self._altitude}\nAzimuth: {"%.2f" % self._azimuth} ({self._direction})\nTime: ' \
               f'{self._time}\nVisibility: {self._visible}'

    def getAltitude(self) -> float:
        """Returns the altitude measured in degrees.."""
        return self._altitude

    def getAzimuth(self) -> float:
        """Returns the azimuth, measured clockwise from north in degrees."""
        return self._azimuth

    def getDirection(self) -> str:
        """Returns the azimuth compass direction."""
        return self._direction

    def getTime(self) -> JulianDate:
        """Returns a JulianDate representing the time of this instance."""
        return self._time

    def getVisibility(self) -> bool:
        """Returns the visibility of the satellite."""
        return self._visible


class Pass:
    """
    A class which contains the information for an overhead satellite pass. The class contains information for the rise
    and set times as well as the time it achieves maximum altitude. If the satellite enters and/or exits the Earths
    shadow during the pass, information regarding the time and position of their occurrences can be held by the Pass
    class as well.

    Attributes
        riseInfo: PositionInfo for the rise time of the satellite pass.
        setInfo: PositionInfo for the set time of the satellite pass.
        maxInfo: PositionInfo for the max time of the satellite pass.
        firstInfo: PositionInfo for when the satellite is illuminated, if not the same as riseInfo.
        lastInfo: PositionInfo for when the satellite is eclipsed, if not the same as setInfo.
        visible: Boolean value for if the pass is visible at any point during the transit.
    """

    def __init__(self, riseInfo: PositionInfo, setInfo: PositionInfo, maxInfo: PositionInfo, *,
                 firstInfo: PositionInfo = None, lastInfo: PositionInfo = None):
        """
        Initializes a pass with the given values. firstInfo and lastInfo should only be set if the satellite enters or
        leaves the Earth's shadow during the transit. If the riseInfo visibility attribute is true or a firstInfo value
        is set, the visibility of the pass is set to true, otherwise it is set to false.

        Args:
            riseInfo: Information for the rise time of the satellite pass.
            setInfo: Information for the set time of the satellite pass.
            maxInfo: Information for the maximum time of the satellite pass.
            firstInfo: Information for the first time the satellite is illuminated if it's not the same as the rise time
                (Default = None).
            lastInfo: Information for the last time the satellite is illuminated if it's not the same as the set time
                (Default = None).
        """

        self._riseInfo = riseInfo
        self._setInfo = setInfo
        self._maxInfo = maxInfo
        self._firstVisible = firstInfo
        self._lastVisible = lastInfo
        if firstInfo is not None or riseInfo.getVisibility():
            self._visible = True
        else:
            self._visible = False

    def __str__(self):
        """Returns a string containing all the pass's information."""
        riseStr = f'Rise:\n{self._riseInfo}'
        setStr = f'Set:\n{self._setInfo}'
        maxStr = f'Max:\n{self._maxInfo}'
        strDict = {info.getTime().value(): string for info, string in
                   zip((self._riseInfo, self._setInfo, self._maxInfo), (riseStr, setStr, maxStr))}
        if self._firstVisible is not None:
            firstStr = f'First:\n{self._firstVisible}'
            strDict[self._firstVisible.getTime().value()] = firstStr
        if self._lastVisible is not None:
            lastStr = f'Last:\n{self._lastVisible}'
            strDict[self._lastVisible.getTime().value()] = lastStr
        rtn = ""
        for i in range(len(strDict)):
            minTm = min(strDict.keys())
            rtn += strDict[minTm] + '\n'
            strDict.pop(minTm)
        return rtn.rstrip()

    def getRiseInfo(self) -> PositionInfo:
        """Returns the rise time information."""
        return self._riseInfo

    def getSetInfo(self) -> PositionInfo:
        """Returns the set time information."""
        return self._setInfo

    def getMaxInfo(self) -> PositionInfo:
        """Returns the max time information."""
        return self._maxInfo

    def getFirstVisibleInfo(self) -> PositionInfo:
        """Returns the first visible time information if set."""
        return self._firstVisible

    def getLastVisibleInfo(self) -> PositionInfo:
        """Returns the last visible time information if set."""
        return self._lastVisible

    def getVisibility(self) -> bool:
        """Returns the visibility of the pass."""
        return self._visible


class PassConstraints:
    """
    Container class to hold values to constrain satellite passes.

    Attributes:
        minAltitude: Minimum altitude the satellite must achieve in degrees.
        minDuration: Minimum duration the satellite must be above the horizon in minutes.
        illuminated: Whether the satellite is illuminated at any point of the pass.
    """

    def __init__(self, *, minAltitude: float = None, minDuration: float = None, illuminated: bool = None):
        """Initializes the values to the argument values."""
        self.minAltitude = minAltitude
        self.minDuration = minDuration
        self.illuminated = illuminated


# todo: check the ordering of this method, see if checks can be made sooner without in order to gain efficiency
def nextPass(sat: Satellite, geo: GeoPosition, time: JulianDate,
             constraints: PassConstraints = None) -> Pass:
    """
    Computes the next overhead pass of a satellite for a geo-position. A PassConstraints object can be set to restrict
    which passes one wishes to find. The default behaviour is any pass is computed, regardless of visibility or height
    (as long as it's greater than zero).

    Args:
        sat: Satellite to find the next pass of.
        geo: GeoPosition to view the pass.
        time: Relative time to find the next pass.
        constraints: A PassConstraints object to constrain allowable passes (Default = 0).

    Returns:
        A Pass object containing the details of the satellite pass.
    """

    nextPassTime = nextPassMax(sat, geo, time)
    maxPos = sat.getState(nextPassTime)[0]
    maxPosSez = toTopocentric(maxPos, nextPassTime, geo)
    maxAlt = degrees(asin(maxPosSez[2] / maxPosSez.mag()))

    # todo: issue with this, if lastTime < maxTime, then should lastTime be the max time?
    if constraints is not None and constraints.minAltitude is not None:
        #   if minimum altitude isn't attained
        if constraints.minAltitude > maxAlt:
            return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    sc = ShadowController(sat.getTle())
    sc.computeValues(nextPassTime)
    enterTime, exitTime = sc.getTimes()
    firstInfo, lastInfo = None, None

    riseTime, setTime = riseSetTimes(sat, geo, nextPassTime)
    if riseTime.value() < enterTime.value() < setTime.value():
        riseIlluminated = True
        firstTime = riseTime
        setIlluminated = False
        lastTime = enterTime
    elif riseTime.future(1 / sat.getTle().getMeanMotion()).value() < exitTime.value() and \
            riseTime.value() < exitTime.value() < setTime.value():
        riseIlluminated = False
        firstTime = exitTime
        setIlluminated = True
        lastTime = setTime
    elif enterTime.value() < riseTime.value() < exitTime.value() \
            and enterTime.value() < setTime.value() < exitTime.value():
        riseIlluminated = False
        firstTime = riseTime  # avoids setting the firstInfo object
        setIlluminated = False
        lastTime = setTime  # avoids setting the lastInfo object
    else:
        riseIlluminated = True
        firstTime = riseTime
        setIlluminated = True
        lastTime = setTime

    risePos = sat.getState(riseTime)[0]
    risePosSez = toTopocentric(risePos, riseTime, geo)
    setPos = sat.getState(setTime)[0]
    setPosSez = toTopocentric(setPos, setTime, geo)

    if constraints is not None:
        if constraints.illuminated is not None:
            #   when constraints is illuminated and sat is not illuminated at rise and set times
            if constraints.illuminated and not riseIlluminated and not setIlluminated:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
            #   when constraints is not illuminated and it is illuminated at some point
            #   todo: can this be constrained using relations of enter and exit times with rise and set?
            if not constraints.illuminated and (riseIlluminated or firstTime != riseTime) \
                    or (setIlluminated or lastTime != setTime):
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
        if constraints.minDuration is not None:
            #   when minimum duration is not met
            if constraints.minDuration > setTime.difference(riseTime) * 1440:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    riseInfo = PositionInfo(degrees(asin(risePosSez[2] / risePosSez.mag())),
                            degrees(atan3(risePosSez[1], -risePosSez[0])),
                            riseTime,
                            riseIlluminated and getTwilightType(riseTime, geo) > TwilightType.Day)

    setInfo = PositionInfo(degrees(asin(setPosSez[2] / setPosSez.mag())),
                           degrees(atan3(setPosSez[1], -setPosSez[0])),
                           setTime,
                           setIlluminated and getTwilightType(setTime, geo) > TwilightType.Day)

    maxIlluminated = not isEclipsed(maxPos, getSunPosition(nextPassTime))
    maxInfo = PositionInfo(degrees(asin(maxPosSez[2] / maxPosSez.mag())),
                           degrees(atan3(maxPosSez[1], -maxPosSez[0])),
                           nextPassTime,
                           maxIlluminated and getTwilightType(nextPassTime, geo) > TwilightType.Day)

    if riseTime != firstTime:
        firstPos = sat.getState(exitTime)[0]
        firstPosSez = toTopocentric(firstPos, firstTime, geo)
        firstInfo = PositionInfo(degrees(asin(firstPosSez[2] / firstPosSez.mag())),
                                 degrees(atan3(firstPosSez[1], -firstPosSez[0])),
                                 firstTime,
                                 getTwilightType(firstTime, geo) > TwilightType.Day)

    if setTime != lastTime:
        lastPos = sat.getState(enterTime)[0]
        lastPosSez = toTopocentric(lastPos, lastTime, geo)
        lastInfo = PositionInfo(degrees(asin(lastPosSez[2] / lastPosSez.mag())),
                                degrees(atan3(lastPosSez[1], -lastPosSez[0])),
                                lastTime,
                                getTwilightType(lastTime, geo) > TwilightType.Day)

    return Pass(riseInfo, setInfo, maxInfo, firstInfo=firstInfo, lastInfo=lastInfo)


def getPassList(sat: Satellite, geo: GeoPosition, start: JulianDate, duration: float,
                constraints: PassConstraints = None) -> tuple[Pass]:
    """
    Generates a list of overhead satellite passes. A PassConstraints object can be used to restrict the types of passes
    allowed in the list.

    Args:
        sat: The satellite to find overhead passes of.
        geo: The GeoPosition which the passes are viewed from.
        start: The time to start finding valid passes.
        duration: The duration of which to find passes, in solar days.
        constraints: A PassConstraints object to constrain allowable passes (Default = None).

    Returns:
        A tuple of Pass objects in chronological order, each containing information for unique passes.
    """

    passList = []
    nTime = start.future(-0.001)
    while nTime.difference(start) < duration:
        nPass = nextPass(sat, geo, nTime.future(0.001), constraints)
        nTime = nPass.getMaxInfo().getTime()
        if nPass.getMaxInfo().getTime().difference(start) < duration:
            passList.append(nPass)
    return tuple(passList)


def nextPassMax(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """
    Computes the maximum time of the next pass from the time parameter.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: The relative time to find the next maximum altitude.

    Returns:
        The time the satellite achieves the next maximum altitude.
    """

    nextMax = nextPassMaxGuess(sat, geo, time)
    return maxPassRefine(sat, geo, nextMax)


def nextPassMaxGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """
    Computes the initial guess the satellite achieves its maximum height.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: The relative time to find the next maximum altitude.

    Returns:
        The approximate time the satellite achieves the next maximum altitude.
    """

    if orbitAltitude(sat, geo, time) < 0:
        t0 = timeToPlane(sat, geo, time)
    else:
        t0 = time

    # rough estimate to next maximum height moving forward
    state = sat.getState(t0)
    pVec = getPVector(geo, *state, t0)
    ma0 = trueToMean(computeTrueAnomaly(state[0], state[1]), sat.getTle().getEcc())
    eccVec = computeEccentricVector(state[0], state[1])
    ta1 = radians(vang(eccVec, pVec))
    if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
        ta1 = TWOPI - ta1
    ma1 = trueToMean(ta1, sat.getTle().getEcc())
    if ma1 < ma0:
        dma0 = ma1 + TWOPI - ma0
    else:
        dma0 = ma1 - ma0
    tn = t0.future(dma0 / (sat.getTle().getMeanMotion() * TWOPI))

    # iterate towards answer moving forward or backward
    state = sat.getState(tn)
    pVec = getPVector(geo, *state, tn)
    while vang(state[0], pVec) > 2.78e-4:
        pVec = getPVector(geo, *state, tn)
        eccVec = computeEccentricVector(state[0], state[1])
        tan = computeTrueAnomaly(state[0], state[1])
        man = trueToMean(tan, sat.getTle().getEcc())
        tan1 = radians(vang(eccVec, pVec))
        if norm(cross(eccVec, pVec)) != norm(cross(state[0], state[1])):
            tan1 = TWOPI - tan1
        man1 = trueToMean(tan1, sat.getTle().getEcc())
        if tan1 <= tan:
            if (tan - tan1) < pi:
                dma = man1 - man
            else:
                dma = man1 + TWOPI - man
        else:
            if (tan1 - tan) > pi:
                dma = man1 - TWOPI - man
            else:
                dma = man1 - man
        tn = tn.future(dma / (sat.getTle().getMeanMotion() * TWOPI))
        state = sat.getState(tn)
        pVec = getPVector(geo, *state, tn)
    if orbitAltitude(sat, geo, tn) < 0:
        return nextPassMax(sat, geo, tn.future(0.001))
    return tn


def maxPassRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """
    Refines the maximum height of a satellite pass to it's exact value.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: The relative time to find the next maximum altitude.

    Returns:
        The refined time the satellite achieves the next maximum altitude.
    """
    # todo: utilize the dadt values to future time increase and make this a more empirically derived guess

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
    """
    Computes the rise and set times of a satellite pass.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: Must be a time that occurs during a pass, i.e. between the rise and set times.

    Returns:
        A tuple with the rise and set times of the pass.
    """

    riseTime, setTime = riseSetGuess(sat, geo, time)
    riseTime = horizonTimeRefine(sat, geo, riseTime)
    setTime = horizonTimeRefine(sat, geo, setTime)
    return riseTime, setTime


def riseSetGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
    """
    Estimates the rise and set times of a satellite pass.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: Must be a time that occurs during a pass, i.e. between the rise and set times.

    Returns:
        A tuple with the estimated rise and set times of the pass.
    """

    a = sat.getTle().getSma()
    c = a * sat.getTle().getEcc()
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
        raise NoPassException(f'No pass during this time: {time}.')
    alpha = atan3(dot(zeta, v), dot(zeta, u))
    w1 = alpha + beta
    w2 = alpha - beta

    rho1 = (u * cos(w1) + v * sin(w1)).mag()
    ta1 = atan3(rho1 * sin(w1), rho1 * cos(w1) - a * sat.getTle().getEcc())
    rho2 = (u * cos(w2) + v * sin(w2)).mag()
    ta2 = atan3(rho2 * sin(w2), rho2 * cos(w2) - a * sat.getTle().getEcc())
    ta10 = computeTrueAnomaly(*sat.getState(time))
    ta20 = computeTrueAnomaly(*sat.getState(time))
    n = sat.getTle().getMeanMotion() * TWOPI / 86400.0
    jd1 = timeToNearestTrueAnomaly(n, sat.getTle().getEcc(), ta10, time, ta1)
    jd2 = timeToNearestTrueAnomaly(n, sat.getTle().getEcc(), ta20, time, ta2)
    if jd1.value() < jd2.value():
        return jd1, jd2
    else:
        return jd2, jd1


def horizonTimeRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """
    Refines a rise or set time of a satellite pass to its correct values.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: The rise or set time to be refined.

    Returns:
        The corrected rise or set time.
    """

    state = sat.getState(time)
    sezPos = toTopocentric(state[0], time, geo)
    alt = asin(sezPos[2] / sezPos.mag())
    while abs(alt) > radians(1 / 3600):
        sezVel = rotateOrderTo(
            ZYX,
            EulerAngles(
                radians(geo.getLongitude()) + earthOffsetAngle(time),
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
    """
    Computes the time until the orbital path begins to rise above the horizon of a given geo-position.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: Relative time of which to find the next time the orbital path is visible.

    Returns:
        The next time the orbital path rises above the horizon.
    """
    jd = time
    state = sat.getState(jd)
    alt = orbitAltitude(sat, geo, jd)
    if alt > 0:
        return jd
    dRaan = raanProcessionRate(sat.getTle())
    while abs(alt) > 2.7e-4:  # one arc-second on either side
        # todo: improve this guess of dt
        dt = -alt / 360.0
        jd = jd.future(dt)
        mat = getMatrix(Axis.Z_AXIS, dRaan * dt)
        state = (mat @ state[0], mat @ state[1])
        alt = orbitAltitude(sat, geo, jd)
    return jd


def orbitAltitude(sat: Satellite, geo: GeoPosition, time: JulianDate) -> float:
    """
    Computes the angle above or below the horizon, of the nearest point along the orbit to a GeoPosition at a given
    time. This in essence tells you if the path of the orbit can be seen above a given horizon. A negative value
    indicates the nearest point on the orbital path is below the horizon, where a positive value indicates both that an
    overhead pass is possible, and tells you the maximum height a pass can achieve.

    Args:
        sat: Satellite to find the orbital path altitude of.
        geo: The GeoPosition which the satellite is viewed from.
        time: Time to find the orbital altitude.

    Returns:
        The maximum orbital path altitude in degrees.
    """

    state = sat.getState(time)
    ecc = sat.getTle().getEcc()
    sma = sat.getTle().getSma()

    # zenith vector for the GeoPosition
    zeta = norm(zenithVector(geo, time))
    # GeoPosition vector in geocentric reference frame
    gamma = geoPositionVector(geo, time)
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
    """
    Determines whether a satellite is eclipsed by the Earth.

    Args:
        satPos: Position vector of the satellite in kilometers.
        sunPos: Position vector of the Sun in kilometers.

    Returns:
        True if the satellite is eclipsed, false otherwise.
    """

    #   vectors are relative to the satellite
    earthPos = -satPos
    sunPos = -satPos + sunPos
    #   semi-diameters of earth and sun
    # todo: calculate real Earth radius from perspective here
    thetaE = asin(6371 / earthPos.mag())  # average earth radius
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


# noinspection PyUnresolvedReferences
class ShadowController:
    """
    Controller class for computing Earth shadow entrance and exit times. The class can be initialized with a TLE or
    a set of orbital elements, and then the computeValues() method can be called with any time to compute the shadow
    times. The class always computes the next shadow exit time, if the time parameter is between the entrance and exit
    times, the previous entrance time is computed, otherwise the next entrance time is computed.
    """

    # todo:
    '''
    when finding zeros, check if the phi's are not None (they've been calculated at least once). If they 
    aren't, then we don't need to do our complicated zero 'finder' method, since using the previous phi's will be very
    close and make a very good first guess. This doesn't just not waste time looking for unneeded zeros, but actually
    decreases the number of Newton iterations by a lot.
    '''

    def __init__(self, orbitData: TwoLineElement | OrbitalElements):
        """Initialize the controller with the object for computing element values, either a TLE or an element object.
        A TLE will yield much more accurate results as the elements can be updated with time."""

        if type(orbitData) == TwoLineElement:
            self._tle = orbitData
            self._elements = [None, None]
        elif type(orbitData) == OrbitalElements:
            self._tle = None
            self._elements = [orbitData, orbitData]
        else:
            raise TypeError('orbitData parameter must be of type TwoLineElement or OrbitalElements')
        self._s = [None, None]
        self._phi = [None, None]
        self._Re = [6371, 6371]
        self._jd = [None, None]
        self._phi0 = None

    def getAnomalies(self) -> tuple[float]:
        """Returns the most recently computed anomalies for the entrance and exit points of Earth's shadow."""
        return tuple(self._phi)

    def getTimes(self) -> tuple[JulianDate]:
        """Returns the most recently computed times the anomalies occur for the entrance and exit points of Earth's
        shadow."""
        return tuple(self._jd)

    def computeValues(self, time: JulianDate):
        """
        Computes/updates the anomalies and times of shadow entrance/exit. If the time is during transit of the shadow,
        the previous entrance anomaly and time will be computed, otherwise the next entrance and exit points will always
        be computed.

        Args:
            time: Time to compute the next(previous) anomalies.
        """

        self.__firstPass(time)
        self.__loop(time)

    def __loop(self, time: JulianDate):
        """
        Main loop which updates time dependent variables between iterations for a more precise answer.

        Args:
            time: Time to compute the next(previous) anomalies.

        Will set/update the elements, s-vectors, Res values and the phi and jd values.
        """

        for i in range(2):
            while True:
                self._elements[i].update(self._jd[i])
                self._s[i] = self.__computes(i)
                self._Re[i] = self.__computeRe(i)
                self._phi[i] = self.__computePhi(i)
                jd = self.__computeTimes(i, time)
                if abs(jd.difference(self._jd[i])) < (0.1 / 86400.0):  # 0.1 of a second
                    self._jd[i] = jd
                    break
                else:
                    self._jd[i] = jd

    def __firstPass(self, time: JulianDate):
        """
        Does the necessary 'first pass' initializations before an 'official' anomaly or time is computed.

        Args:
            time: Time to compute the next(previous) anomalies.

        Will set the jds, elements, s-vectors, and phis
        """

        self._jd = [time, time]
        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
            self._elements = [elements, elements]
        s = self.__computes(0)
        self._s = [s, s]
        self._phi[0] = self.__computePhi(0)
        self._phi[1] = self.__computePhi(1)
        self._phi0 = self._elements[0].trueAnomalyAt(time)
        self._jd[0] = self.__computeTimes(0, time)
        self._jd[1] = self.__computeTimes(1, time)

    def __computeTimes(self, index: int, time: JulianDate):
        """Computes the time to the specified anomaly signaled by the index value."""
        if index == 0 and self._phi[0] < self._phi0 <= self._phi[1]:
            return self._elements[0].timeToPrevTrueAnomaly(self._phi[0], time)
        else:
            return self._elements[index].timeToNextTrueAnomaly(self._phi[index], time)

    def __computePhi(self, index: int):
        """Computes the true anomaly by finding the zeros of Escobal's G function."""
        zeros = self.__getZeros(index)
        checks = [self._s[index][0] * cos(z) + self._s[index][1] * sin(z) for z in zeros]
        phis = [zeros[i] for i in range(4) if checks[i] > 0]
        gamma = self.__getGamma(self._s[index])
        if index == 0:
            if phis[0] + gamma < pi:
                return phis[0]
            else:
                return phis[1]
        else:
            if phis[0] + gamma < pi:
                return phis[1]
            else:
                return phis[0]

    def __getZeros(self, index: int):
        """Finds the zeros of Escobal's G function."""
        frac = 0.5
        ls = [self.__escobalNewtonMethod(pi * i * frac, index) for i in range(int(2 / frac))]
        for i in range(len(ls)):
            if ls[i] >= TWOPI:
                ls[i] %= TWOPI
        frac /= 2.0
        while len(ls) < 4:
            tmp = [self.__escobalNewtonMethod(pi * i * frac, index) for i in range(int(2 / frac))]
            for i in range(len(tmp)):
                if tmp[i] >= TWOPI:
                    tmp[i] %= TWOPI
            for i in tmp:
                for j in ls:
                    if abs(i - j) > 1e-5:
                        ls.append(i)
        return tuple(ls)

    def __escobalGFunction(self, phi: float, index: int) -> float:
        """Escobal's G function whose zeros represent the true anomalies of shadow entrance/exit."""
        ecc = self._elements[index].getEcc()
        c = (self._elements[index].getSma() * (1 - ecc * ecc)) ** 2
        term1 = self._Re[index] * self._Re[index] * ((1 + ecc * cos(phi)) ** 2)
        term2 = (-self._s[index][0] * cos(phi) - self._s[index][1] * sin(phi)) ** 2
        return term1 + c * term2 - c

    def __gPrime(self, phi: float, index: int) -> float:
        """Derivative of Escobal's G function, used while iterating Newton's method."""
        ecc = self._elements[index].getEcc()
        term1 = 2 * (1 + ecc * cos(phi)) * (-ecc * sin(phi))
        term2 = 2 * (-self._s[index][0] * cos(phi) - self._s[index][1] * sin(phi)) * (
                self._s[index][0] * sin(phi) - self._s[index][1] * cos(phi))
        return self._Re[index] * self._Re[index] * term1 + ((self._elements[index].getSma() * (
                1 - ecc * ecc)) ** 2) * term2

    def __escobalNewtonMethod(self, guess, index) -> float:
        """Solves Escobal's G function by iterating with the Newton-raphson method."""
        phi = guess
        gi = self.__escobalGFunction(phi, index)
        while abs(gi) > 1e-5:
            phi = phi - gi / self.__gPrime(phi, index)
            gi = self.__escobalGFunction(phi, index)

        '''gi = self.__escobalGFunction(guess, index)
        phi = phi - gi / self.__gPrime(phi, index)
        while abs(gi) > 1e-5:
            gi = self.__escobalGFunction(phi, index)
            phi = phi - gi / self.__gPrime(phi, index)'''
        return phi

    def __computeRe(self, index: float):
        """Computes the radius of the earth that the sun (dis)appears from."""
        a = self._elements[index].getSma()
        ecc = self._elements[index].getEcc()
        inc = self._elements[index].getInc()
        aop = self._elements[index].getAop()
        f = EARTH_FLATTENING
        ae = EARTH_EQUITORIAL_RADIUS
        phi = self._phi[index]

        Rz = (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
                sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - self._s[index][2] * (
                    self._s[index][0] * cos(phi) + self._s[index][1] * sin(phi)))
        fTerm = 2 * f - f * f
        mainTerm = ae * ae * (1 - fTerm)
        cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
        nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
        while abs(nextR - self._Re[index]) > 1e-7:
            self._Re[index] = nextR
            phi = self.__computePhi(index)
            Rz = (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
                    sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - self._s[index][2] * (
                        self._s[index][0] * cos(phi) + self._s[index][1] * sin(phi)))
            fTerm = 2 * f - f * f
            mainTerm = ae * ae * (1 - fTerm)
            cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
            nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
        return nextR

    @staticmethod
    def __getGamma(s):
        """Computes the gamma angle used in the algorithm, the angle between the Sun projection unit vector and the
        eccentric vector."""
        if s[0] < 0:
            if s[1] == 0:
                return 0
            elif s[1] > 0:
                return atan(-s[1] / s[0])
            else:
                return atan(s[1] / s[0])
        elif s[0] == 0:
            return pi / 2.0
        elif s[0] > 0:
            if s[1] > 0:
                return pi - atan(s[1] / s[0])
            elif s[1] == 0:
                return pi
            else:
                return pi - atan(-s[1] / s[0])

    def __computes(self, index: int):
        """Computes s, the unit projection vector of the Sun's position on the orbital plane."""
        S = -norm(getSunPosition2(self._jd[index]))
        return rotateOrderTo(ZXZ, EulerAngles(self._elements[index].getRaan(), self._elements[index].getInc(),
                                              self._elements[index].getAop()), S)
