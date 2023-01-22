import json
from math import radians, asin, cos, acos, sqrt, sin, degrees, pi, tan, atan
from operator import index

from pyevspace import Vector, vang, norm, cross, dot, ZYX, getMatrixEuler, Angles, rotateOffsetFrom, vxcl
from sattrack.sgp4 import computeEccentricVector

from sattrack._coordinates import _computePositionVector, _computeZenithVector
from sattrack._orbit import _trueAnomalyFromState, _trueToMeanAnomaly, \
    _smaToMeanMotion, _nearestTrueAnomaly, _vectorAlmostEqual
from sattrack._topocentric import _toTopocentricOffset, _toTopocentric
from sattrack.coordinates import GeoPosition
from sattrack.eclipse import getShadowTimes, Shadow
from sattrack.exceptions import PassConstraintException, NoPassException
from sattrack.orbit import Orbitable, Elements
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.sun import getSunTimes
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import atan3

__all__ = ('PositionInfo', 'SatellitePass', 'SatellitePassConstraints', 'getNextPass', 'getPassList', 'toTopocentric',
           'fromTopocentric', 'getAltitude', 'getAzimuth', 'azimuthAngleString', 'timeToHorizon')


def _azimuthAngleString(azimuth):
    """Converts an azimuth angle in degrees to a compass direction string."""
    if azimuth > 348.25 or azimuth <= 11.25:
        return 'N'
    elif azimuth <= 33.75:
        return 'NNE'
    elif azimuth <= 56.25:
        return 'NE'
    elif azimuth <= 78.75:
        return 'ENE'
    elif azimuth <= 101.25:
        return 'E'
    elif azimuth <= 123.75:
        return 'ESE'
    elif azimuth <= 146.25:
        return 'SE'
    elif azimuth <= 168.75:
        return 'SSE'
    elif azimuth <= 191.25:
        return 'S'
    elif azimuth <= 213.75:
        return 'SSW'
    elif azimuth <= 236.25:
        return 'SW'
    elif azimuth <= 258.75:
        return 'WSW'
    elif azimuth <= 281.25:
        return 'W'
    elif azimuth <= 303.75:
        return 'WNW'
    elif azimuth <= 326.25:
        return 'NW'
    else:
        return 'NNW'


def _checkNumericType(param, paramName: str):
    """Checks if a parameter type is numeric"""
    if not isinstance(param, (int, float)):
        raise TypeError(f'{paramName} parameter must be a numeric type', type(param))
    return param


def _checkJdType(param, paramName: str):
    """Checks if a parameter type is a JulianDate type."""
    if not isinstance(param, JulianDate):
        raise TypeError(f'{paramName} parameter must be a JulianDate type', type(param))
    return param


class PositionInfo:
    """A class to hold information regarding a single point in on overhead orbital pass."""

    __slots__ = '_altitude', '_azimuth', '_direction', '_time', '_illuminated', '_unobscured', '_visible'

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, illuminated: bool = False,
                 unobscured: bool = False):
        """Initializes the PositionInfo at a specific time with an altitude and azimuth in degrees. Illuminated refers
        to if the satellite is in sunlight and unobscured refers to the viewer being unobscured by sunlight from
        the view geo-position. Both values default to zero and the visible attribute is derived from a combination of
        these values."""

        self._altitude = _checkNumericType(altitude, 'altitude')
        self._azimuth = _checkNumericType(azimuth, 'azimuth')
        self._direction = _azimuthAngleString(azimuth)
        self._time = _checkJdType(time, 'time')
        self._illuminated = index(illuminated)
        self._unobscured = index(unobscured)
        self._visible = self._illuminated and self._unobscured

    def __str__(self):
        """Creates a string representation of the PositionInfo."""
        return ' {{:^17}} | {:^12} | {:^{altW}.{p}f} | {:^{azw}.{p}f} {:^5} | {!s:^11} | {!s:^10} | {!s:^7} ' \
            .format(self._time.time(), self._altitude, self._azimuth, '({})'.format(self._direction),
                    bool(self._illuminated), bool(self._unobscured), bool(self._visible), altW=8, azw=6, p=2)

    def __repr__(self):
        """Creates a string representation of the PositionInfo."""
        return f'PositionInfo({self._altitude}, {self._azimuth}, {repr(self._time)}, {self._illuminated},' \
               f'{self._unobscured})'

    def toJson(self):
        """Returns a string of the PositionInfo in json format."""
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        """Returns a dictionary of the PositionInfo to create json formats of other types containing a GeoPosition."""
        return {"altitude": self._altitude, "azimuth": self._azimuth, "direction": self._direction,
                "time": self._time.toDict(), "illuminated": self._illuminated, "unobscured": self._unobscured,
                "visible": self._visible}

    @property
    def altitude(self):
        """Return the altitude of the position in degrees.."""
        return self._altitude

    @property
    def azimuth(self):
        """Return the azimuth of the position in degrees."""
        return self._azimuth

    @property
    def direction(self):
        """Return the compass direction of the position as a string."""
        return self._direction

    @property
    def time(self):
        """Return the time of the position."""
        return self._time

    @property
    def illuminated(self):
        """Return a boolean value of if the position is illuminated by sunlight."""
        return self._illuminated

    @property
    def unobscured(self):
        """Return a boolean value of if the viewer of the position is unobscured by daylight."""
        return self._unobscured

    @property
    def visible(self):
        """Return if the position is visible i.e. illuminated and unobscured."""
        return self._visible


def _checkInfoType(param, paramName: str):
    """Checks if the param type is a PositionInfo."""
    if not isinstance(param, PositionInfo):
        raise TypeError(f'{paramName} parameter must be PositionInfo type')
    return param


def _checkInfoTypeNone(param, paramName: str):
    """Checks if the param type is a PositionInfo if it's not None."""
    if param is not None:
        return _checkInfoType(param, paramName)
    return param


class SatellitePass:
    """An object that holds multiple PositionInfos used to describe an entire satellite pass. All passes have a rise,
    set and maximum altitude PositionInfos. Some passes have other information significant to its viewing like first and
    last visible, first and last illuminated and first and last unobscured. These later PositionInfos are not required,
    and shouldn't be used when another PositionInfo shares the same information, e.g. the last visible PositionInfo
    shouldn't be set if it's the same as setInfo."""

    __slots__ = '_riseInfo', '_setInfo', '_maxInfo', '_firstUnobscured', '_lastUnobscured', '_firstIlluminated', \
                '_lastIlluminated', '_illuminated', '_unobscured', '_visible', '_name'

    def __init__(self, riseInfo: PositionInfo, setInfo: PositionInfo, maxInfo: PositionInfo, *,
                 firstUnobscuredInfo: PositionInfo = None, lastUnobscuredInfo: PositionInfo = None,
                 firstIlluminatedInfo: PositionInfo = None, lastIlluminatedInfo: PositionInfo = None, name: str = ''):
        """Initialize a SatellitePass with the respective PositionInfos. All non-required infos are keyword only and
        default to None."""

        self._riseInfo = _checkInfoType(riseInfo, 'riseInfo')
        self._setInfo = _checkInfoType(setInfo, 'setInfo')
        self._maxInfo = _checkInfoType(maxInfo, 'maxInfo')
        self._firstUnobscured = _checkInfoTypeNone(firstUnobscuredInfo, 'firstUnobscuredInfo')
        self._lastUnobscured = _checkInfoTypeNone(lastUnobscuredInfo, 'lastUnobscuredInfo')
        self._firstIlluminated = _checkInfoTypeNone(firstIlluminatedInfo, 'firstIlluminatedInfo')
        self._lastIlluminated = _checkInfoTypeNone(lastIlluminatedInfo, 'lastIlluminatedInfo')
        if firstIlluminatedInfo is not None or riseInfo.illuminated:
            self._illuminated = True
        else:
            self._illuminated = False
        if firstUnobscuredInfo is not None or riseInfo.unobscured:
            self._unobscured = True
        else:
            self._unobscured = False
        if self._riseInfo.visible \
                or (self._firstIlluminated is not None and self._firstIlluminated.visible) \
                or (self._firstUnobscured is not None and self._firstUnobscured.visible):
            self._visible = True
        else:
            self._visible = False
        self._name = name

    def _getInfos(self):
        """Returns the PositionInfos as a list of tuples with their names and values."""
        return [
            ('rise', self._riseInfo),
            ('set', self._setInfo),
            ('max', self._maxInfo),
            ('first unobscured', self._firstUnobscured),
            ('last unobscured', self._lastUnobscured),
            ('first illuminated', self._firstIlluminated),
            ('last illuminated', self._lastIlluminated)
        ]

    def __str__(self):
        """Returns a human readably text table with the pass information."""
        width = 97
        if self._name == '':
            title = '{:^{w}}'.format('Pass details at {}'.format(self._maxInfo.time.date()), w=width)
        else:
            title = '{:^{w}}'.format('Pass details for {}, at {}'.format(self._name, self._maxInfo.time.date()),
                                     w=width)
        heading = ' {:^17} | {:^12} | altitude | {:^12} | illuminated | unobscured | visible ' \
            .format('instance', 'time', 'azimuth')
        lineBreak = '-' * width
        string = title + '\n' + heading + '\n'

        infos = [i for i in self._getInfos() if i[1] is not None]
        infos.sort(key=lambda o: o[1].time.value)

        for name, info in infos:
            string += '{}\n{}\n'.format(lineBreak, str(info).format(name))
        return string

    def __repr__(self):
        """Creates a string representation of the pass."""
        return self.__str__()

    def toJson(self):
        """Returns a string of the pass in json format."""
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        """Returns a dictionary of the pass to create json formats of other types containing a GeoPosition."""
        return {"name": self._name, "riseInfo": self._riseInfo, "setInfo": self._setInfo, "maxInfo": self._maxInfo,
                "firstUnobscured": self._firstUnobscured, "lastUnobscured": self._lastUnobscured,
                "firstIlluminated": self._firstIlluminated, "lastIlluminated": self._lastIlluminated,
                "illuminated": self._illuminated, "unobscured": self._unobscured, "visible": self._visible}

    @property
    def riseInfo(self):
        """Returns the PositionInfo for the rise time of the pass."""
        return self._riseInfo

    @property
    def setInfo(self):
        """Returns the PositionInfo for the set time of the pass."""
        return self._setInfo

    @property
    def maxInfo(self):
        """Returns the PositionInfo for the time the pass reaches maximum altitude."""
        return self._maxInfo

    @property
    def firstUnobscuredInfo(self):
        """Returns the PositionInfo for the first unobscured time of the pass."""
        return self._firstUnobscured

    @property
    def lastUnobscuredInfo(self):
        """Returns the PositionInfo for the last unobscured time of the pass."""
        return self._lastUnobscured

    @property
    def firstIlluminatedInfo(self):
        """Returns the PositionInfo for the first illuminated time of the pass."""
        return self._firstIlluminated

    @property
    def lastIlluminatedInfo(self):
        """Returns the PositionInfo for the last illuminated time of the pass."""
        return self._lastIlluminated

    @property
    def unobscured(self):
        """Returns True if the pass is unobscured by daylight, False otherwise."""
        return self._unobscured

    @property
    def illuminated(self):
        """Returns True if the satellite is illuminated by sunlight during the pass, False otherwise."""
        return self._illuminated

    @property
    def visible(self):
        """Returns True if the satellite pass is visible, False otherwise."""
        return self._visible


def _checkAltitudes(minAlt, maxAlt):
    """Checks if minimum and maximum altitude constraints are valid."""

    if minAlt is not None and maxAlt is not None and minAlt > maxAlt:
        raise PassConstraintException('Cannot set a minimum altitude constraint less than a maximum altitude '
                                      'constraint.')
    if minAlt is not None and minAlt > 90:
        raise ValueError('minAltitude cannot be above 90 degrees')
    if maxAlt is not None and maxAlt < 0:
        raise ValueError('maxAltitude cannot be below 0 degrees')


def _checkDuration(minDur, maxDur):
    """Checks if minimum and maximum pass duration constraints are valid."""

    if minDur is not None and maxDur is not None and minDur > maxDur:
        raise PassConstraintException('Cannot set a minimum duration constraint longer than a maximum duration'
                                      'constraint.')
    if maxDur is not None and maxDur < 0:
        raise ValueError('maxDuration cannot be less than 0')


class SatellitePassConstraints:
    """A set of constraints for filtering satellite passes based on certain criteria. The possible constraints are
    minium or maximum altitudes or durations, and whether or not passes are illuminated, unobscured or visible."""

    __slots__ = '_minAltitude', '_maxAltitude', '_minDuration', '_maxDuration', '_illuminated', '_unobscured', \
                '_visible'

    def __init__(self, *, minAltitude: float = None, maxAltitude: float = None, minDuration: float = None,
                 maxDuration: float = None, illuminated: bool = None, unobscured: bool = None, visible: bool = None):
        """Initializes the SatellitePassConstraints object with the given constraints, all of which are keyword only.
        Altitudes are in degrees and durations are in minutes."""

        _checkNumericType(minAltitude, 'minAltitude')
        _checkNumericType(maxAltitude, 'maxAltitude')
        _checkAltitudes(minAltitude, maxAltitude)
        _checkDuration(minDuration, maxDuration)
        self._minAltitude = minAltitude
        self._maxAltitude = maxAltitude
        self._minDuration = minDuration
        self._maxDuration = maxDuration
        self._illuminated = index(illuminated)
        self._unobscured = index(unobscured)
        self._visible = index(visible)

    @property
    def minAltitude(self):
        """Returns the minimum altitude constraint in degrees."""
        return self._minAltitude

    @minAltitude.setter
    def minAltitude(self, value):
        """Sets the minimum altitude constraint in degrees."""
        _checkNumericType(value, 'value')
        _checkAltitudes(value, self._maxAltitude)
        self._minAltitude = value

    @minAltitude.deleter
    def minAltitude(self):
        """Removes the minimum altitude constraint."""
        self._minAltitude = None

    @property
    def maxAltitude(self):
        """Returns the maximum altitude constraint in degrees."""
        return self._maxAltitude

    @maxAltitude.setter
    def maxAltitude(self, value):
        """Sets the maximum altitude constraint in degrees."""
        _checkNumericType(value, 'value')
        _checkAltitudes(self._minAltitude, value)
        self._maxAltitude = value

    @maxAltitude.deleter
    def maxAltitude(self):
        """Removes the maximum altitude constraint."""
        self._maxAltitude = None

    @property
    def minDuration(self):
        """Returns the minimum pass duration in minutes."""
        return self._minDuration

    @minDuration.setter
    def minDuration(self, value):
        """Sets the minimum pass duration in minutes."""
        _checkNumericType(value, 'value')
        _checkDuration(value, self._maxDuration)
        self._minDuration = value

    @minDuration.deleter
    def minDuration(self):
        """Removes the minimum pass duration constraint."""
        self._minDuration = None

    @property
    def maxDuration(self):
        """Returns the maximum pass duration in minutes."""
        return self._maxDuration

    @maxDuration.setter
    def maxDuration(self, value):
        """Sets the maximum pass duration in minutes."""
        _checkNumericType(value, 'value')
        _checkDuration(self._minDuration, value)
        self._maxDuration = value

    @maxDuration.deleter
    def maxDuration(self):
        """Removes the maximum pass duration constraint."""
        self._maxDuration = None

    @property
    def illuminated(self):
        """Returns the illuminated pass constraint."""
        return self._illuminated

    @illuminated.setter
    def illuminated(self, value):
        """Sets the illuminated pass constraint."""
        self._illuminated = index(value)

    @illuminated.deleter
    def illuminated(self):
        """Removes the illuminated pass constraint."""
        self._illuminated = None

    @property
    def unobscured(self):
        """Returns the unobscured pass constraint."""
        return self._unobscured

    @unobscured.setter
    def unobscured(self, value):
        """Sets the unobscured pass constraint."""
        self._unobscured = index(value)

    @unobscured.deleter
    def unobscured(self):
        """Removes the unobscured pass constraint."""
        self._unobscured = None

    @property
    def visible(self):
        """Returns the visible pass constraint."""
        return self._visible

    @visible.setter
    def visible(self, value):
        """Sets the visible pass constraint."""
        self._visible = index(value)

    @visible.deleter
    def visible(self):
        """Removes the visible pass constraint."""
        self._unobscured = None


def __maxPassAnomalies(position, velocity, mu, pVector, ecc):
    """Compute the anomalies as they are used in the max pass finder methods."""
    # eccentricVector = _compute_eccentric_vector(position, velocity, mu)
    eccentricVector = computeEccentricVector(position, velocity, mu)
    trueAnomaly_n = _trueAnomalyFromState(position, velocity, mu)
    # todo: use eccentricVector.mag() if we know it's accurate
    meanAnomaly_n = _trueToMeanAnomaly(trueAnomaly_n, ecc)
    trueAnomaly_n1 = vang(eccentricVector, pVector)
    if norm(cross(eccentricVector, pVector)) != norm(cross(position, velocity)):
        trueAnomaly_n1 = TWOPI - trueAnomaly_n1
    meanAnomaly_n1 = _trueToMeanAnomaly(trueAnomaly_n1, ecc)
    return trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1


def _getPVector(geo: GeoPosition, position: Vector, velocity: Vector, time: JulianDate) -> Vector:
    """Computes the 'P-Vector', the shortest vector between a geo-position and the line formed by the intersection of
    the horizontal plane of the geo-position and the orbital plane. This is the vector that points towards the highest
    point of the orbital path if it's above a viewer's horizon."""
    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    elevation = geo.elevation

    # zenith vector for geo-position
    zeta = _computeZenithVector(latitude, longitude, elevation, time)
    # position vector for geo-position
    gamma = _computePositionVector(latitude, longitude, elevation, time)
    # normalized angular momentum, vector equation for orbital plane
    lamb = norm(cross(position, velocity))

    # compute intermediate values to find solution to parameterized vector intersection
    if lamb[1] != 0:
        x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
        y = -lamb[0] * x / lamb[1]
    else:
        y = dot(zeta, gamma) / (zeta[1] - (zeta[0] * lamb[1] / lamb[0]))
        x = -lamb[1] * y / lamb[0]
    r = Vector(x, y, 0)
    v = cross(zeta, lamb)

    # compute exact solution for parameterized vector which yields the nearest point to intersection
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    return r + v * t


def _orbitAltitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    """Computes the altitude of the orbital path above or below the geo-position's horizontal reference frame."""

    state = satellite.getState(time)
    elements = satellite.getElements(time)

    gamma = _computePositionVector(radians(geo.latitude), radians(geo.longitude), geo.elevation, time)
    pVector = _getPVector(geo, *state, time)

    eccentricVector = computeEccentricVector(*state, satellite.body.mu)
    trueAnomaly = vang(eccentricVector, pVector)
    pSat = norm(pVector) * ((elements.sma * (1 - elements.ecc * elements.ecc)) / (1 + elements.ecc * cos(trueAnomaly)))

    angle = vang(pVector - gamma, pSat - gamma)
    return angle if pSat.mag2() > pVector.mag2() else -angle


def _timeToPlane(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """Estimates the time it takes for the orbital path to be level with a geo-position's horizontal reference frame,
    which begins the opportunity for overhead satellite passes."""

    jd = time
    alt = _orbitAltitude(satellite, geo, time)
    if alt > 0:
        return jd

    while abs(alt) > 2.7e-4:  # one arc-second
        # todo: find an analytical guess of this dt
        dt = -alt / 360
        jd = jd.future(dt)
        alt = _orbitAltitude(satellite, geo, jd)
    return jd


# def _nextPassMaxApprox(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
#     """Approximates the time of the next satellite pass. For circular orbits this usually approximates the point the
#     satellite achieves its maximum altitude above the horizon."""
#
#     # if the orbital plane is not above, find the time it appears overhead to start looking for passes
#     if _orbitAltitude(satellite, geo, time) < 0:
#         tn = _timeToPlane(satellite, geo, time)
#     else:
#         tn = time
#
#     # get variables needed for computations
#     mu = satellite.body.mu
#     elements = satellite.getElements(time)
#     ecc = elements.ecc
#     # mean motion in radians / day
#     meanMotion = _smaToMeanMotion(elements.sma, mu) * 86400
#
#     # find the first guess of when the satellites position occupies the 'p-vector'
#     state = satellite.getState(tn)
#     pVector = _getPVector(geo, *state, tn)
#     _, _, meanAnomaly_n, meanAnomaly_n1 = __maxPassAnomalies(*state, mu, pVector, ecc)
#     dma = meanAnomaly_n1 - meanAnomaly_n
#     if meanAnomaly_n1 < meanAnomaly_n:
#         dma += TWOPI
#     tn = tn.future(dma / meanMotion)
#
#     # update state and pVector for the while loop condition, this should be a do-while =(
#     state = satellite.getState(tn)
#     pVector = _getPVector(geo, *state, tn)
#     while vang(state[0], pVector) > 4.85e-6:  # one arc-second
#         trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1 = __maxPassAnomalies(*state, mu, pVector, ecc)
#         dma = meanAnomaly_n1 - meanAnomaly_n
#         if trueAnomaly_n1 <= trueAnomaly_n:
#             if (trueAnomaly_n - trueAnomaly_n1) >= pi:
#                 dma += TWOPI
#         else:
#             if (trueAnomaly_n1 - trueAnomaly_n) > pi:
#                 dma -= TWOPI
#         tn = tn.future(dma / meanMotion)
#         state = satellite.getState(tn)
#         pVector = _getPVector(geo, *state, tn)
#
#     return tn


# def _vector_almost_equal(vec1, vec2, epsilon=1e-5):
#     return dot(vec1, vec2) > (1 - epsilon)

    # for v1, v2 in zip(vec1, vec2):
    #     if abs(v1 - v2) > epsilon:
    #         return False
    #
    # return True


def _nextPassMaxApproxFix(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """Approximates the time of the next satellite pass. For circular orbits this usually approximates the point the
    satellite achieves its maximum altitude above the horizon."""
    '''
    move the time forward until geo is under the plane
    
    need to find time (angle) until sat occupies the pVector position moving forward in time
        - find eccentric vector
            - gotten from state vectors
        - *this only works if a change in true anomaly is constant change in mean anomaly*
          find angle between position vector to pVector moving forward in time
            - compute angular momentum vector (normalize it)
            - find angle between pVector and current position vector
                - use vang() to get angle.
                - if norm of cross product of pVector to position is the same direction
                  as angular momentum angle to travel is TWOPI - angle from pVector to position
                - else (norm of cross of pVector to position is opposite of angular momentum)
                  angle to travel is angle from position to pVector
            - MAKE SURE!!!! this works for prograde AND retrograde orbits
        - find the current anomaly and anomaly at pVector position
        - find the angle between these two angle (a0, a1)
            - a1 - a0, if answer is < 0 just add TWOPI to it
        - find the time to travel this angle
            - compute mean motion in radians / days
            - divide the angle by mean motion for time to occupy pVector in solar days
        - extra: account for the earth's rotation while satellite moves to pVector position
            - would need to understand how pVector changes relative to surface
    '''

    # todo: how do we incorporate earth's rotation/what happens when we do?

    # if the orbital plane is not above, find the time it appears overhead to start looking for passes
    if _orbitAltitude(satellite, geo, time) < 0:
        tn = _timeToPlane(satellite, geo, time)
    else:
        tn = time

    # get variables needed for computations
    mu = satellite.body.mu
    state = satellite.getState(tn)
    elements = Elements.fromState(*state, tn)
    # mean motion in radians / day
    meanMotion = _smaToMeanMotion(elements.sma, mu) * 86400
    pVector = _getPVector(geo, *state, tn)
    eccVector = computeEccentricVector(*state, mu)
    angularMomentumNorm = norm(cross(*state))

    # find the first guess of when the satellites position occupies the 'p-vector'
    # todo: replace the methods when the vector equality works
    taPVector = vang(eccVector, pVector)
    if not _vectorAlmostEqual(norm(cross(eccVector, pVector)), angularMomentumNorm):
        taPVector = TWOPI - taPVector

    maPVector = _trueToMeanAnomaly(taPVector, elements.ecc)
    maPosition = elements.meanAnomaly

    dma = maPVector - maPosition
    if dma < 0:
        dma += TWOPI

    tn = tn.future(dma / meanMotion)

    # update the state and pVector for the while loop condition, this should be a do-while =(
    state = satellite.getState(tn)
    pVector = _getPVector(geo, *state, tn)
    while vang(state[0], pVector) > 4.58e-6:  # 1 arc second
        elements = Elements.fromState(*state, tn)
        # mean motion in radians / day
        meanMotion = _smaToMeanMotion(elements.sma, mu) * 86400
        eccVector = computeEccentricVector(*state, mu)
        angularMomentumNorm = norm(cross(*state))

        # find the position when the satellite occupies the 'p-vector'
        taPVector = vang(eccVector, pVector)
        if not _vectorAlmostEqual(norm(cross(eccVector, pVector)), angularMomentumNorm):
            taPVector = TWOPI - taPVector

        maPVector = _trueToMeanAnomaly(taPVector, elements.ecc)
        maPosition = elements.meanAnomaly

        # adjust the time based on the anomalies found above
        dma = maPVector - maPosition
        if dma < 0:
            dma += TWOPI

        tn = tn.future(dma / meanMotion)
        state = satellite.getState(tn)
        pVector = _getPVector(geo, *state, tn)

    return tn


def _computeAltitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    """Computes the altitude in radians of the highest point of an orbital path above a geo-position's horizontal
    reference frame."""

    position = satellite.getState(time)[0]
    topocentricPosition = _toTopocentricOffset(position, geo, time)
    return asin(topocentricPosition[2] / topocentricPosition.mag())


def _getAltitudeDerivative(topoPosition: Vector, topoVelocity: Vector) -> float:
    """Computes the derivative of altitude with respect to time from topocentric position and velocity vectors."""

    rDotV = dot(topoPosition, topoVelocity)
    rMag2 = topoPosition.mag2()

    rightNumeratorTerm = topoPosition[2] * rDotV / rMag2
    rightDenominatorTerm = topoPosition[2] * topoPosition[2] / rMag2
    numerator = topoVelocity[2] - rightNumeratorTerm
    denominator = sqrt(1 - rightDenominatorTerm)

    return numerator / denominator


def _getAltitude(topoPosition: Vector) -> float:
    """Computes the altitude of a topocentric position vector."""
    return asin(topoPosition[2] / topoPosition.mag())


def _maxPassRefine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """Refines the original time approximation of the maximum altitude of a satellite pass."""

    # this works but the times are inaccurate, either the equation implemented in _getAltitudeDerivative is wrong,
    # or maybe the earth rotation needs to be accounted for.
    '''state = satellite.getState(time)
    topoPosition = _toTopocentricOffset(state[0], geo, time)
    topoVelocity = _toTopocentric(state[1], geo, time)
    dadt = _getAltitudeDerivative(topoPosition, topoVelocity)

    # move forward or back 1 second and see how what the derivative is
    # todo: can we analytical explain why 1 second should be used here?
    if dadt < 0:
        direction = -1
    elif dadt > 0:
        direction = 1
    else:
        # if dadt is zero where at maximum altitude
        return time

    dadt2 = direction
    time2 = time
    topoPosition2 = 0   # to please the IDE
    topoVelocity2 = 0   # to please the IDE
    # keep moving forward or back until time and time2 are not both increasing or both decreasing
    while dadt * dadt2 > 0:
        time2 = time2.future(direction / 86400)
        state2 = satellite.getState(time2)
        topoPosition2 = _toTopocentricOffset(state2[0], geo, time2)
        topoVelocity2 = _toTopocentric(state2[1], geo, time2)
        dadt2 = _getAltitudeDerivative(topoPosition2, topoVelocity2)

    # use bisecting logic to find the maximum altitude
    if time < time2:
        beforeTime = time
        afterTime = time2
    else:
        beforeTime = time2
        afterTime = time
        beforeDadt = _getAltitudeDerivative(topoPosition2, topoVelocity2)
        afterDadt = _getAltitudeDerivative(topoPosition, topoVelocity)
    while (afterTime - beforeTime) > (0.001 / 86400):     # 1/1000th of a second
        halftime = (afterTime - beforeTime) / 2
        middleTime = beforeTime.future(halftime)
        state = satellite.getState(middleTime)
        topoPos = _toTopocentricOffset(state[0], geo, middleTime)
        topoVel = _toTopocentric(state[1], geo, middleTime)
        dadt = _getAltitudeDerivative(topoPos, topoVel)
        if dadt < 0:
            afterTime = middleTime
        elif dadt > 0:
            beforeTime = middleTime
        else:
            return middleTime

    return beforeTime'''

    # todo: utilize the dadt values to future time increase and make this a more analytical guess
    oneSecond = 1 / 86400
    altitude = _computeAltitude(satellite, geo, time)
    futureAltitude = _computeAltitude(satellite, geo, time.future(oneSecond))
    pastAltitude = _computeAltitude(satellite, geo, time.future(-oneSecond))

    parity = 0  # to please the editor
    if altitude >= futureAltitude and altitude >= pastAltitude:
        return time
    elif altitude < futureAltitude:
        parity = 1
    elif altitude < pastAltitude:
        parity = -1

    jd = time
    nextAltitude = _computeAltitude(satellite, geo, time.future(parity * oneSecond))
    while altitude < nextAltitude:
        jd = jd.future(parity * oneSecond)
        altitude = _computeAltitude(satellite, geo, jd)
        nextAltitude = _computeAltitude(satellite, geo, jd.future(parity * oneSecond))

    return jd


def _nextPassMax(satellite: Orbitable, geo: GeoPosition, time: JulianDate, timeout: float) -> JulianDate:
    """Computes the time of maximum altitude of the next satellite pass."""

    nextMax = _nextPassMaxApproxFix(satellite, geo, time)
    altitude = _orbitAltitude(satellite, geo, nextMax)
    while altitude < 0:
        if nextMax - time < timeout:
            nextMax = _nextPassMaxApproxFix(satellite, geo, nextMax.future(0.001))
            altitude = _orbitAltitude(satellite, geo, nextMax)
        else:
            # todo: make a more specific exception type for this
            raise NoPassException('timed out looking for next pass')
    return _maxPassRefine(satellite, geo, nextMax)


def _riseSetTimesApprox(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    """Computes the approximate rise and set times of a satellite above a geo-position's horizon."""

    elements = satellite.getElements(time)
    mu = satellite.body.mu
    ecc = elements.ecc
    a = elements.sma
    c = a * ecc
    b = sqrt(a * a - c * c)

    state = satellite.getState(time)
    hNorm = norm(cross(*state))
    u = norm(computeEccentricVector(*state, mu)) * a
    v = norm(cross(hNorm, u)) * b
    ce = -norm(u) * c

    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    elevation = geo.elevation
    zeta = _computeZenithVector(latitude, longitude, elevation, time)
    gamma = _computePositionVector(latitude, longitude, elevation, time)

    dotZU = dot(zeta, u)
    dotZV = dot(zeta, v)
    R = sqrt(dotZU * dotZU + dotZV * dotZV)
    try:
        arg = dot(zeta, gamma - ce) / R
        beta = acos(arg)
    except ValueError:
        raise NoPassException(f'No pass during this time: {time}.')
    alpha = atan3(dotZV, dotZU)
    w1 = alpha + beta
    w2 = alpha - beta

    rho1 = (u * cos(w1) + v * sin(w1)).mag()
    trueAnomaly_1 = atan3(rho1 * sin(w1), rho1 * cos(w1) - a * ecc)
    rho2 = (u * cos(w2) + v * sin(w2)).mag()
    trueAnomaly_2 = atan3(rho2 * sin(w2), rho2 * cos(w2) - a * ecc)
    # todo: find out why we did this?
    trueAnomaly_10 = _trueAnomalyFromState(*state, mu)
    trueAnomaly_20 = _trueAnomalyFromState(*state, mu)
    # meanMotion in rev / day
    meanMotion = _smaToMeanMotion(a, mu) * 86400 / TWOPI
    jd1 = _nearestTrueAnomaly(meanMotion, ecc, trueAnomaly_10, time, trueAnomaly_1)
    jd2 = _nearestTrueAnomaly(meanMotion, ecc, trueAnomaly_20, time, trueAnomaly_2)
    if jd1 < jd2:
        return jd1, jd2
    return jd2, jd1


def __timeToHorizonIteration(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    state = satellite.getState(time)
    topoPosition = _toTopocentricOffset(state[0], geo, time)

    vel = geo.getRadius() * TWOPI / satellite.body.period
    geoVelocity = Vector(0, vel, 0)
    topoVelocity = _toTopocentric(state[1], geo, time) - geoVelocity

    velExclude = vxcl(topoVelocity, topoPosition)

    a = asin(topoPosition[2] / topoPosition.mag())
    C = pi / 2
    # use topoVelocity dot Vector.e3 to get sign of Vector.e3 here
    B = vang(velExclude, -Vector.e3)

    # using cotangent four part formulae
    temp = (cos(a) * cos(B) + (1 / tan(C)) * sin(B)) / sin(a)
    c = atan(1 / temp)

    alpha = velExclude.mag() / topoPosition.mag()
    dt = c / alpha

    return time.future(dt / 86400)


def _horizonTimeRefine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    updatedTime = time
    alt = getAltitude(satellite, geo, updatedTime)
    while alt > 4.848e-6:   # 1 arc-second
        updatedTime = __timeToHorizonIteration(satellite, geo, updatedTime)
        alt = getAltitude(satellite, geo, updatedTime)

    return updatedTime


def _riseSetTimes(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    """Computes the rise and set time of the satellite pass where time is a time during the pass."""

    riseTime, setTime = _riseSetTimesApprox(satellite, geo, time)
    riseTime = _horizonTimeRefine(satellite, geo, riseTime)
    setTime = _horizonTimeRefine(satellite, geo, setTime)
    return riseTime, setTime


def _getSpecialTimes(riseTime: JulianDate, setTime: JulianDate, beginTime: JulianDate, finishTime: JulianDate) \
        -> (JulianDate, JulianDate, bool, bool):
    """Computes the first and last times for the 'special' events in a pass (illuminated or unobscured)."""
    firstTime = riseTime
    lastTime = setTime
    #   rise and set sandwich begin time
    if riseTime.value < beginTime.value < setTime.value:
        isFirst = True
        isLast = False
        lastTime = beginTime
    #   rise and set sandwich finish time
    elif riseTime.value < finishTime.value < setTime.value:
        isFirst = False
        firstTime = finishTime
        isLast = True
    #   rise and set are between begin and finish
    elif beginTime.value < riseTime.value < finishTime.value \
            and beginTime.value < setTime.value < finishTime.value:
        isFirst = False
        firstTime = None
        isLast = False
        lastTime = None
    #   rise and set are before begin
    else:
        isFirst = True
        isLast = True

    return firstTime, lastTime, isFirst, isLast


def _deriveBasicInfo(satellite: Orbitable, geo: GeoPosition, time: JulianDate, isIlluminated: bool,
                     isUnobscured: bool) -> PositionInfo:
    """Creates a PositionInfo object for a given event during a pass."""

    position = satellite.getState(time)[0]
    topocentricPosition = _toTopocentricOffset(position, geo, time)
    altitude = degrees(asin(topocentricPosition[2] / topocentricPosition.mag()))
    azimuth = degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))
    return PositionInfo(altitude, azimuth, time, isIlluminated, isUnobscured)


def getNextPass(satellite: Orbitable, geo: GeoPosition, timeOrPass: JulianDate | SatellitePass,
                timeout: float = 7) -> SatellitePass:
    """Computes the next satellite pass after a given time or another SatellitePass object. In the case of the latter
    the pass times in the SatellitePass are used to find the next pass. For multiple passes over a specific duration
    see getPassList()."""

    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if isinstance(timeOrPass, JulianDate):
        time = timeOrPass
    elif isinstance(timeOrPass, SatellitePass):
        time = timeOrPass.maxInfo.time.future(0.001)
    else:
        raise TypeError('timeOrPass parameter must be JulianDate or SatellitePass type')
    nextPassTime = _nextPassMax(satellite, geo, time, timeout)

    riseTime, setTime = _riseSetTimes(satellite, geo, nextPassTime)
    enterTime, exitTime = getShadowTimes(satellite, nextPassTime, Shadow.PENUMBRA)
    # todo: if timezone of riseTime is different than the timezone of geo, then the rise/set times may be for the
    #   wrong day and the unobscured values may be wrong.
    # todo: idea: check sun altitude and always find the next sunset based on the angle, then find the previous rise,
    #   then the above issue shouldn't exist
    sunRiseTime, sunSetTime = getSunTimes(riseTime, geo)

    firstIlluminatedTime, lastIlluminatedTime, riseIlluminated, setIlluminated \
        = _getSpecialTimes(riseTime, setTime, enterTime, exitTime)
    firstUnobscuredTime, lastUnobscuredTime, riseUnobscured, setUnobscured \
        = _getSpecialTimes(riseTime, setTime, sunRiseTime, sunSetTime)
    maxIlluminated = not (enterTime <= nextPassTime <= exitTime)
    maxUnobscured = (nextPassTime < sunRiseTime) or (nextPassTime >= sunSetTime)

    riseInfo = _deriveBasicInfo(satellite, geo, riseTime, riseIlluminated, riseUnobscured)
    setInfo = _deriveBasicInfo(satellite, geo, setTime, setIlluminated, setUnobscured)
    maxInfo = _deriveBasicInfo(satellite, geo, nextPassTime, maxIlluminated, maxUnobscured)

    firstIlluminatedInfo = lastIlluminatedInfo = firstUnobscuredInfo = lastUnobscuredInfo = None
    if firstIlluminatedTime is not None and firstIlluminatedTime != riseTime:
        firstIlluminatedInfo = _deriveBasicInfo(satellite, geo, firstIlluminatedTime, True,
                                                firstIlluminatedTime < sunRiseTime
                                                or firstIlluminatedTime >= sunSetTime)
    if lastIlluminatedTime is not None and lastIlluminatedTime != setTime:
        lastIlluminatedInfo = _deriveBasicInfo(satellite, geo, lastIlluminatedTime, True,
                                               lastIlluminatedTime < sunRiseTime
                                               or lastIlluminatedTime >= sunSetTime)
    if firstUnobscuredTime is not None and firstUnobscuredTime == sunSetTime:
        firstUnobscuredInfo = _deriveBasicInfo(satellite, geo, firstUnobscuredTime,
                                               firstUnobscuredTime < enterTime
                                               or firstUnobscuredTime >= exitTime, True)
    if lastUnobscuredTime is not None and lastUnobscuredTime == sunRiseTime:
        lastUnobscuredInfo = _deriveBasicInfo(satellite, geo, lastUnobscuredTime,
                                              lastUnobscuredTime < enterTime
                                              or lastUnobscuredTime >= exitTime, True)

    np = SatellitePass(riseInfo, setInfo, maxInfo,
                       firstUnobscuredInfo=firstUnobscuredInfo,
                       lastUnobscuredInfo=lastUnobscuredInfo,
                       firstIlluminatedInfo=firstIlluminatedInfo,
                       lastIlluminatedInfo=lastIlluminatedInfo)

    # todo: filter the pass
    return np


# todo: pass a PassConstraints into this?
def getPassList(satellite: Orbitable, geo: GeoPosition, start: JulianDate, duration: float) -> tuple[SatellitePass]:
    """Generate a list of SatellitePass objects over a duration of time in solar days."""

    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype', type(satellite))
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type', type(geo))
    _checkJdType(start, 'start')
    _checkNumericType(duration, 'duration')
    if duration < 0:
        raise ValueError('duration parameter must not be negative', duration)

    passList = []
    # offset the initial small step forward between computations
    nextTime = start.future(-0.001)

    remainingTime = duration + 0.001
    while remainingTime > 0:
        try:
            nextPass = getNextPass(satellite, geo, nextTime.future(0.001), remainingTime)
        except NoPassException:
            break
        passList.append(nextPass)
        nextTime = nextPass.maxInfo.time
        remainingTime = duration - (nextTime - start)
    return tuple(passList)


def toTopocentric(vector: Vector, geo: GeoPosition, time: JulianDate) -> Vector:
    """Converts a vector to a topocentric reference frame defined by a geo-position at a specified time."""

    if not isinstance(vector, Vector):
        raise TypeError('vector parameter must be EVector type')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    return _toTopocentricOffset(vector, geo, time)


def fromTopocentric(vector: Vector, geo: GeoPosition, time: JulianDate) -> Vector:
    """Converts a vector from a topocentric reference frame defined by a geo-position at a specified time."""
    if not isinstance(vector, Vector):
        raise TypeError('vector parameter must be EVector type')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    geoVector = geo.getPositionVector(time)
    matrix = getMatrixEuler(
        ZYX,
        Angles(
            radians(geo.longitude) + earthOffsetAngle(time),
            radians(90 - geo.latitude),
            0.0
        )
    )

    return rotateOffsetFrom(matrix, geoVector, vector)


def getAltitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    """Returns the altitude of an orbitable from a geo-position at a specified time."""

    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    position = satellite.getState(time)[0]
    topocentricPosition = _toTopocentricOffset(position, geo, time)
    arg = topocentricPosition[2] / topocentricPosition.mag()
    return degrees(asin(arg))


def getAzimuth(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    """Returns the azimuth of an orbitable from a geo-position at a specified time."""

    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    position = satellite.getState(time)[0]
    topocentricPosition = _toTopocentricOffset(position, geo, time)
    return degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))


def azimuthAngleString(azimuth: float) -> str:
    """Returns a compass direction string from an azimuth in degrees."""

    # azimuth in degrees
    _checkNumericType(azimuth, 'azimuth')
    if not (0 <= azimuth <= 360):
        raise ValueError('azimuth angle must be between 0 and 360 degrees')
    return _azimuthAngleString(azimuth)

# todo: create alt-az type and methods for it


def timeToHorizon(sat: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    state = sat.getState(time)
    topoPosition = _toTopocentricOffset(state[0], geo, time)
    topoVelocity = _toTopocentric(state[1], geo, time)
    velExclude = vxcl(topoVelocity, topoPosition)

    a = asin(topoPosition[2] / topoPosition.mag())
    C = pi / 2
    # use topoVelocity dot Vector.e3 to get sign of Vector.e3 here
    B = vang(velExclude, -Vector.e3)

    # using cotangent four part formulae
    temp = (cos(a) * cos(B) + (1 / tan(C)) * sin(B)) / sin(a)
    c = atan(1 / temp)

    alpha = velExclude.mag() / topoPosition.mag()
    dt = c / alpha

    return time.future(dt / 86400)
