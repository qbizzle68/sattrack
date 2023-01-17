import json
from math import radians, pi, asin, cos, acos, sqrt, sin, degrees
from operator import index

from pyevspace import Vector, vang, norm, cross, dot, ZYX, getMatrixEuler, Angles, rotateEulerTo, rotateOffsetFrom
from sattrack._topocentric import _to_topocentric
from sattrack.eclipse import getShadowTimes, Shadow
from sattrack.exceptions import PassConstraintException, NoPassException
# from sattrack.rotation.order import ZYX
# from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateOrderTo, \
#     undoRotateToThenOffset
from sattrack.sun import getSunTimes
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack._coordinates import _compute_position_vector, _compute_zenith_vector
from sattrack._orbit import _true_anomaly_from_state, _true_to_mean_anomaly, \
    _sma_to_mean_motion, _nearest_true_anomaly
from sattrack.coordinates import GeoPosition
from sattrack.orbit import Orbitable, Elements
from sattrack.util.constants import TWOPI, EARTH_MU
from sattrack.util.conversions import atan3
from sattrack.sgp4 import computeEccentricVector, elementsFromState

__all__ = ('PositionInfo', 'SatellitePass', 'SatellitePassConstraints', 'getNextPass', 'getPassList', 'toTopocentric',
           'fromTopocentric', 'getAltitude', 'getAzimuth', 'azimuthAngleString')


def _azimuth_angle_string(azimuth):
    # azimuth in degrees
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


def _check_real_type(param, paramName: str):
    if not isinstance(param, (int, float)):
        raise TypeError(f'{paramName} parameter must be a real number')
    return param


class PositionInfo:

    __slots__ = '_altitude', '_azimuth', '_direction', '_time', '_illuminated', '_unobscured', '_visible'

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, illuminated: bool = False,
                 unobscured: bool = False):
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')

        self._altitude = _check_real_type(altitude, 'altitude')
        self._azimuth = _check_real_type(azimuth, 'azimuth')
        self._direction = _azimuth_angle_string(azimuth)
        self._time = time
        self._illuminated = index(illuminated)
        self._unobscured = index(unobscured)
        self._visible = self._illuminated and self._unobscured

    def __iter__(self):
        yield from {
            "altitude": self._altitude,
            "azimuth": self._azimuth,
            "direction": self._direction,
            "time": dict(self._time),
            "illuminated": self._illuminated,
            "unobscured": self._unobscured,
            "visible": self._visible
        }.items()

    def __str__(self):
        return ' {{:^17}} | {:^12} | {:^{altW}.{p}f} | {:^{azw}.{p}f} {:^5} | {!s:^11} | {!s:^10} | {!s:^7} ' \
            .format(self._time.time(), self._altitude, self._azimuth, '({})'.format(self._direction),
                    bool(self._illuminated), bool(self._unobscured), bool(self._visible), altW=8, azw=6, p=2)

    def __repr__(self):
        return json.dumps(self, default=lambda o: dict(o))

    @property
    def altitude(self):
        return self._altitude

    @property
    def azimuth(self):
        return self._azimuth

    @property
    def direction(self):
        return self._direction

    @property
    def time(self):
        return self._time

    @property
    def illuminated(self):
        return self._illuminated

    @property
    def unobscured(self):
        return self._unobscured

    @property
    def visible(self):
        return self._visible


def _check_info_type(param, paramName: str):
    if not isinstance(param, PositionInfo):
        raise TypeError(f'{paramName} parameter must be PositionInfo type')
    return param


def _check_info_type_none(param, paramName: str):
    if param is not None:
        return _check_info_type(param, paramName)
    return param


class SatellitePass:

    __slots__ = '_riseInfo', '_setInfo', '_maxInfo', '_firstUnobscured', '_lastUnobscured', '_firstIlluminated', \
                '_lastIlluminated', '_illuminated', '_unobscured', '_visible', '_name'

    def __init__(self, riseInfo: PositionInfo, setInfo: PositionInfo, maxInfo: PositionInfo, *,
                 firstUnobscuredInfo: PositionInfo = None, lastUnobscuredInfo: PositionInfo = None,
                 firstIlluminatedInfo: PositionInfo = None, lastIlluminatedInfo: PositionInfo = None, name: str = ''):
        self._riseInfo = _check_info_type(riseInfo, 'riseInfo')
        self._setInfo = _check_info_type(setInfo, 'setInfo')
        self._maxInfo = _check_info_type(maxInfo, 'maxInfo')
        self._firstUnobscured = _check_info_type_none(firstUnobscuredInfo, 'firstUnobscuredInfo')
        self._lastUnobscured = _check_info_type_none(lastUnobscuredInfo, 'lastUnobscuredInfo')
        self._firstIlluminated = _check_info_type_none(firstIlluminatedInfo, 'firstIlluminatedInfo')
        self._lastIlluminated = _check_info_type_none(lastIlluminatedInfo, 'lastIlluminatedInfo')
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

    def __iter__(self):
        yield from {
            'riseInfo': dict(self._riseInfo),
            'setInfo': dict(self._setInfo),
            'maxInfo': dict(self._maxInfo),
            'firstUnobscuredInfo': dict(self._firstUnobscured) if self._firstUnobscured is not None else None,
            'lastUnobscuredInfo': dict(self._lastUnobscured) if self._lastUnobscured is not None else None,
            'firstIlluminatedInfo': dict(self._firstIlluminated) if self._firstIlluminated is not None else None,
            'lastIlluminatedInfo': dict(self._lastIlluminated) if self._lastIlluminated is not None else None,
            'illuminated': self._illuminated,
            'unobscured': self._unobscured,
            'visible': self._visible
        }.items()

    def _get_infos(self):
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
        width = 97
        if self._name == '':
            title = '{:^{w}}'.format('Pass details at {}'.format(self._maxInfo.time.date()), w=width)
        else:
            title = '{:^{w}}'.format('Pass details for {}, at {}'.format(self._name, self._maxInfo.time.date()),
                                     w=width)
        heading = ' {:^17} | {:^12} | altitude | {:^12} | illuminated | unobscured | visible '\
            .format('instance', 'time', 'azimuth')
        lineBreak = '-' * width
        string = title + '\n' + heading + '\n'

        infos = [i for i in self._get_infos() if i[1] is not None]
        infos.sort(key=lambda o: o[1].time.value)

        for name, info in infos:
            string += '{}\n{}\n'.format(lineBreak, str(info).format(name))
        return string

    def __repr__(self):
        return json.dumps(self, default=lambda o: dict(o))

    @property
    def riseInfo(self):
        return self._riseInfo

    @property
    def setInfo(self):
        return self._setInfo

    @property
    def maxInfo(self):
        return self._maxInfo

    @property
    def firstUnobscuredInfo(self):
        return self._firstUnobscured

    @property
    def lastUnobscuredInfo(self):
        return self._lastUnobscured

    @property
    def firstIlluminatedInfo(self):
        return self._firstIlluminated

    @property
    def lastIlluminatedInfo(self):
        return self._lastIlluminated

    @property
    def unobscured(self):
        return self._unobscured

    @property
    def illuminated(self):
        return self._illuminated

    @property
    def visible(self):
        return self._visible


def _check_altitudes(minAlt, maxAlt):
    if minAlt is not None and maxAlt is not None and minAlt > maxAlt:
        raise PassConstraintException('Cannot set a minimum altitude constraint less than a maximum altitude '
                                      'constraint.')
    if minAlt is not None and minAlt > 90:
        raise ValueError('minAltitude cannot be above 90 degrees')
    if maxAlt is not None and maxAlt < 0:
        raise ValueError('maxAltitude cannot be below 0 degrees')


def _check_duration(minDur, maxDur):
    if minDur is not None and maxDur is not None and minDur > maxDur:
        raise PassConstraintException('Cannot set a minimum duration constraint longer than a maximum duration'
                                      'constraint.')
    if maxDur is not None and maxDur < 0:
        raise ValueError('maxDuration cannot be less than 0')


class SatellitePassConstraints:

    __slots__ = '_minAltitude', '_maxAltitude', '_minDuration', '_maxDuration', '_illuminated', '_unobscured', \
                '_visible'

    def __init__(self, *, minAltitude: float = None, maxAltitude: float = None, minDuration: float = None,
                 maxDuration: float = None, illuminated: bool = None, unobscured: bool = None, visible: bool = None):
        _check_real_type(minAltitude, 'minAltitude')
        _check_real_type(maxAltitude, 'maxAltitude')
        _check_altitudes(minAltitude, maxAltitude)
        _check_duration(minDuration, maxDuration)
        self._minAltitude = minAltitude
        self._maxAltitude = maxAltitude
        self._minDuration = minDuration
        self._maxDuration = maxDuration
        self._illuminated = index(illuminated)
        self._unobscured = index(unobscured)
        self._visible = index(visible)

    @property
    def minAltitude(self):
        return self._minAltitude

    @minAltitude.setter
    def minAltitude(self, value):
        _check_real_type(value, 'value')
        _check_altitudes(value, self._maxAltitude)
        self._minAltitude = value

    @minAltitude.deleter
    def minAltitude(self):
        self._minAltitude = None

    @property
    def maxAltitude(self):
        return self._maxAltitude

    @maxAltitude.setter
    def maxAltitude(self, value):
        _check_real_type(value, 'value')
        _check_altitudes(self._minAltitude, value)
        self._maxAltitude = value

    @maxAltitude.deleter
    def maxAltitude(self):
        self._maxAltitude = None

    @property
    def minDuration(self):
        return self._minDuration

    @minDuration.setter
    def minDuration(self, value):
        _check_real_type(value, 'value')
        _check_duration(value, self._maxDuration)
        self._minDuration = value

    @minDuration.deleter
    def minDuration(self):
        self._minDuration = None

    @property
    def maxDuration(self):
        return self._maxDuration

    @maxDuration.setter
    def maxDuration(self, value):
        _check_real_type(value, 'value')
        _check_duration(self._minDuration, value)
        self._maxDuration = value

    @maxDuration.deleter
    def maxDuration(self):
        self._maxDuration = None

    @property
    def illuminated(self):
        return self._illuminated

    @illuminated.setter
    def illuminated(self, value):
        self._illuminated = index(value)

    @illuminated.deleter
    def illuminated(self):
        self._illuminated = None

    @property
    def unobscured(self):
        return self._unobscured

    @unobscured.setter
    def unobscured(self, value):
        self._unobscured = index(value)

    @unobscured.deleter
    def unobscured(self):
        self._unobscured = None

    @property
    def visible(self):
        return self._visible

    @visible.setter
    def visible(self, value):
        self._visible = index(value)

    @visible.deleter
    def visible(self):
        self._unobscured = None


def __max_pass_anomalies(position, velocity, mu, pVector, ecc):
    # eccentricVector = _compute_eccentric_vector(position, velocity, mu)
    eccentricVector = computeEccentricVector(position, velocity, mu)
    trueAnomaly_n = _true_anomaly_from_state(position, velocity, mu)
    # todo: use eccentricVector.mag() if we know it's accurate
    meanAnomaly_n = _true_to_mean_anomaly(trueAnomaly_n, ecc)
    trueAnomaly_n1 = vang(eccentricVector, pVector)
    if norm(cross(eccentricVector, pVector)) != norm(cross(position, velocity)):
        trueAnomaly_n1 = TWOPI - trueAnomaly_n1
    meanAnomaly_n1 = _true_to_mean_anomaly(trueAnomaly_n1, ecc)
    return trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1


def _get_p_vector(geo: GeoPosition, position: Vector, velocity: Vector, time: JulianDate) -> Vector:
    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    elevation = geo.elevation

    # zenith vector for geo-position
    zeta = _compute_zenith_vector(latitude, longitude, elevation, time)
    # position vector for geo-position
    gamma = _compute_position_vector(latitude, longitude, elevation, time)
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


def _orbit_altitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    state = satellite.getState(time)
    elements = satellite.getElements(time)

    gamma = _compute_position_vector(radians(geo.latitude), radians(geo.longitude), geo.elevation, time)
    pVector = _get_p_vector(geo, *state, time)

    # eccentricVector = _compute_eccentric_vector(*state, satellite.body.mu)
    eccentricVector = computeEccentricVector(*state, satellite.body.mu)
    trueAnomaly = vang(eccentricVector, pVector)
    pSat = norm(pVector) * ((elements.sma * (1 - elements.ecc * elements.ecc)) / (1 + elements.ecc * cos(trueAnomaly)))

    angle = vang(pVector - gamma, pSat - gamma)
    return angle if pSat.mag2() > pVector.mag2() else -angle


def _time_to_plane(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    jd = time
    alt = _orbit_altitude(satellite, geo, time)
    if alt > 0:
        return jd

    while abs(alt) > 2.7e-4:  # one arc-second
        # todo: find an analytical guess of this dt
        dt = -alt / 360
        jd = jd.future(dt)
        alt = _orbit_altitude(satellite, geo, jd)
    return jd


def _next_pass_max_approx(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    if _orbit_altitude(satellite, geo, time) < 0:
        tn = _time_to_plane(satellite, geo, time)
    else:
        tn = time

    mu = satellite.body.mu
    elements = satellite.getElements(time)
    ecc = elements.ecc
    # mean motion in radians / day
    meanMotion = _sma_to_mean_motion(elements.sma, mu) * 86400

    state = satellite.getState(tn)
    pVector = _get_p_vector(geo, *state, tn)
    _, _, meanAnomaly_n, meanAnomaly_n1 = __max_pass_anomalies(*state, mu, pVector, ecc)
    dma = meanAnomaly_n1 - meanAnomaly_n
    if meanAnomaly_n1 < meanAnomaly_n:
        dma += TWOPI
    tn = tn.future(dma / meanMotion)

    state = satellite.getState(tn)
    pVector = _get_p_vector(geo, *state, tn)
    while vang(state[0], pVector) > 4.85e-6:  # one arc-second
        trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1 = __max_pass_anomalies(*state, mu, pVector, ecc)
        dma = meanAnomaly_n1 - meanAnomaly_n
        if trueAnomaly_n1 <= trueAnomaly_n:
            if (trueAnomaly_n - trueAnomaly_n1) >= pi:
                dma += TWOPI
        else:
            if (trueAnomaly_n1 - trueAnomaly_n) > pi:
                dma -= TWOPI
        tn = tn.future(dma / meanMotion)
        print("dma: ", dma, "tn: ", tn)
        state = satellite.getState(tn)
        pVector = _get_p_vector(geo, *state, tn)

    return tn

def _vector_almost_equal(vec1, vec2, epsilon=1e-5):
    return dot(vec1, vec2) > (1 - epsilon)

    # for v1, v2 in zip(vec1, vec2):
    #     if abs(v1 - v2) > epsilon:
    #         return False
    #
    # return True


def _next_pass_max_approx_fix(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    '''this line is to count towards the doc'''
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

    if _orbit_altitude(satellite, geo, time) < 0:
        tn = _time_to_plane(satellite, geo, time)
    else:
        tn = time

    mu = satellite.body.mu
    state = satellite.getState(tn)
    # elements = elementsFromState(*state)
    elements = Elements.fromState(*state, tn)
    # mean motion in radians / day
    # meanMotion = _sma_to_mean_motion(elements[4], mu) * 86400
    meanMotion = _sma_to_mean_motion(elements.sma, mu) * 86400

    pVector = _get_p_vector(geo, *state, tn)
    eccVector = computeEccentricVector(*state, mu)
    angularMomentumNorm = norm(cross(*state))

    #todo: replace the methods when the vector equality works
    taPVector = vang(eccVector, pVector)
    if not _vector_almost_equal(norm(cross(eccVector, pVector)), angularMomentumNorm):
        taPVector = TWOPI - taPVector

    # todo: replace this with values from elements variable
    # maPVector = _true_to_mean_anomaly(taPVector, elements[3])
    maPVector = _true_to_mean_anomaly(taPVector, elements.ecc)
    # maPosition = elements[-2]
    maPosition = elements.meanAnomaly

    dma = maPVector - maPosition
    if dma < 0:
        dma += TWOPI

    tn = tn.future(dma / meanMotion)
    state = satellite.getState(tn)
    pVector = _get_p_vector(geo, *state, tn)

    while vang(state[0], pVector) > 4.58e-6:    # 1 arc second
        # elements = elementsFromState(*state)
        elements = Elements.fromState(*state, tn)
        # mean motion in radians / day
        # meanMotion = _sma_to_mean_motion(elements[4], mu) * 86400
        meanMotion = _sma_to_mean_motion(elements.sma, mu) * 86400

        eccVector = computeEccentricVector(*state, mu)
        angularMomentumNorm = norm(cross(*state))

        taPVector = vang(eccVector, pVector)
        if not _vector_almost_equal(norm(cross(eccVector, pVector)), angularMomentumNorm):
            taPVector = TWOPI - taPVector

        # maPVector = _true_to_mean_anomaly(taPVector, elements[3])
        maPVector = _true_to_mean_anomaly(taPVector, elements.ecc)
        # maPosition = elements[-2]
        maPosition = elements.meanAnomaly

        dma = maPVector - maPosition
        if dma < 0:
            dma += TWOPI

        tn = tn.future(dma / meanMotion)
        state = satellite.getState(tn)
        pVector = _get_p_vector(geo, *state, tn)

    return tn


def _compute_altitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    position = satellite.getState(time)[0]
    topocentricPosition = _to_topocentric(position, geo, time)
    return asin(topocentricPosition[2] / topocentricPosition.mag())


def _max_pass_refine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    # todo: utilize the dadt values to future time increase and make this a more analytical guess
    oneSecond = 1 / 86400
    altitude = _compute_altitude(satellite, geo, time)
    futureAltitude = _compute_altitude(satellite, geo, time.future(oneSecond))
    pastAltitude = _compute_altitude(satellite, geo, time.future(-oneSecond))

    parity = 0  # to please the editor
    if altitude >= futureAltitude and altitude >= pastAltitude:
        return time
    elif altitude < futureAltitude:
        parity = 1
    elif altitude < pastAltitude:
        parity = -1

    jd = time
    nextAltitude = _compute_altitude(satellite, geo, time.future(parity * oneSecond))
    while altitude < nextAltitude:
        jd = jd.future(parity * oneSecond)
        altitude = _compute_altitude(satellite, geo, jd)
        nextAltitude = _compute_altitude(satellite, geo, jd.future(parity * oneSecond))

    return jd


def _next_pass_max(satellite: Orbitable, geo: GeoPosition, time: JulianDate, timeout: float) -> JulianDate:
    nextMax = _next_pass_max_approx_fix(satellite, geo, time)
    altitude = _orbit_altitude(satellite, geo, nextMax)
    while altitude < 0:
        if nextMax - time < timeout:
            nextMax = _next_pass_max_approx_fix(satellite, geo, nextMax.future(0.001))
            altitude = _orbit_altitude(satellite, geo, nextMax)
        else:
            # todo: make a more specific exception type for this
            raise NoPassException('timed out looking for next pass')
    return _max_pass_refine(satellite, geo, nextMax)


def _rise_set_times_approx(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    elements = satellite.getElements(time)
    mu = satellite.body.mu
    ecc = elements.ecc
    a = elements.sma
    c = a * ecc
    b = sqrt(a * a - c * c)

    state = satellite.getState(time)
    hNorm = norm(cross(*state))
    # u = norm(_compute_eccentric_vector(*state, mu)) * a
    u = norm(computeEccentricVector(*state, mu)) * a
    v = norm(cross(hNorm, u)) * b
    ce = -norm(u) * c

    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    elevation = geo.elevation
    zeta = _compute_zenith_vector(latitude, longitude, elevation, time)
    gamma = _compute_position_vector(latitude, longitude, elevation, time)

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
    trueAnomaly_10 = _true_anomaly_from_state(*state, mu)
    trueAnomaly_20 = _true_anomaly_from_state(*state, mu)
    # meanMotion in rev / day
    meanMotion = _sma_to_mean_motion(a, mu) * 86400 / TWOPI
    jd1 = _nearest_true_anomaly(meanMotion, ecc, trueAnomaly_10, time, trueAnomaly_1)
    jd2 = _nearest_true_anomaly(meanMotion, ecc, trueAnomaly_20, time, trueAnomaly_2)
    if jd1 < jd2:
        return jd1, jd2
    return jd2, jd1


'''this method gives us some issues. it is not uncommon for an infinite loop to occur
while approaching the horizon. the altitude oscillates around 0, never converging. A
brute force fix could be to take a proportion of the approximated time, but would 
converge much more slowly. Current idea is to interpolate the 'right' dt term after
actually computing the future value before looping to the next iteration.'''
# todo: can we account for earth's rotation here? will that be more accurate/fix our issues?

def _horizon_time_refine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    # needed to add a proportional governor because of cases oscillating around the zero altitude
    # todo: look into a pid loop for this
    p = 1.0
    iterCount = 0

    state = satellite.getState(time)
    topocentricPosition = _to_topocentric(state[0], geo, time)
    altitude = asin(topocentricPosition[2] / topocentricPosition.mag())
    while abs(altitude) > 4.85e-6:  # one arc-second
        iterCount += 1
        if iterCount % 10 == 0:
            p /= 2

        topocentricVelocity = rotateEulerTo(
            ZYX,
            Angles(
                radians(geo.longitude) + earthOffsetAngle(time),
                radians(90 - geo.latitude),
                0.0
            ),
            state[1]
        )

        dz = topocentricPosition.mag() * altitude
        dt = (-dz / topocentricVelocity[2]) / 86400
        time = time.future(dt * p)

        state = satellite.getState(time)
        topocentricPosition = _to_topocentric(state[0], geo, time)
        altitude = asin(topocentricPosition[2] / topocentricPosition.mag())
    return time


def _rise_set_times(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    riseTime, setTime = _rise_set_times_approx(satellite, geo, time)
    riseTime = _horizon_time_refine(satellite, geo, riseTime)
    setTime = _horizon_time_refine(satellite, geo, setTime)
    return riseTime, setTime


def _get_special_times(riseTime: JulianDate, setTime: JulianDate, beginTime: JulianDate, finishTime: JulianDate)\
        -> (JulianDate, JulianDate, bool, bool):
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


def _derive_basic_info(satellite: Orbitable, geo: GeoPosition, time: JulianDate, isIlluminated: bool,
                       isUnobscured: bool) -> PositionInfo:
    position = satellite.getState(time)[0]
    topocentricPosition = _to_topocentric(position, geo, time)
    altitude = degrees(asin(topocentricPosition[2] / topocentricPosition.mag()))
    azimuth = degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))
    return PositionInfo(altitude, azimuth, time, isIlluminated, isUnobscured)


def getNextPass(satellite: Orbitable, geo: GeoPosition, timeOrPass: JulianDate | SatellitePass,
                timeout: float = 7) -> SatellitePass:
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
    nextPassTime = _next_pass_max(satellite, geo, time, timeout)

    riseTime, setTime = _rise_set_times(satellite, geo, nextPassTime)
    enterTime, exitTime = getShadowTimes(satellite, nextPassTime, Shadow.PENUMBRA)
    # todo: if timezone of riseTime is different than the timezone of geo, then the rise/set times may be for the
    #   wrong day and the unobscured values may be wrong.
    # todo: idea: check sun altitude and always find the next sunset based on the angle, then find the previous rise,
    #   then the above issue shouldn't exist
    sunRiseTime, sunSetTime = getSunTimes(riseTime, geo)

    firstIlluminatedTime, lastIlluminatedTime, riseIlluminated, setIlluminated \
        = _get_special_times(riseTime, setTime, enterTime, exitTime)
    firstUnobscuredTime, lastUnobscuredTime, riseUnobscured, setUnobscured \
        = _get_special_times(riseTime, setTime, sunRiseTime, sunSetTime)
    maxIlluminated = not (enterTime <= nextPassTime <= exitTime)
    maxUnobscured = (nextPassTime < sunRiseTime) or (nextPassTime >= sunSetTime)

    riseInfo = _derive_basic_info(satellite, geo, riseTime, riseIlluminated, riseUnobscured)
    setInfo = _derive_basic_info(satellite, geo, setTime, setIlluminated, setUnobscured)
    maxInfo = _derive_basic_info(satellite, geo, nextPassTime, maxIlluminated, maxUnobscured)

    firstIlluminatedInfo = lastIlluminatedInfo = firstUnobscuredInfo = lastUnobscuredInfo = None
    if firstIlluminatedTime is not None and firstIlluminatedTime != riseTime:
        firstIlluminatedInfo = _derive_basic_info(satellite, geo, firstIlluminatedTime, True,
                                                  firstIlluminatedTime < sunRiseTime
                                                  or firstIlluminatedTime >= sunSetTime)
    if lastIlluminatedTime is not None and lastIlluminatedTime != setTime:
        lastIlluminatedInfo = _derive_basic_info(satellite, geo, lastIlluminatedTime, True,
                                                 lastIlluminatedTime < sunRiseTime
                                                 or lastIlluminatedTime >= sunSetTime)
    if firstUnobscuredTime is not None and firstUnobscuredTime == sunSetTime:
        firstUnobscuredInfo = _derive_basic_info(satellite, geo, firstUnobscuredTime,
                                                 firstUnobscuredTime < enterTime
                                                 or firstUnobscuredTime >= exitTime, True)
    if lastUnobscuredTime is not None and lastUnobscuredTime == sunRiseTime:
        lastUnobscuredInfo = _derive_basic_info(satellite, geo, lastUnobscuredTime,
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
    if not isinstance(vector, Vector):
        raise TypeError('vector parameter must be EVector type')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    return _to_topocentric(vector, geo, time)


def fromTopocentric(vector: Vector, time: JulianDate, geo: GeoPosition) -> Vector:
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
    # matrix = getEulerMatrix(
    #     ZYX,
    #     EulerAngles(
    #         radians(geo.longitude) + earthOffsetAngle(time),
    #         radians(90 - geo.latitude),
    #         0.0
    #     )
    # )

    # return undoRotateToThenOffset(matrix, geoVector, vector)
    return rotateOffsetFrom(matrix, geoVector, vector)


def getAltitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    position = satellite.getState(time)[0]
    topocentricPosition = _to_topocentric(position, geo, time)
    arg = topocentricPosition[2] / topocentricPosition.mag()
    return degrees(asin(arg))


def getAzimuth(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    if not isinstance(satellite, Orbitable):
        raise TypeError('satellite parameter must be Orbitable subtype')
    if not isinstance(geo, GeoPosition):
        raise TypeError('geo parameter must be GeoPosition type')
    if not isinstance(time, JulianDate):
        raise TypeError('time parameter must be JulianDate type')

    position = satellite.getState(time)[0]
    topocentricPosition = _to_topocentric(position, geo, time)
    return degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))


def azimuthAngleString(azimuth: float) -> str:
    # azimuth in degrees
    _check_real_type(azimuth, 'azimuth')
    if not (0 <= azimuth <= 360):
        raise ValueError('azimuth angle must be between 0 and 360 degrees')
    return _azimuth_angle_string(azimuth)


# todo: create alt-az type and methods for it
