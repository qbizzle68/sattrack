import json
import re
from math import cos, radians, pi, sqrt, acos, sin, degrees, asin, atan2
from operator import index as _index

from pyevspace import EVector, cross, dot, norm, vang
from sattrack.eclipse import getShadowTimes, _PENUMBRA, Shadow

from sattrack.exceptions import NoPassException, TLEException, PositiveZeroException, PassConstraintException
from sattrack.rotation.order import Axis, ZYX, ZXZ
from sattrack.rotation.rotation import getMatrix, rotateOrderTo, EulerAngles
from sattrack.sampa import getSunTimes
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.structures._coordinates import _compute_zenith_vector, _compute_position_vector
from sattrack.structures._orbit import _compute_eccentric_vector, _true_anomaly_from_state, _true_to_mean_anomaly, \
    _sma_to_mean_motion, _nearest_true_anomaly
from sattrack.structures.coordinates import GeoPosition
from sattrack.structures.elements import computeEccentricVector, raanProcessionRate, OrbitalElements
from sattrack.structures.orbit import Orbitable
# from sattrack.structures.satellite import Satellite
from sattrack.structures.tle import TwoLineElement
from sattrack.sun import TwilightType, getTwilightType, getSunPosition, getSunRiseSetTimes2
from sattrack.topos import getPVector, getAltitude, azimuthAngleString
from sattrack.topocentric import toTopocentric
from sattrack.util.anomalies import trueToMean, timeToNearestTrueAnomaly, computeTrueAnomaly
from sattrack.util.constants import SUN_RADIUS, TWOPI, EARTH_FLATTENING, EARTH_EQUITORIAL_RADIUS
from sattrack.util.conversions import atan3


class PositionInfo:

    __slots__ = '_altitude', '_azimuth', '_direction', '_time', '_illuminated', '_unobscured', '_visible'

    def __init__(self, altitude: float, azimuth: float, time: JulianDate, illuminated: bool = False,
                 unobscured: bool = False):
        self._altitude = altitude
        self._azimuth = azimuth
        self._direction = azimuthAngleString(azimuth)
        self._time = time
        self._illuminated = illuminated
        self._unobscured = unobscured
        self._visible = illuminated and unobscured

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
            .format(self._time.time(), self._altitude, self._azimuth, '({})'.format(self._direction), self._illuminated,
                    self._unobscured, self._visible, altW=8, azw=6, p=2)


        # return f'Altitude: {"%.2f" % self._altitude}\nAzimuth: {"%.2f" % self._azimuth} ({self._direction})\nTime: ' \
        #        f'{self._time}\nIlluminated: {self._illuminated}\nUnobscured: {self._unobscured}\nVisibility: ' \
        #        f'{self._visible}'

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


class SatellitePass:
    def __init__(self, riseInfo: PositionInfo, setInfo: PositionInfo, maxInfo: PositionInfo, *,
                 firstUnobscuredInfo: PositionInfo = None, lastUnobscuredInfo: PositionInfo = None,
                 firstIlluminatedInfo: PositionInfo = None, lastIlluminatedInfo: PositionInfo = None):
        self._riseInfo = riseInfo
        self._setInfo = setInfo
        self._maxInfo = maxInfo
        self._firstUnobscured = firstUnobscuredInfo
        self._lastUnobscured = lastUnobscuredInfo
        self._firstIlluminated = firstIlluminatedInfo
        self._lastIlluminated = lastIlluminatedInfo
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
        #todo: actually add this
        self._name = 'iss'

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
            # title = f'Pass details at {self._maxInfo.time.date()}'
            # title = ':^{width}'.format('Pass details at {}', self._maxInfo.time.date(), width=width)
            title = '{:^{width}}'.format('Pass details at {}'.format(self._maxInfo.time.date()), width=width)
        else:
            # title = f'Pass details for {self._name} at {self._maxInfo.time.date()}'
            # title = ':^{width}'.format('Pass details for {}, at {}', self._name, self._maxInfo.time.date(), width=width)
            title = '{:^{w}}'.format('Pass details for {}, at {}'.format(self._name, self._maxInfo.time.date()), w=width)
        heading = ' {:^17} | {:^12} | altitude | {:^12} | illuminated | unobscured | visible ' \
            .format('instance', 'time', 'azimuth')
        lineBreak = '-' * width
        string = title + '\n' + heading + '\n'  # + lineBreak

        infos = [i for i in self._get_infos() if i[1] is not None]
        infos.sort(key=lambda o: o[1].time.value)

        # infos = [(name, info) for name, info in dict(self).items() if
        #          info is not None and not isinstance(info, bool)]
        # print(infos)
        # infos.sort(key=lambda o: o[1]['time']['dayNumber'] + o[1]['time']['dayFraction'])
        # print(infos)

        for name, info in infos:
            string += '{}\n{}\n'.format(lineBreak, str(info).format(name))
        return string


        # # prettify this using new formatting learned
        # riseStr = f'Rise:\n{self._riseInfo}'
        # setStr = f'Set:\n{self._setInfo}'
        # maxStr = f'Max:\n{self._maxInfo}'
        # strDict = {info.time.value: string for info, string in
        #            zip((self._riseInfo, self._setInfo, self._maxInfo), (riseStr, setStr, maxStr))}
        # if self._firstIlluminated is not None:
        #     firstIllStr = f'First Illuminated:\n{self._firstIlluminated}'
        #     strDict[self._firstIlluminated.time.value] = firstIllStr
        # if self._lastIlluminated is not None:
        #     lastIllStr = f'Last Illuminated:\n{self._lastIlluminated}'
        #     strDict[self._lastIlluminated.time.value] = lastIllStr
        # if self._firstUnobscured is not None:
        #     firstStr = f'First Unobscured:\n{self._firstUnobscured}'
        #     strDict[self._firstUnobscured.time.value] = firstStr
        # if self._lastUnobscured is not None:
        #     lastStr = f'Last Unobscured:\n{self._lastUnobscured}'
        #     strDict[self._lastUnobscured.time.value] = lastStr
        # rtn = ""
        # for i in range(len(strDict)):
        #     minTm = min(strDict.keys())
        #     rtn += strDict[minTm] + '\n'
        #     strDict.pop(minTm)
        # rtn.rstrip()
        # rtn += f'---------------------------------\nIlluminated: {self._illuminated}\nUnobscured: {self._unobscured}' \
        #        f'\nVisible: {self._visible}'
        # return rtn

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
    def firstIlluminated(self):
        return self._firstIlluminated

    @property
    def lastIlluminated(self):
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
        _check_altitudes(minAltitude, maxAltitude)
        _check_duration(minDuration, maxDuration)
        self._minAltitude = minAltitude
        self._maxAltitude = maxAltitude
        self._minAltitude = minAltitude
        self._maxAltitude = maxAltitude
        self._illuminated = illuminated
        self._unobscured = unobscured
        self._visible = visible

    @property
    def minAltitude(self):
        return self._minAltitude

    @minAltitude.setter
    def minAltitude(self, value):
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
        self._illuminated = _index(value)

    @illuminated.deleter
    def illuminated(self):
        self._illuminated = None

    @property
    def unobscured(self):
        return self._unobscured

    @unobscured.setter
    def unobscured(self, value):
        self._unobscured = _index(value)

    @unobscured.deleter
    def unobscured(self):
        self._unobscured = None

    @property
    def visible(self):
        return self._visible

    @visible.setter
    def visible(self, value):
        self._visible = _index(value)

    @visible.deleter
    def visible(self):
        self._unobscured = None


def _orbit_altitude(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> float:
    state = satellite.getState(time)
    elements = satellite.getElements(time)

    gamma = _compute_position_vector(radians(geo.latitude), radians(geo.longitude), geo.elevation, time)
    pVector = _get_p_vector(geo, *state, time)

    # find the true anomaly of this vector if it were a position vector
    eccentricVector = _compute_eccentric_vector(*state, satellite.body.mu)
    trueAnomaly = vang(eccentricVector, pVector)
    pSat = norm(pVector) * ((elements.sma * (1 - elements.ecc * elements.ecc)) / (1 + elements.ecc * cos(trueAnomaly)))

    angle = degrees(vang(pVector - gamma, pSat - gamma))
    return angle if pSat.mag2() > pVector.mag2() else - angle


def _time_to_plane(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    jd = time
    alt = _orbit_altitude(satellite, geo, time)
    if alt > 0:
        return jd

    while abs(alt) > 2.7e-4:  # one arc-second on either side
        # todo: find an analytical guess of this dt
        dt = -alt / 360
        jd = jd.future(dt)
        alt = _orbit_altitude(satellite, geo, jd)
    return jd


def _get_p_vector(geo: GeoPosition, position: EVector, velocity: EVector, time: EVector) -> EVector:
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
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)

    # compute exact solution for parameterized vector which yields the nearest point to intersection
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    return r + v * t


def __max_pass_anomalies(position, velocity, geo, time, mu, pVector, ecc):
    # pVector = _get_p_vector(geo, position, velocity, time)
    eccentricVector = _compute_eccentric_vector(position, velocity, mu)
    trueAnomaly_n = _true_anomaly_from_state(position, velocity, mu)
    # todo: use eccentricVector.mag() if we know it's accurate
    meanAnomaly_n = _true_to_mean_anomaly(trueAnomaly_n, ecc)
    trueAnomaly_n1 = vang(eccentricVector, pVector)
    if norm(cross(eccentricVector, pVector)) != norm(cross(position, velocity)):
        trueAnomaly_n1 = TWOPI - trueAnomaly_n1
    meanAnomaly_n1 = _true_to_mean_anomaly(trueAnomaly_n1, ecc)
    return trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1


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
    _, _, meanAnomaly_n, meanAnomaly_n1 = __max_pass_anomalies(*state, geo, tn, mu, pVector, ecc)
    dma = meanAnomaly_n1 - meanAnomaly_n
    if meanAnomaly_n1 < meanAnomaly_n:
        dma += TWOPI
    tn = tn.future(dma / meanMotion)

    state = satellite.getState(tn)
    pVector = _get_p_vector(geo, *state, tn)
    while vang(state[0], pVector) > 4.85e-6:  # one arc-second
        trueAnomaly_n, trueAnomaly_n1, meanAnomaly_n, meanAnomaly_n1 = __max_pass_anomalies(*state, geo, tn, mu, pVector, ecc)
        dma = meanAnomaly_n1 - meanAnomaly_n
        if trueAnomaly_n1 <= trueAnomaly_n:
            if (trueAnomaly_n - trueAnomaly_n1) >= pi:
                dma += TWOPI
        else:
            if (trueAnomaly_n1 - trueAnomaly_n) > pi:
                dma -= TWOPI
        tn = tn.future(dma / meanMotion)
        state = satellite.getState(tn)
        pVector = _get_p_vector(geo, *state, tn)

    return tn

    # # todo: come back and clean this up
    # if _orbit_altitude(satellite, geo, time) < 0:
    #     t0 = _time_to_plane(satellite, geo, time)
    # else:
    #     t0 = time
    #
    # MU = satellite.body.mu
    # state = satellite.getState(t0)
    # pVector = _get_p_vector(geo, *state, t0)
    # trueAnomaly_0 = _true_anomaly_from_state(*state, MU)
    # eccentricity = satellite.getElements(time).ecc
    # meanAnomaly_0 = _true_to_mean_anomaly(trueAnomaly_0, eccentricity)
    # eccentricVector = _compute_eccentric_vector(*state, MU)
    # trueAnomaly_1 = vang(eccentricVector, pVector)
    # if norm(cross(eccentricVector, pVector)) != norm(cross(*state)):
    #     trueAnomaly_1 = TWOPI - trueAnomaly_1
    # meanAnomaly_1 = _true_to_mean_anomaly(trueAnomaly_1, eccentricity)
    # if meanAnomaly_1 < meanAnomaly_0:
    #     dma0 = meanAnomaly_1 - meanAnomaly_0 + TWOPI
    # else:
    #     dma0 = meanAnomaly_1 - meanAnomaly_0
    # meanMotion = _sma_to_mean_motion(satellite.getElements(time).sma, satellite.body.mu) * 86400
    # tn = t0.future(dma0 / meanMotion)
    #
    # # iterate towards answer moving forward or backward
    # state = satellite.getState(tn)
    # pVector = _get_p_vector(geo, *state, tn)
    # while vang(state[0], pVector) > 4.85e-6: # 1 arc second
    #     pVector = _get_p_vector(geo, *state, tn)
    #     eccentricVector = _compute_eccentric_vector(*state, MU)
    #     trueAnomaly_n = _true_anomaly_from_state(*state, MU)
    #     meanAnomaly_n = _true_to_mean_anomaly(trueAnomaly_n, eccentricity)
    #     trueAnomaly_n1 = vang(eccentricVector, pVector)
    #     if norm(cross(eccentricVector, pVector)) != norm(cross(*state)):
    #         trueAnomaly_n1 = TWOPI - trueAnomaly_n1
    #     meanAnomaly_n1 = _true_to_mean_anomaly(trueAnomaly_n1, eccentricity)
    #     if trueAnomaly_n1 <= trueAnomaly_n:
    #         if (trueAnomaly_n - trueAnomaly_n1) < pi:
    #             dma = meanAnomaly_n1 - meanAnomaly_n
    #         else:
    #             dma = meanAnomaly_n1 + TWOPI - meanAnomaly_n
    #     else:
    #         if (trueAnomaly_n1 - trueAnomaly_n) > pi:
    #             dma = meanAnomaly_n1 - TWOPI - meanAnomaly_n
    #         else:
    #             dma = meanAnomaly_n1 - meanAnomaly_n
    #     tn = tn.future(dma / meanMotion)
    #     state = satellite.getState(tn)
    #     pVector = getPVector(geo, *state, tn)
    # # if _orbit_altitude(satellite, geo, tn) < 0:
    # #     return _next_pass_max_approx(satellite, geo, tn.future(0.001))
    # return tn


def _max_pass_refine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    # todo: utilize the dadt values to future time increase and make this a more analytical guess
    oneSecond = 1 / 86400
    # todo: see if we put this in a helper method
    alt = getAltitude(satellite, time, geo)
    futureAlt = getAltitude(satellite, time.future(oneSecond), geo)
    pastAlt = getAltitude(satellite, time.future(-oneSecond), geo)

    parity = 0  # to please the editor
    if alt >= futureAlt and alt >= pastAlt:
        return time
    elif alt < futureAlt:
        parity = 1
    elif alt < pastAlt:
        parity = -1

    jd = time
    nextAlt = getAltitude(satellite, time.future(parity * oneSecond), geo)
    while alt < nextAlt:
        jd = jd.future(parity * oneSecond)
        alt = getAltitude(satellite, jd, geo)
        nextAlt = getAltitude(satellite, jd.future(parity * oneSecond), geo)

    return jd


def _next_pass_max(satellite: Orbitable, geo: GeoPosition, time: JulianDate, timeout: float) -> JulianDate:
    nextMax = _next_pass_max_approx(satellite, geo, time)
    # todo: see if we put this in a helper method
    altitude = _orbit_altitude(satellite, geo, nextMax)
    while altitude < 0:
        if nextMax - time < timeout:
            nextMax = _next_pass_max_approx(satellite, geo, nextMax.future(0.001))
            altitude = _orbit_altitude(satellite, geo, nextMax)
        else:
            # todo: make a more specific exception type for this
            raise NoPassException('timed out looking for next pass')
    return _max_pass_refine(satellite, geo, nextMax)


def _rise_set_times_approx(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    elements = satellite.getElements(time)
    MU = satellite.body.mu
    a = elements.sma
    ecc = elements.ecc
    c = a * ecc
    b = sqrt(a * a - c * c)

    state = satellite.getState(time)
    hNorm = norm(cross(*state))
    u = norm(_compute_eccentric_vector(*state, MU)) * a
    v = norm(cross(hNorm, u)) * b
    ce = -norm(u) * c

    latitude = radians(geo.latitude)
    longitude = radians(geo.longitude)
    elevation = geo.elevation
    zeta = _compute_zenith_vector(latitude, longitude, elevation, time)
    gamma = _compute_position_vector(latitude, longitude, elevation, time)

    R = sqrt((dot(zeta, u) ** 2) + (dot(zeta, v) ** 2))
    try:
        beta = acos(dot(zeta, gamma - ce) / R)
    except ValueError:
        raise NoPassException(f'No pass during this time: {time}.')
    alpha = atan3(dot(zeta, v), dot(zeta, u))
    w1 = alpha + beta
    w2 = alpha - beta

    rho1 = (u * cos(w1) + v * sin(w1)).mag()
    trueAnomaly_1 = atan3(rho1 * sin(w1), rho1 * cos(w1) - a * ecc)
    rho2 = (u * cos(w2) + v * sin(w2)).mag()
    trueAnomaly_2 = atan3(rho2 * sin(w2), rho2 * cos(w2) - a * ecc)
    trueAnomaly_10 = _true_anomaly_from_state(*state, MU)
    trueAnomaly_20 = _true_anomaly_from_state(*state, MU)
    meanMotion = _sma_to_mean_motion(elements.sma, MU) * 86400 / TWOPI
    jd1 = _nearest_true_anomaly(meanMotion, ecc, trueAnomaly_10, time, trueAnomaly_1)
    jd2 = _nearest_true_anomaly(meanMotion, ecc, trueAnomaly_20, time, trueAnomaly_2)
    if jd1 < jd2:
        return jd1, jd2
    else:
        return jd2, jd1


def _horizon_time_refine(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> JulianDate:
    # needed to add a sort of proportional governor because of cases oscillating around the zero altitude
    p = 1.0
    iterCount = 0

    state = satellite.getState(time)
    # todo: look at replacing this
    topocentricPosition = toTopocentric(state[0], geo, time)
    altitude = asin(topocentricPosition[2] / topocentricPosition.mag())
    while abs(altitude) > 4.85e-6: # 1 arc-second
        iterCount += 1
        if iterCount % 10 == 0:
            p /= 2.0

        topocentricVelocity = rotateOrderTo(
            ZYX,
            EulerAngles(
                radians(geo.longitude) + earthOffsetAngle(time),
                radians(90 - geo.latitude),
                0.0),
            state[1])

        dz = topocentricPosition.mag() * altitude
        dt = (-dz / topocentricVelocity[2]) / 86400
        time = time.future(dt * p)  # p term is linear, so we can put it in the dz or dt term too

        state = satellite.getState(time)
        topocentricPosition = toTopocentric(state[0], geo, time)
        altitude = asin(topocentricPosition[2] / topocentricPosition.mag())
    return time


def _rise_set_times(satellite: Orbitable, geo: GeoPosition, time: JulianDate) -> (JulianDate, JulianDate):
    riseTime, setTime = _rise_set_times_approx(satellite, geo, time)
    riseTime = _horizon_time_refine(satellite, geo, riseTime)
    setTime = _horizon_time_refine(satellite, geo, setTime)
    return riseTime, setTime


# todo: better name
def _get_special_times(riseTime: JulianDate, setTime: JulianDate, beginTime: JulianDate, finishTime: JulianDate) \
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
    # todo: look at replacing
    topocentricPosition = toTopocentric(position, geo, time)
    altitude = degrees(asin(topocentricPosition[2] / topocentricPosition.mag()))
    azimuth = degrees(atan3(topocentricPosition[1], -topocentricPosition[0]))
    return PositionInfo(altitude, azimuth, time, isIlluminated, isUnobscured)


def getNextPass(satellite: Orbitable, geo: GeoPosition, timeOrPass: JulianDate | SatellitePass, timeout: float = 7) \
        -> SatellitePass:
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
    # if not isinstance(time, JulianDate):
    #     raise TypeError('time parameter must be JulianDate type')
    nextPassTime = _next_pass_max(satellite, geo, time, timeout)

    riseTime, setTime = _rise_set_times(satellite, geo, nextPassTime)
    enterTime, exitTime = getShadowTimes(satellite, nextPassTime, Shadow.PENUMBRA)
    # todo: if timezone of riseTime is different than the timezone of geo, then the rise/set times may be for the
    #   wrong day and the unobscured values may be wrong.
    # todo: idea: check sun altitude and always find the next sunset based on the angle, then find the previous rise,
    #   then the above issue shouldn't exist
    sunRiseTime, sunSetTime, _ = getSunTimes(riseTime, geo)

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
    # offset the initial small step forward between pass computations
    nTime = start.future(-0.001)

    while nTime - start < duration:
        nPass = getNextPass(satellite, geo, nTime.future(0.001))
        nTime = nPass.maxInfo.time
        if nPass.maxInfo.time - start < duration:
            passList.append(nPass)
    return tuple(passList)


def filterPassList(plist: list[SatellitePass], constraints: SatellitePassConstraints) -> list[SatellitePass]:
    filteredList = []
    for p in plist:
        altitudes = {}
        times = []
        infos = [i[1] for i in p if i[1] is not None and re.search('[a-zA-Z]+Info$', i[0])]
        if constraints.visible is not None:
            if constraints.visible is True:
                if not p.visible:
                    continue
                # do alt / time computations here
                for i in infos:
                    if i['visible'] is True:
                        altitudes[i['altitude']] = i
                        times.append(i['time']['day_number'] + i['time']['day_fraction'])
            elif constraints.visible is False and p.visible:
                continue
        if constraints.illuminated is not None:
            if constraints.illuminated is True:
                if not p.illuminated:
                    continue
                if len(altitudes) == 0:  # only do this if we haven't already
                    for i in infos:
                        if i['illuminated'] is True:
                            altitudes[i['altitude']] = i
                            times.append(i['time']['day_number'] + i['time']['day_fraction'])
            elif constraints.illuminated is False and p.illuminated:
                continue
        if constraints.unobscured is not None:
            if constraints.unobscured is True:
                if not p.unobscured:
                    continue
                if len(altitudes) == 0:  # only do this if we haven't already
                    for i in infos:
                        if i['unobscured'] is True:
                            altitudes[i['altitude']] = i
                            times.append(i['time']['day_number'] + i['time']['day_fraction'])
            elif constraints.unobscured is False and p.unobscured:
                continue

        if constraints.visible is True or constraints.illuminated is True or constraints.unobscured is True:
            highestAlt = max(altitudes.keys())
            earliestTime = min(times)
            latestTime = max(times)

            if constraints.minAltitude is not None and highestAlt < constraints.minAltitude:
                continue
            if constraints.maxAltitude is not None and highestAlt > constraints.maxAltitude:
                continue
            duration = (latestTime - earliestTime) * 1440.0
            if constraints.minDuration is not None and duration < constraints.minDuration:
                continue
            if constraints.maxDuration is not None and duration > constraints.maxDuration:
                continue

        filteredList.append(p)

    return filteredList


'''class PassController:

    def __init__(self, satellite: Satellite, geoPosition: GeoPosition, start: JulianDate, *,
                 constraints: SatellitePassConstraints = None, timeout=7, tle: TwoLineElement = None):
        self._sat = satellite
        self._geo = geoPosition
        self._initTime = start
        self._constraints = constraints
        self._timeout = timeout
        self._time = start
        if self._sat.hasTle():
            self._shadowController = ShadowController(self._sat.getTle())
        elif tle is not None:
            self._shadowController = ShadowController(tle)
        else:
            raise TLEException('No TLE given to instantiate shadow controller')

    def reset(self):
        self._time = self._initTime

    def getInitialTime(self):
        return self._initTime

    def getCurrentTime(self):
        return self._time

    def getSatellite(self):
        return self._sat

    def getGeoPosition(self):
        return self._geo

    def getConstraints(self):
        return self._constraints

    def setConstraints(self, constraints):
        self._constraints = constraints

    def getTimeout(self):
        return self._timeout

    def setTimeout(self, timeout):
        self._timeout = timeout

    def getNextPass(self) -> Pass | None:
        """Computes the next pass based on the time attribute of the PassController."""
        nextPassTime = nextPassMax(self._sat, self._geo, self._time)
        if nextPassTime is None:
            return None
        # stop searching for a pass after a time-out period
        if self._time is not None:
            # if nextPassTime.difference(self._initTime) > self._timeout:
            if nextPassTime - self._initTime > self._timeout:
                return None
        futurePassTime = nextPassTime.future(0.001)

        riseTime, setTime = riseSetTimes(self._sat, self._geo, nextPassTime)
        self._shadowController.computeValues(riseTime)
        # enterTime, exitTime = self._shadowController.getTimes()
        enterTime, exitTime = getShadowTimes(self._sat, nextPassTime, _PENUMBRA)
        # todo: if timezone of riseTime is different than the timezone of geo, then the rise/set times may be for the
        #   wrong day and the unobscured values may be wrong.
        # todo: idea: check sun altitude and always find the next sunset based on the angle, then find the previous rise,
        #   then the above issue shouldn't exist
        sunRiseTime, sunSetTime = getSunRiseSetTimes2(riseTime, self._geo)

        firstIlluminatedTime = riseTime
        lastIlluminatedTime = setTime
        #   rise and set sandwich shadow entrance
        if riseTime.value < enterTime.value < setTime.value:
            riseIlluminated = True
            setIlluminated = False
            lastIlluminatedTime = enterTime
        #   rise and set sandwich shadow exit
        elif riseTime.value < exitTime.value < setTime.value:
            riseIlluminated = False
            firstIlluminatedTime = exitTime
            setIlluminated = True
        #   rise and set are in shadow
        elif enterTime.value < riseTime.value < exitTime.value \
                and enterTime.value < setTime.value < exitTime.value:
            riseIlluminated = False
            firstIlluminatedTime = None
            setIlluminated = False
            lastIlluminatedTime = None
        #   rise and set are in sunlight
        else:
            riseIlluminated = True
            setIlluminated = True
        # check for illumination constraints
        # if self._constraints is not None:
        #     if self._constraints.illuminated is True and not (riseIlluminated or setIlluminated):
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     elif self._constraints.illuminated is False and (riseIlluminated or setIlluminated):
        #         self._time = futurePassTime
        #         return self.getNextPass()

        firstUnobscuredTime = riseTime
        lastUnobscuredTime = setTime
        #   rise and set sandwich sunrise
        if riseTime.value < sunRiseTime.value < setTime.value:
            riseUnobscured = True
            setUnobscured = False
            lastUnobscuredTime = sunRiseTime
        #   rise and set sandwich sunset
        elif riseTime.value < sunSetTime.value < setTime.value:
            riseUnobscured = False
            firstUnobscuredTime = sunSetTime
            setUnobscured = True
        #   rise and set are during the day
        elif sunRiseTime.value < riseTime.value < sunSetTime.value \
                and sunRiseTime.value < setTime.value < sunSetTime.value:
            riseUnobscured = False
            firstUnobscuredTime = None
            setUnobscured = False
            lastUnobscuredTime = None
        #   rise and set are during the night
        else:
            riseUnobscured = True
            setUnobscured = True
        # check for visibility and unobscured constraints
        # if self._constraints is not None:
        #     if self._constraints.unobscured is True and not (riseUnobscured or setUnobscured):
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     elif self._constraints.unobscured is False and (riseIlluminated or setIlluminated):
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     if self._constraints.visible is True:
        #         if firstUnobscuredTime is None or firstIlluminatedTime is None:
        #             self._time = futurePassTime
        #             return self.getNextPass()
        #         elif not (firstIlluminatedTime is firstUnobscuredTime or lastIlluminatedTime is lastUnobscuredTime) or \
        #             not (lastUnobscuredTime.value > firstIlluminatedTime.value or lastIlluminatedTime.value > firstUnobscuredTime.value):
        #             self._time = futurePassTime
        #             return self.getNextPass()
        #     elif self._constraints.visible is False:
        #         if firstUnobscuredTime is not None or firstUnobscuredTime is not None:
        #             self._time = futurePassTime
        #             return self.getNextPass()
        #         elif (firstIlluminatedTime is firstUnobscuredTime or lastIlluminatedTime is lastUnobscuredTime) or \
        #                 (lastUnobscuredTime.value > firstIlluminatedTime.value or lastIlluminatedTime.value > firstUnobscuredTime.value):
        #             self._time = futurePassTime
        #             return self.getNextPass()

        risePos = self._sat.getState(riseTime)[0]
        risePosSez = toTopocentric(risePos, riseTime, self._geo)
        riseAlt = degrees(asin(risePosSez[2] / risePosSez.mag()))
        riseInfo = PositionInfo(riseAlt,
                                degrees(atan3(risePosSez[1], -risePosSez[0])),
                                riseTime,
                                riseIlluminated,
                                riseUnobscured)
        setPos = self._sat.getState(setTime)[0]
        setPosSez = toTopocentric(setPos, setTime, self._geo)
        setAlt = degrees(asin(setPosSez[2] / setPosSez.mag()))
        setInfo = PositionInfo(setAlt,
                               degrees(atan3(setPosSez[1], -setPosSez[0])),
                               setTime,
                               setIlluminated,
                               setUnobscured)
        maxPos = self._sat.getState(nextPassTime)[0]
        maxPosSez = toTopocentric(maxPos, nextPassTime, self._geo)
        maxAlt = degrees(asin(maxPosSez[2] / maxPosSez.mag()))
        maxIlluminated = not (enterTime.value <= nextPassTime.value <= exitTime.value)
        maxUnobscured = nextPassTime.value < sunRiseTime.value or nextPassTime.value >= sunSetTime.value
        maxInfo = PositionInfo(maxAlt,
                               degrees(atan3(maxPosSez[1], -maxPosSez[0])),
                               nextPassTime,
                               maxIlluminated,
                               maxUnobscured)

        firstIlluminatedInfo, lastIlluminatedInfo = None, None
        firstIllAlt, lastIllAlt = 0, 0
        if firstIlluminatedTime is not None and firstIlluminatedTime != riseTime:
            firstIllPos = self._sat.getState(firstIlluminatedTime)[0]
            firstIllPosSez = toTopocentric(firstIllPos, firstIlluminatedTime, self._geo)
            firstIllAlt = degrees(asin(firstIllPosSez[2] / firstIllPosSez.mag()))
            firstIlluminatedInfo = PositionInfo(firstIllAlt,
                                                degrees(atan3(firstIllPosSez[1], -firstIllPosSez[0])),
                                                firstIlluminatedTime,
                                                True,
                                                firstIlluminatedTime.value < sunRiseTime.value
                                                or firstIlluminatedTime.value >= sunSetTime.value)
        if lastIlluminatedTime is not None and lastIlluminatedTime != setTime:
            lastIllPos = self._sat.getState(lastIlluminatedTime)[0]
            lastIllPosSez = toTopocentric(lastIllPos, lastIlluminatedTime, self._geo)
            lastIllAlt = degrees(asin(lastIllPosSez[2] / lastIllPosSez.mag()))
            lastIlluminatedInfo = PositionInfo(lastIllAlt,
                                               degrees(atan3(lastIllPosSez[1], -lastIllPosSez[0])),
                                               lastIlluminatedTime,
                                               True,
                                               lastIlluminatedTime.value < sunRiseTime.value
                                               or lastIlluminatedTime.value >= sunSetTime.value)
        firstUnobscuredInfo, lastUnobscuredInfo = None, None
        firstUnobAlt, lastUnobAlt = 0, 0
        if firstUnobscuredTime is not None and firstUnobscuredTime == sunSetTime:
            firstUnobPos = self._sat.getState(firstUnobscuredTime)[0]
            firstUnobPosSez = toTopocentric(firstUnobPos, firstUnobscuredTime, self._geo)
            firstUnobAlt = degrees(asin(firstUnobPosSez[2] / firstUnobPosSez.mag()))
            firstUnobscuredInfo = PositionInfo(firstUnobAlt,
                                               degrees(atan3(firstUnobPosSez[1], -firstUnobPosSez[0])),
                                               firstUnobscuredTime,
                                               firstUnobscuredTime.value < enterTime.value
                                               or firstUnobscuredTime.value >= exitTime.value,
                                               True)
        if lastUnobscuredTime is not None and lastUnobscuredTime == sunRiseTime:
            lastUnobPos = self._sat.getState(lastUnobscuredTime)[0]
            lastUnobPosSez = toTopocentric(lastUnobPos, lastUnobscuredTime, self._geo)
            lastUnobAlt = degrees(asin(lastUnobPosSez[2] / lastUnobPosSez.mag()))
            lastUnobscuredInfo = PositionInfo(lastUnobAlt,
                                              degrees(atan3(lastUnobPosSez[1], -lastUnobPosSez[0])),
                                              lastUnobscuredTime,
                                              lastUnobscuredTime.value < enterTime.value
                                              or lastUnobscuredTime.value >= exitTime.value,
                                              True)

        # if self._constraints is not None and (self._constraints.minAltitude is not None
        #                                       or self._constraints.maxAltitude is not None
        #                                       or self._constraints.minDuration is not None
        #                                       or self._constraints.maxDuration is not None):
        #     # want to find max height, then check if it's between constraint boundaries
        #     # prioritize visible, illuminated then unobscured
        # 
        #     altitudes = {}
        #     times = {}
        #     if self._constraints.visible is True:
        #         if riseIlluminated and riseUnobscured:
        #             altitudes[riseAlt] = riseInfo
        #             times[riseTime.value] = riseTime
        #         if maxIlluminated and maxUnobscured:
        #             altitudes[maxAlt] = maxInfo
        #             times[nextPassTime.value] = nextPassTime
        #         if setIlluminated and setUnobscured:
        #             altitudes[setAlt] = setInfo
        #             times[setTime.value] = setTime
        #         if firstIlluminatedInfo is not None and firstIlluminatedInfo.getVisibility():
        #             altitudes[firstIllAlt] = firstIlluminatedInfo
        #             times[firstIlluminatedTime.value] = firstIlluminatedTime
        #         if lastIlluminatedInfo is not None and lastIlluminatedInfo.getVisibility():
        #             altitudes[lastIllAlt] = lastIlluminatedInfo
        #             times[lastIlluminatedTime.value] = lastIlluminatedTime
        #         if firstUnobscuredInfo is not None and firstUnobscuredInfo.getVisibility():
        #             altitudes[firstUnobAlt] = firstUnobscuredInfo
        #             times[firstUnobscuredTime.value] = firstUnobscuredTime
        #         if lastUnobscuredInfo is not None and lastUnobscuredInfo.getVisibility():
        #             altitudes[lastUnobAlt] = lastUnobscuredInfo
        #             times[lastUnobscuredTime.value] = lastUnobscuredTime
        #         # todo: if len(altitudes) == 0, this isn't visible. maybe this is easier to check for visibility
        #         highestAlt = max(altitudes.keys())
        #         earliestTime = times[min(times.keys())]
        #         latestTime = times[max(times.keys())]
        #     elif self._constraints.illuminated is True:
        #         if riseIlluminated:
        #             altitudes[riseAlt] = riseInfo
        #             times[riseTime.value] = riseTime
        #         if maxIlluminated:
        #             altitudes[maxAlt] = maxInfo
        #             times[nextPassTime.value] = nextPassTime
        #         if setIlluminated:
        #             altitudes[setAlt] = setInfo
        #             times[setTime.value] = setTime
        #         if firstIlluminatedInfo is not None:
        #             altitudes[firstIllAlt] = firstIlluminatedInfo
        #             times[firstIlluminatedTime.value] = firstIlluminatedTime
        #         if lastIlluminatedInfo is not None:
        #             altitudes[lastIllAlt] = lastIlluminatedInfo
        #             times[lastIlluminatedTime.value] = lastIlluminatedTime
        #         if firstUnobscuredInfo is not None and firstUnobscuredInfo.getIlluminated():
        #             altitudes[firstUnobAlt] = firstUnobscuredInfo
        #             times[firstUnobscuredTime.value] = firstUnobscuredTime
        #         if lastUnobscuredInfo is not None and lastUnobscuredInfo.getIlluminated():
        #             altitudes[lastUnobAlt] = lastUnobscuredInfo
        #             times[lastUnobscuredTime.value] = lastUnobscuredTime
        #         highestAlt = max(altitudes.keys())
        #         earliestTime = times[min(times.keys())]
        #         latestTime = times[max(times.keys())]
        #     elif self._constraints.unobscured is True:
        #         if riseUnobscured:
        #             altitudes[riseAlt] = riseInfo
        #             times[riseTime.value] = riseTime
        #         if maxUnobscured:
        #             altitudes[maxAlt] = maxInfo
        #             times[nextPassTime.value] = nextPassTime
        #         if setUnobscured:
        #             altitudes[setAlt] = setInfo
        #             times[setTime.value] = setTime
        #         if firstIlluminatedInfo is not None and firstIlluminatedInfo.getUnobscured():
        #             altitudes[firstIllAlt] = firstIlluminatedInfo
        #             times[firstIlluminatedTime.value] = firstIlluminatedTime
        #         if lastIlluminatedInfo is not None and lastIlluminatedInfo.getUnobscured():
        #             altitudes[lastIllAlt] = lastIlluminatedInfo
        #             times[lastIlluminatedTime.value] = lastIlluminatedTime
        #         if firstUnobscuredInfo is not None:
        #             altitudes[firstUnobAlt] = firstUnobscuredInfo
        #             times[firstUnobscuredTime.value] = firstUnobscuredTime
        #         if lastUnobscuredInfo is not None:
        #             altitudes[lastUnobAlt] = lastUnobscuredInfo
        #             times[lastUnobscuredTime.value] = lastUnobscuredTime
        #         highestAlt = max(altitudes.keys())
        #         earliestTime = times[min(times.keys())]
        #         latestTime = times[max(times.keys())]
        #     else:
        #         highestAlt = maxAlt
        #         earliestTime = riseTime
        #         latestTime = setTime
        #     print('highestAlt:', highestAlt)
        #     print('earliestTime:', earliestTime)
        #     print('latestTime:', latestTime)
        # 
        #     if self._constraints.minAltitude is not None and highestAlt < self._constraints.minAltitude:
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     if self._constraints.maxAltitude is not None and highestAlt > self._constraints.maxAltitude:
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     duration = (latestTime.value - earliestTime.value) * 1440.0
        #     if self._constraints.minDuration is not None and duration < self._constraints.minDuration:
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     if self._constraints.maxDuration is not None and duration > self._constraints.maxDuration:
        #         self._time = futurePassTime
        #         return self.getNextPass()

        # if self._constraints.visible is True:
        #     maxAltList = [firstUnobAlt, lastUnobAlt]
        #     if maxInfo.getVisibility():
        #         maxAltList.append(maxAlt)
        #     if firstUnobscuredInfo is not None:
        #         earliestTime = firstUnobscuredTime
        #     elif self._constraints.illuminated is True and firstIlluminatedInfo is not None:
        #         earliestTime = firstIlluminatedTime
        #     else:
        #         earliestTime = riseTime
        #     if lastUnobscuredInfo is not None:
        #         latestTime = lastUnobscuredTime
        #     elif self._constraints.illuminated is True and lastIlluminatedInfo is not None:
        #         latestTime = lastIlluminatedTime
        #     else:
        #         latestTime = setTime
        # elif self._constraints.illuminated is True:
        #     maxAltList = [firstIllAlt, lastIllAlt]
        #     if maxIlluminated:
        #         maxAltList.append(maxAlt)
        #     if firstIlluminatedInfo is not None:
        #         earliestTime = firstIlluminatedTime
        #     else:
        #         earliestTime = riseTime
        #     if lastIlluminatedInfo is not None:
        #         latestTime = lastIlluminatedTime
        #     else:
        #         latestTime = setTime
        # else:
        #     maxAltList = [maxAlt]
        #     earliestTime = riseTime
        #     latestTime = setTime

        # if self._constraints.minAltitude is not None:
        #     if self._constraints.minAltitude > max(maxAltList):
        #         self._time = futurePassTime
        #         return self.getNextPass()
        # if self._constraints.maxAltitude is not None:
        #     if self._constraints.maxAltitude < max(maxAltList):
        #         self._time = futurePassTime
        #         return self.getNextPass()
        # if self._constraints.minDuration is not None:
        #     if self._constraints.minDuration < latestTime.difference(earliestTime) * 1440:
        #         self._time = futurePassTime
        #         return self.getNextPass()
        #     if self._constraints.maxDuration > latestTime.difference(earliestTime) * 1440:
        #         self._time = futurePassTime
        #         return self.getNextPass()

        # self._time = futurePassTime
        np = Pass(riseInfo, setInfo, maxInfo, firstUnobscuredInfo=firstUnobscuredInfo,
                  lastUnobscuredInfo=lastUnobscuredInfo,
                  firstIlluminatedInfo=firstIlluminatedInfo, lastIlluminatedInfo=lastIlluminatedInfo)
        self._time = futurePassTime
        if self._constraints is not None:
            filterPass = filterPassList([np], self._constraints)
            if not filterPass:
                return self.getNextPass()
            return filterPass[0]
        return np

    def getPassList(self, duration: float, start: JulianDate = None):
        tmp = self._time
        if start is None:
            self._time = self._initTime
        else:
            self._time = start

        passList = []
        nTime = self._initTime.future(-0.001)
        # while nTime.difference(self._initTime) < duration:
        while nTime - self._initTime < duration:
            nPass = self.getNextPass()
            if nPass is None:
                self._time = tmp
                return tuple(passList)
            nTime = nPass.getMaxInfo().time
            # if nPass.getMaxInfo().time.difference(self._initTime) < duration:
            if nPass.getMaxInfo().time - self._initTime < duration:
                passList.append(nPass)

        self._time = tmp
        return tuple(passList)'''


'''def nextPass(sat: Satellite, geo: GeoPosition, time: JulianDate,
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

    # this does not treat a last time before max time as the max time
    if constraints is not None and constraints.minAltitude is not None:
        #   if minimum altitude isn't attained
        if constraints.minAltitude > maxAlt:
            # print(nextPassTime, maxAlt)
            return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    sc = ShadowController(sat.getTle())
    sc.computeValues(nextPassTime)
    enterTime, exitTime = sc.getTimes()
    firstInfo, lastInfo = None, None

    riseTime, setTime = riseSetTimes(sat, geo, nextPassTime)
    # rise and set sandwich shadow entrance
    if riseTime.value < enterTime.value < setTime.value:
        riseIlluminated = True
        firstTime = riseTime
        setIlluminated = False
        lastTime = enterTime
    # rise and set sandwich shadow exit
    elif riseTime.future(1 / sat.getTle().getMeanMotion()).value < exitTime.value and \
            riseTime.value < exitTime.value < setTime.value:
        riseIlluminated = False
        firstTime = exitTime
        setIlluminated = True
        lastTime = setTime
    # rise and set are in shadow
    elif enterTime.value < riseTime.value < exitTime.value \
            and enterTime.value < setTime.value < exitTime.value:
        riseIlluminated = False
        firstTime = riseTime  # avoids setting the firstInfo object
        setIlluminated = False
        lastTime = setTime  # avoids setting the lastInfo object
    # rise and set are in sunlight
    else:
        riseIlluminated = True
        firstTime = riseTime
        setIlluminated = True
        lastTime = setTime

    # todo: find the time of sunrise/sunset, then compute when the sat is first visible including this knowledge

    if constraints is not None:
        if constraints.illuminated is not None:
            #   when constraints is illuminated and sat is not illuminated at rise and set times
            if constraints.illuminated and not riseIlluminated and not setIlluminated:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
            #   when constraints is not illuminated and sat is illuminated at some point
            if not constraints.illuminated and ((riseIlluminated or firstTime != riseTime)
                                                or (setIlluminated or lastTime != setTime)):
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)
        if constraints.minDuration is not None:
            #   when minimum duration is not met
            # if constraints.minDuration > setTime.difference(riseTime) * 1440:
            if constraints.minDuration > (setTime - riseTime) * 1440:
                return nextPass(sat, geo, nextPassTime.future(0.001), constraints)

    risePos = sat.getState(riseTime)[0]
    risePosSez = toTopocentric(risePos, riseTime, geo)
    riseInfo = PositionInfo(degrees(asin(risePosSez[2] / risePosSez.mag())),
                            degrees(atan3(risePosSez[1], -risePosSez[0])),
                            riseTime,
                            riseIlluminated and getTwilightType(riseTime, geo) > TwilightType.Day)

    setPos = sat.getState(setTime)[0]
    setPosSez = toTopocentric(setPos, setTime, geo)
    setInfo = PositionInfo(degrees(asin(setPosSez[2] / setPosSez.mag())),
                           degrees(atan3(setPosSez[1], -setPosSez[0])),
                           setTime,
                           setIlluminated and getTwilightType(setTime, geo) > TwilightType.Day)

    maxIlluminated = not isEclipsed(sat, nextPassTime)
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

    return Pass(riseInfo, setInfo, maxInfo, firstInfo=firstInfo, lastInfo=lastInfo)'''


'''def getPassList(sat: Satellite, geo: GeoPosition, start: JulianDate, duration: float,
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
    # while nTime.difference(start) < duration:
    while nTime - start < duration:
        nPass = nextPass(sat, geo, nTime.future(0.001), constraints)
        nTime = nPass.getMaxInfo().time
        # if nPass.getMaxInfo().time.difference(start) < duration:
        if nPass.getMaxInfo().time - start < duration:
            passList.append(nPass)
    return tuple(passList)'''


'''def nextPassMax(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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
    return maxPassRefine(sat, geo, nextMax)'''


'''def nextPassMaxGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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
    ta1 = vang(eccVec, pVec)
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
    while vang(state[0], pVec) > 4.85e-06:  # 1 arc second
        pVec = getPVector(geo, *state, tn)
        eccVec = computeEccentricVector(state[0], state[1])
        tan = computeTrueAnomaly(state[0], state[1])
        man = trueToMean(tan, sat.getTle().getEcc())
        tan1 = vang(eccVec, pVec)
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
    return tn'''


'''def maxPassRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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

    return jd'''


'''def riseSetTimes(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
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
    return riseTime, setTime'''


'''def riseSetGuess(sat: Satellite, geo: GeoPosition, time: JulianDate) -> tuple[JulianDate]:
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

    # zeta = norm(zenithVector(geo, time))
    zeta = geo.getZenithVector(time)
    # gamma = geoPositionVector(geo, time)
    gamma = geo.getPositionVector(time)

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
    if jd1.value < jd2.value:
        return jd1, jd2
    else:
        return jd2, jd1'''


'''def horizonTimeRefine(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
    """
    Refines a rise or set time of a satellite pass to its correct values.

    Args:
        sat: Satellite to find the overhead pass of.
        geo: The GeoPosition which the pass is viewed from.
        time: The rise or set time to be refined.

    Returns:
        The corrected rise or set time.
    """

    # needed to add a sort of proportional governor because there was a case of oscillating around the zero altitude
    p = 1.0
    iterCount = 0

    state = sat.getState(time)
    sezPos = toTopocentric(state[0], time, geo)
    alt = asin(sezPos[2] / sezPos.mag())
    while abs(alt) > radians(1 / 3600):
        iterCount += 1
        if iterCount % 10 == 0:
            p /= 2.0

        sezVel = rotateOrderTo(
            ZYX,
            EulerAngles(
                radians(geo.longitude) + earthOffsetAngle(time),
                radians(90 - geo.latitude),
                0.0
            ),
            state[1]
        )

        dz = sezPos.mag() * alt
        dt = (-dz / sezVel[2]) / 86400.0
        time = time.future(dt * p)  # p term is linear, so we can put it in the dz or dt term too

        state = sat.getState(time)
        sezPos = toTopocentric(state[0], time, geo)
        alt = asin(sezPos[2] / sezPos.mag())
    return time'''


'''def timeToPlane(sat: Satellite, geo: GeoPosition, time: JulianDate) -> JulianDate:
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
    return jd'''


'''def orbitAltitude(sat: Satellite, geo: GeoPosition, time: JulianDate) -> float:
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
    # zeta = norm(zenithVector(geo, time))
    zeta = geo.getZenithVector(time)
    # GeoPosition vector in geocentric reference frame
    # gamma = geoPositionVector(geo, time)
    gamma = geo.getPositionVector(time)
    # normalized angular momentum, vector equation for orbital plane
    lamb = cross(state[0], state[1])

    # compute intermediate values to find solution to parameterized vector intersection
    if lamb[1] != 0:
        x = dot(zeta, gamma) / (zeta[0] - (zeta[1] * lamb[0] / lamb[1]))
        y = -lamb[0] * x / lamb[1]
    else:
        y = dot(zeta, gamma) / (zeta[1] - (zeta[0] * lamb[1] / lamb[0]))
        x = -lamb[1] * y / lamb[0]
    r = EVector(x, y, 0)
    v = cross(zeta, lamb)

    # compute exact solution for parametrized vector which yields the nearest point to intersection
    t = (dot(v, gamma) - dot(v, r)) / v.mag2()
    p = v * t + r

    # find the true anomaly of this vector if it were a position vector
    trueAnom = vang(computeEccentricVector(state[0], state[1]), p)
    pSat = norm(p) * ((sma * (1 - ecc * ecc)) / (1 + ecc * cos(trueAnom)))

    ang = degrees(vang(p - gamma, pSat - gamma))
    return ang if pSat.mag2() > p.mag2() else -ang'''


'''def isEclipsed(sat: Satellite, time: JulianDate) -> bool:
    """
    Determines whether a satellite is eclipsed by the Earth.

    Args:
        sat: Satellite that might be eclipsed.
        time: Time to see if sat is eclipsed.

    Returns:
        True if the satellite is eclipsed, false otherwise.
    """

    satPos = sat.getState(time)[0]
    sunPos = getSunPosition(time)
    #   vectors relative to the satellite
    earthPos = -satPos
    relSunPos = -satPos + sunPos

    #   semi-diameters of earth and sun
    thetaE = asin(__getPerspectiveRadius(sat, time, sunPos) / earthPos.mag())
    thetaS = asin(SUN_RADIUS / relSunPos.mag())
    #   angle between earth and sun centers relative to the satellite
    theta = vang(earthPos, relSunPos)

    #   umbral eclipse
    if thetaE > thetaS and theta < (thetaE - thetaS):
        return True
    else:
        return False
    # penumbral eclipse: abs(thetaE - thetaS) < theta < (thetaE + thetaS)
    # annular eclipse: thetaS > thetaE and theta < (thetaS - thetaE)'''


'''def __getPerspectiveRadius(sat: Satellite, time: JulianDate, sunPos: EVector) -> float:
    """Computes the Earth's radius from the perspective of a satellite at a given time."""
    elements = OrbitalElements.fromTle(sat.getTle(), time)
    a = elements.sma
    ecc = elements.ecc
    inc = elements.inc
    aop = elements.aop
    f = EARTH_FLATTENING
    ae = EARTH_EQUITORIAL_RADIUS
    phi = elements.trueAnomalyAt(time)
    s = rotateOrderTo(ZXZ, EulerAngles(elements.raan, inc, aop), -norm(sunPos))

    Rz = (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
            sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - s[2] * (s[0] * cos(phi) + s[1] *
                                                                                      sin(phi)))
    fTerm = 2 * f - f * f
    mainTerm = ae * ae * (1 - fTerm)
    cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
    return ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)'''


# # noinspection PyUnresolvedReferences
# class ShadowController:
#     """
#     Controller class for computing Earth shadow entrance and exit times. The class can be initialized with a TLE or
#     a set of orbital elements, and then the computeValues() method can be called with any time to compute the shadow
#     times. The class always computes the next shadow exit time, if the time parameter is between the entrance and exit
#     times, the previous entrance time is computed, otherwise the next entrance time is computed.
#     """
#
#     def __init__(self, orbitData: TwoLineElement | OrbitalElements):
#         """Initialize the controller with the object for computing element values, either a TLE or an element object.
#         A TLE will yield much more accurate results as the elements can be updated with time."""
#
#         if type(orbitData) == TwoLineElement:
#             self._tle = orbitData
#             self._elements = [None, None]
#         elif type(orbitData) == OrbitalElements:
#             self._tle = None
#             self._elements = [orbitData, orbitData]
#         else:
#             raise TypeError('orbitData parameter must be of type TwoLineElement or OrbitalElements')
#         self._s = [None, None]
#         self._phi = [None, None]
#         self._Re = [6371, 6371]
#         self._jd = [None, None]
#         self._phi0 = None
#         self._originalTime = None
#
#     def getAnomalies(self) -> tuple[float]:
#         """Returns the most recently computed anomalies for the entrance and exit points of Earth's shadow."""
#         return tuple(self._phi)
#
#     def getTimes(self) -> tuple[JulianDate]:
#         """Returns the most recently computed times the anomalies occur for the entrance and exit points of Earth's
#         shadow."""
#         return tuple(self._jd)
#
#     def computeValues(self, time: JulianDate):
#         """
#         Computes/updates the anomalies and times of shadow entrance/exit. If the time is during transit of the shadow,
#         the previous entrance anomaly and time will be computed, otherwise the next entrance and exit points will always
#         be computed.
#
#         Args:
#             time: Time to compute the next(previous) anomalies.
#         """
#
#         if self._originalTime is None:
#             self._originalTime = time
#         self.__firstPass(time)
#         self.__loop(time)
#
#     def __loop(self, time: JulianDate):
#         """
#         Main loop which updates time dependent variables between iterations for a more precise answer.
#
#         Args:
#             time: Time to compute the next(previous) anomalies.
#
#         Will set/update the elements, s-vectors, Res values and the phi and jd values.
#         """
#
#         for i in range(2):
#             while True:
#                 self._elements[i].update(self._jd[i])
#                 self._s[i] = self.__computes(i)
#                 self._Re[i] = self.__computeRe(i)
#                 self._phi[i] = self.__computePhi(i)
#                 jd = self.__computeTimes(i)
#                 # if abs(jd.difference(self._jd[i])) < (0.1 / 86400.0):  # 0.1 of a second
#                 if abs(jd - self._jd[i]) < (0.1 / 86400.0):  # 0.1 of a second
#                     self._jd[i] = jd
#                     break
#                 else:
#                     self._jd[i] = jd
#
#         if self._originalTime.future(1 / self._tle.getMeanMotion()).value <= self._jd[1].value:
#             # middleTime = self._jd[0].future(self._jd[1].difference(self._jd[0]) / 2.0)
#             middleTime = self._jd[0].future((self._jd[1] - self._jd[0]) / 2.0)
#             self.computeValues(middleTime.future(-2.0 / self._tle.getMeanMotion()))
#         else:
#             self._originalTime = None
#         if self._jd[1].value < time.value:
#             self.computeValues(time.future((1 / self._tle.getMeanMotion()) * 0.5))
#         else:
#             self._originalTime = None
#
#     def __firstPass(self, time: JulianDate):
#         """
#         Does the necessary 'first pass' initializations before an 'official' anomaly or time is computed.
#
#         Args:
#             time: Time to compute the next(previous) anomalies.
#
#         Will set the jds, elements, s-vectors, and phis
#         """
#
#         self._jd = [time, time]
#         if self._tle:
#             elements = OrbitalElements.fromTle(self._tle, time)
#             self._elements = [elements, elements]
#         s = self.__computes(0)
#         self._s = [s, s]
#         self._phi[0] = self.__computePhi(0)
#         self._phi[1] = self.__computePhi(1)
#         self._phi0 = self._elements[0].trueAnomalyAt(time)
#         '''self._jd[0] = self.__computeTimes(0, time)
#         self._jd[1] = self.__computeTimes(1, time)'''
#         self._jd[0] = self._elements[0].timeToNextTrueAnomaly(self._phi[0], time)
#         self._jd[1] = self._elements[1].timeToNextTrueAnomaly(self._phi[1], self._jd[0])
#
#     def __computeTimes(self, index: int):
#         """Computes the time to the specified anomaly signaled by the index value."""
#         if index == 0:
#             return self._elements[0].timeToPrevTrueAnomaly(self._phi[0], self._jd[1])
#         elif index == 1:
#             return self._elements[1].timeToNextTrueAnomaly(self._phi[1], self._jd[0])
#
#     def __inShadow(self):
#         if self._phi[0] < self._phi[1]:
#             return self._phi[0] < self._phi0 <= self._phi[1]
#         else:
#             return (self._phi[1] < self._phi[0] < self._phi0) or (self._phi0 <= self._phi[1] < self._phi[0])
#
#     def __computePhi(self, index: int):
#         """Computes the true anomaly by finding the zeros of Escobal's G function."""
#         # if the phi has been computed before, use it as initial guess instead of searching for it
#         if self._phi[index] is not None:
#             phi = self.__escobalNewtonMethod(self._phi[index], index)
#             if self._s[index][0] * cos(phi) + self._s[index][1] * sin(phi) > 0:
#                 return phi
#
#         zeros = self.__getZeros(index)
#         checks = [self._s[index][0] * cos(z) + self._s[index][1] * sin(z) for z in zeros]
#         phis = [zeros[i] for i in range(4) if checks[i] > 0]
#         if len(phis) != 2:
#             raise PositiveZeroException(f'Number of positive zeros, which is {len(phis)}, should be 2')
#
#         gamma = self.__getGamma(self._s[index])
#         if index == 0:
#             if (phis[0] + gamma) % TWOPI < pi:
#                 return phis[0]
#             else:
#                 return phis[1]
#         else:
#             if (phis[0] + gamma) % TWOPI < pi:
#                 return phis[1]
#             else:
#                 return phis[0]
#
#     def __getZeros(self, index: int):
#         """Finds the zeros of Escobal's G function."""
#         frac = 0.5
#         ls = []
#         while len(ls) != 4:
#             tmp = [self.__escobalNewtonMethod(pi * i * frac, index) for i in range(int(2 / frac))]
#             for i in range(len(tmp)):
#                 if tmp[i] < 0 or tmp[i] >= TWOPI:
#                     tmp[i] %= TWOPI
#             ls = ShadowController.__removeDuplicates(tmp, 1e-4)
#             frac /= 2.0
#         return tuple(ls)
#
#     @staticmethod
#     def __removeDuplicates(list0, epsilon):
#         rtn = []
#         for l in list0:
#             addToList = True
#             for r in rtn:
#                 if abs(l - r) < epsilon:
#                     addToList = False
#                     break
#             if addToList:
#                 rtn.append(l)
#         return rtn
#
#     def __escobalGFunction(self, phi: float, index: int) -> float:
#         """Escobal's G function whose zeros represent the true anomalies of shadow entrance/exit."""
#         ecc = self._elements[index].ecc
#         c = (self._elements[index].sma * (1 - ecc * ecc)) ** 2
#         term1 = self._Re[index] * self._Re[index] * ((1 + ecc * cos(phi)) ** 2)
#         term2 = (-self._s[index][0] * cos(phi) - self._s[index][1] * sin(phi)) ** 2
#         return term1 + c * term2 - c
#
#     def __gPrime(self, phi: float, index: int) -> float:
#         """Derivative of Escobal's G function, used while iterating Newton's method."""
#         ecc = self._elements[index].ecc
#         term1 = 2 * (1 + ecc * cos(phi)) * (-ecc * sin(phi))
#         term2 = 2 * (-self._s[index][0] * cos(phi) - self._s[index][1] * sin(phi)) * (
#                 self._s[index][0] * sin(phi) - self._s[index][1] * cos(phi))
#         return self._Re[index] * self._Re[index] * term1 + ((self._elements[index].sma * (
#                 1 - ecc * ecc)) ** 2) * term2
#
#     def __escobalNewtonMethod(self, guess, index) -> float:
#         """Solves Escobal's G function by iterating with the Newton-raphson method."""
#         phi = guess
#         gi = self.__escobalGFunction(phi, index)
#         while abs(gi) > 1e-5:
#             phi = phi - gi / self.__gPrime(phi, index)
#             gi = self.__escobalGFunction(phi, index)
#
#         '''gi = self.__escobalGFunction(guess, index)
#         phi = phi - gi / self.__gPrime(phi, index)
#         while abs(gi) > 1e-5:
#             gi = self.__escobalGFunction(phi, index)
#             phi = phi - gi / self.__gPrime(phi, index)'''
#         return phi
#
#     def __computeRe(self, index: float):
#         """Computes the radius of the earth that the sun (dis)appears from."""
#         a = self._elements[index].sma
#         ecc = self._elements[index].ecc
#         inc = self._elements[index].inc
#         aop = self._elements[index].aop
#         f = EARTH_FLATTENING
#         ae = EARTH_EQUITORIAL_RADIUS
#         phi = self._phi[index]
#
#         Rz = (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
#                 sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - self._s[index][2] * (
#                     self._s[index][0] * cos(phi) + self._s[index][1] * sin(phi)))
#         fTerm = 2 * f - f * f
#         mainTerm = ae * ae * (1 - fTerm)
#         cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
#         nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
#         while abs(nextR - self._Re[index]) > 1e-7:
#             self._Re[index] = nextR
#             phi = self.__computePhi(index)
#             Rz = (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
#                     sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - self._s[index][2] * (
#                         self._s[index][0] * cos(phi) + self._s[index][1] * sin(phi)))
#             fTerm = 2 * f - f * f
#             mainTerm = ae * ae * (1 - fTerm)
#             cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
#             nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
#         return nextR
#
#     @staticmethod
#     def __getGamma(s):
#         """Computes the gamma angle used in the algorithm, the angle between the Sun projection unit vector and the
#         eccentric vector."""
#         # if s[0] < 0:
#         #     if s[1] == 0:
#         #         return 0
#         #     elif s[1] > 0:
#         #         return atan(-s[1] / s[0])
#         #     else:
#         #         return atan(s[1] / s[0])
#         # elif s[0] == 0:
#         #     return pi / 2.0
#         # elif s[0] > 0:
#         #     if s[1] > 0:
#         #         return pi - atan(s[1] / s[0])
#         #     elif s[1] == 0:
#         #         return pi
#         #     else:
#         #         return pi - atan(-s[1] / s[0])
#         return atan2(s[1], -s[0])
#
#     def __computes(self, index: int):
#         """Computes s, the unit projection vector of the Sun's position on the orbital plane."""
#         S = -norm(getSunPosition(self._jd[index]))
#         return rotateOrderTo(ZXZ, EulerAngles(self._elements[index].raan, self._elements[index].inc,
#                                               self._elements[index].aop), S)
