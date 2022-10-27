from inspect import Parameter, signature
from math import radians, cos, sin, pi, acos, sqrt, atan2, floor
from typing import Callable, Union

from pyevspace import EVector, norm, cross, dot
from sattrack.rotation.order import ZXZ, Axis
from sattrack.rotation.rotation import EulerAngles, getEulerMatrix, rotateMatrixFrom, getMatrix, \
    ReferenceFrame
from sattrack.sampa import getSunPosition
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import siderealTime
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import trueToMean
from sattrack.util.constants import TWOPI, EARTH_MU, EARTH_POLAR_RADIUS, EARTH_EQUITORIAL_RADIUS, SUN_MU, SUN_RADIUS
from sattrack.util.conversions import atan3


class Elements:

    __slots__ = '_raan', '_inc', '_aop', '_ecc', '_sma', '_meanAnomaly', '_epoch'

    def __init__(self, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                 semiMajorAxis: float, meanAnomaly: float, epoch: JulianDate = None):
        _check_type(raan, 'raan')
        _check_type(inclination, 'inclination')
        if not 0 <= inclination <= pi:
            raise ValueError('inclination parameter must be between 0 and 180')
        _check_type(argumentOfPeriapsis, 'argumentOfPeriapsis')
        _check_type(eccentricity, 'eccentricity')
        _check_type(semiMajorAxis, 'semi-major axis')
        _check_type(meanAnomaly, 'mean anomaly')

        self._raan = raan % TWOPI
        self._inc = inclination
        self._aop = argumentOfPeriapsis % TWOPI
        self._ecc = eccentricity
        self._sma = semiMajorAxis
        self._meanAnomaly = meanAnomaly % TWOPI
        self._epoch = epoch

        # # in degrees for readability when constructing
        # self._raan, self._inc, self._aop, self._meanAnomaly \
        #     = _check_elements(raan, inclination, argumentOfPeriapsis, eccentricity, semiMajorAxis, meanAnomaly)
        # self._epoch = epoch

    @classmethod
    def fromDegrees(cls, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                    semiMajorAxis: float, meanAnomaly: float, epoch: JulianDate = None):
        return cls(radians(raan), radians(inclination), radians(argumentOfPeriapsis), eccentricity, semiMajorAxis,
                   radians(meanAnomaly), epoch)
        #
        # _check_type(raan, 'raan')
        # _check_type(inclination, 'inclination')
        # if not 0 <= inclination <= pi:
        #     raise ValueError('inclination parameter must be between 0 and 180')
        # _check_type(argumentOfPeriapsis, 'argumentOfPeriapsis')
        # _check_type(eccentricity, 'eccentricity')
        # _check_type(semiMajorAxis, 'semi-major axis')
        # _check_type(meanAnomaly, 'mean anomaly')
        # rtn = cls

    @classmethod
    def fromTle(cls, tle: TwoLineElement, time: JulianDate = None):
        if not isinstance(tle, TwoLineElement):
            raise TypeError('tle parameter must be a TwoLineElement type')
        if time is not None and not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')

        dt = time - tle.getEpoch()
        inc = radians(tle.getInc())
        # incRad = radians(inc)
        n0 = tle.getMeanMotion()
        dM = (dt * (n0 + dt * (tle.getMeanMotionDot() + dt * tle.getMeanMotionDDot()))) * 360
        meanAnomaly = radians(tle.getMeanAnomaly() + dM) % TWOPI
        ecc0 = tle.getEcc()
        n0dot = tle.getMeanMotionDot() * 2
        a0 = tle.getSma()
        aDot = -2 * a0 * n0dot / (3 * n0)
        sma = a0 + aDot * dt
        eDot = -2 * (1 - ecc0) * n0dot / (3 * n0)
        ecc = ecc0 + eDot * dt
        temp = (a0 ** -3.5) / ((1 - (ecc0 * ecc0)) ** 2)

        # perturbations
        # non-spherical earth
        lanJ2Dot = -2.06474e14 * temp * cos(inc)
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * (sin(inc) ** 2))
        # third-body
        lanMoon = -0.00338 * cos(inc) / n0
        lanSun = -0.00154 * cos(inc) / n0
        aopMoon = 0.00169 * (4 - (5 * sin(inc) ** 2)) / n0
        aopSun = 0.00077 * (4 - (5 * sin(inc) ** 2)) / n0
        raan = (tle.getRaan() + (lanJ2Dot + lanMoon + lanSun) * dt) % 360
        aop = (tle.getAop() + (aopJ2Dot + aopMoon + aopSun) * dt) % 36

        return cls(radians(raan), inc, radians(aop), ecc, sma, meanAnomaly, time)

    @classmethod
    def fromState(cls, position: EVector, velocity: EVector, time: JulianDate = None, MU: float = EARTH_MU):
        if not isinstance(position, EVector):
            raise TypeError('position parameter must be an EVector type')
        if not isinstance(velocity, EVector):
            raise TypeError('velocity parameter must be an EVector type')
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')
        if not isinstance(MU, (int, float)):
            raise TypeError('MU parameter must be an int float type')

        angularMomentum = cross(position, velocity)
        lineOfNodes = norm(cross(EVector.e3, angularMomentum))
        eccentricityVector = _compute_eccentric_vector(position, velocity, MU)

        ecc = eccentricityVector.mag()
        inc = acos(angularMomentum[2] / angularMomentum.mag())
        raan = acos(lineOfNodes[0])
        if lineOfNodes[1] < 0:
            raan = TWOPI - raan
        aop = acos(dot(lineOfNodes, eccentricityVector) / eccentricityVector.mag())
        if eccentricityVector[2] < 0:
            aop = TWOPI - aop
        tAnom = acos(dot(eccentricityVector, position) / (eccentricityVector.mag() * position.mag()))
        if dot(position, velocity) < 0:
            tAnom = 2 * pi - tAnom
        mAnom = trueToMean(tAnom, ecc)
        sma = (angularMomentum.mag() ** 2) / ((1 - (ecc * ecc)) * MU)

        return cls(raan, inc, aop, ecc, sma, mAnom, time)

    @property
    def raan(self):
        return self._raan

    @raan.setter
    def raan(self, value):
        _check_type(value, 'raan')
        self._raan = radians(value) % TWOPI

    @property
    def inc(self):
        return self._inc

    @inc.setter
    def inc(self, value):
        _check_type(value, 'inclination')
        if not 0 <= value <= 180:
            raise ValueError('inclination parameter must be between 0 and 180')
        self._inc = radians(value)

    @property
    def aop(self):
        return self._aop

    @aop.setter
    def aop(self, value):
        _check_type(value, 'argumentOfPeriapsis')
        self._aop = radians(value) % TWOPI

    @property
    def ecc(self):
        return self._ecc

    @ecc.setter
    def ecc(self, value):
        _check_type(value, 'eccentricity')
        self._ecc = value

    @property
    def sma(self):
        return self._sma

    @sma.setter
    def sma(self, value):
        _check_type(value, 'semi-major axis')
        self._sma = value

    @property
    def meanAnomaly(self):
        return self._meanAnomaly

    @meanAnomaly.setter
    def meanAnomaly(self, value):
        _check_type(value, 'mean anomaly')
        self._meanAnomaly = radians(value)

    @property
    def epoch(self):
        return self._epoch

    @epoch.setter
    def epoch(self, value: JulianDate):
        if not isinstance(value, JulianDate):
            raise TypeError('value parameter must be JulianDate type')
        self._epoch = value

    def setMeanAnomaly(self, anomaly: float, epoch: JulianDate):
        self.meanAnomaly = anomaly
        self.epoch = epoch


def _check_type(value, paramName):
    if not isinstance(value, (int, float)):
        raise TypeError(f'{paramName} parameter must be an int or float type')


def _check_elements(raan, inc, aop, ecc, sma, ma):
    _check_type(raan, 'raan')
    _check_type(inc, 'inclination')
    if not 0 <= inc <= 180:
        raise ValueError('inclination parameter must be between 0 and 180')
    _check_type(aop, 'argumentOfPeriapsis')
    _check_type(ecc, 'eccentricity')
    _check_type(sma, 'semi-major axis')
    _check_type(ma, 'mean anomaly')
    raanRad = radians(raan) % TWOPI
    incRad = radians(inc)
    aopRad = radians(aop) % TWOPI
    maRad = radians(ma) % TWOPI
    return raanRad, incRad, aopRad, maRad


def _compute_eccentric_vector(position: EVector, velocity: EVector, MU: float) -> EVector:
    lhs = position * ((velocity.mag2() / MU) - (1 / position.mag()))
    rhs = velocity * (dot(position, velocity) / MU)
    return lhs - rhs


Func = Callable[[JulianDate], float]


class Body:
    __slots__ = '_name', '_MU', '_Re', '_Rp', '_flattening', '_revPeriod', '_orbit', '_parent', \
                '_offsetFunc', '_positionFunc'

    def __init__(self, name: str, MU: float, Re: float, revPeriod: float, offsetFunc: Func, positionFunc: Func,
                 *, Rp: float = 0, orbit: 'Orbit' = None, parent: 'Body' = None):
        self._name = name
        self._MU = MU
        self._Re = Re
        if not Rp:
            self._flattening = (Re-Rp) / Re
            self._Rp = Rp
        else:
            self._flattening = 0
            self._Rp = Re
        self._revPeriod = revPeriod
        if orbit is not None:
            if not isinstance(orbit, Orbit):
                raise TypeError('orbit parameter must be an Orbit type')
        elif parent is not None:
            if not isinstance(parent, Body):
                raise TypeError('parent parameter must be a Body type')
            self._parent = parent
        self._orbit = orbit
        # todo: bottle these into a callable method
        if isinstance(offsetFunc, Callable):
            params = signature(offsetFunc).parameters
            if len(params) != 1:
                ls = [p for p in params.values() if p.default is Parameter.empty and p.kind in
                      (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)]
                if len(ls) != 1:
                    raise TypeError('offsetFunc must be a Callable with 1 non-default positional parameter')
            self._offsetFunc = offsetFunc
        else:
            raise TypeError('offsetFunc must be a Callable type')
        if isinstance(positionFunc, Callable):
            params = signature(positionFunc).parameters
            if len(params) != 1:
                ls = [p for p in params.values() if p.default is Parameter.empty and p.kind in
                      (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)]
                if len(ls) != 1:
                    raise TypeError('positionFunc must be a Callable with 1 non-default positional parameter')
            self._positionFunc = positionFunc
        else:
            raise TypeError('positionFunc must be a Callable type')

    @property
    def name(self):
        return self._name

    @property
    def mu(self):
        return self._MU

    @property
    def equitorialRadius(self):
        return self._Re

    @property
    def polarRadius(self):
        return self._Rp

    @property
    def flattening(self):
        return self._flattening

    @property
    def period(self):
        return self._revPeriod

    @property
    def orbit(self):
        return self._orbit

    @property
    def parent(self):
        return self._parent

    def getOffsetAngle(self, time: JulianDate) -> float:
        return self._offsetFunc(time)

    def getPosition(self, time: JulianDate) -> EVector:
        return self._positionFunc(time)


SUN_BODY = Body('Sun', SUN_MU, SUN_RADIUS, 0, lambda time: 0, getSunPosition)

EARTH_BODY = Body('Earth', EARTH_MU, EARTH_EQUITORIAL_RADIUS, 86164.090531, siderealTime,
                  lambda time: EVector(0, 0, 0), Rp=EARTH_POLAR_RADIUS, parent=SUN_BODY)

AnomalyDirection = Union[int]
NEXT = AnomalyDirection(0)
PREVIOUS = AnomalyDirection(1)
NEAREST = AnomalyDirection(2)
Anomaly = Union[int]
MEAN = Anomaly(0)
TRUE = Anomaly(1)


class Orbit:
    __slots__ = '_elements', '_name', '_body', '_periapsis', '_apoapsis'

    def __init__(self, elements: Elements, name: str = "", body: Body = EARTH_BODY):
        if not isinstance(elements, Elements):
            raise TypeError('elements parameter must be an Elements type')
        if not isinstance(name, str):
            raise TypeError('name parameter must be a str type')
        if not isinstance(body, Body):
            raise TypeError('body parameter must be a Body type')
        self._elements = elements
        self._name = name
        self._body = body
        periapsis, apoapsis = self._get_apsides()
        self._periapsis = periapsis
        self._apoapsis = apoapsis

    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, value):
        if not isinstance(value, Elements):
            raise TypeError('value parameter must be an Elements type')
        self._elements = value

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, str):
            raise TypeError('value parameter must be a str type')
        self._name = value

    @property
    def body(self):
        return self._body

    @body.setter
    def body(self, value):
        if not isinstance(value, Body):
            raise TypeError('value parameter must be a Body type')

    @property
    def periapsis(self):
        return self._periapsis

    @property
    def apoapsis(self):
        return self._apoapsis

    def anomalyAt(self, time: JulianDate, anomalyType: Anomaly = TRUE) -> float:
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')
        if not isinstance(anomalyType, Anomaly):
            raise TypeError('anomalyType parameter must be an Anomaly type')
        if self._elements.epoch is None:
            raise ValueError('epoch value in the elements field must be set to compute anomalies')

        meanMotion = _sma_to_mean_motion(self._elements.sma, self._body.mu)
        if anomalyType is MEAN:
            anomaly = _mean_anomaly_at(meanMotion, self._elements.meanAnomaly, self._elements.epoch, time)
        elif anomalyType is TRUE:
            t0 = _mean_to_true_anomaly_newton(self._elements.meanAnomaly, self._elements.ecc)
            anomaly = _true_anomaly_at(meanMotion, self._elements.ecc, t0, self._elements.epoch, time)
        else:
            raise ValueError('anomalyType parameter must be MEAN or TRUE')
        return anomaly

    def timeToAnomaly(self, anomaly: float, time: JulianDate, anomalyType: Anomaly, direction: AnomalyDirection) \
            -> JulianDate:
        # anomaly - radians
        if self._elements.epoch is None:
            raise ValueError('epoch value in the elements field must be set to compute anomalies')
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')
        if not isinstance(anomalyType, Anomaly):
            raise TypeError('anomalyType parameter must be an Anomaly type')
        if not isinstance(direction, AnomalyDirection):
            raise TypeError('direction parameter must be an AnomalyDirection type')
        if direction not in (NEXT, PREVIOUS):
            raise ValueError('direction parameter must be NEXT or PREVIOUS')

        meanMotion = _sma_to_mean_motion(self._elements.sma, self._body.mu)
        if anomalyType is TRUE:
            t0 = _mean_to_true_anomaly_newton(self._elements.meanAnomaly, self._elements.ecc)
            if direction is NEXT:
                func = _next_true_anomaly
            else:
                func = _previous_true_anomaly
            return func(meanMotion, self._elements.ecc, t0, self._elements.epoch, anomaly, time)
        elif anomalyType is MEAN:
            if direction is NEXT:
                func = _next_mean_anomaly
            else:
                func = _previous_mean_anomaly
            return func(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly, time)
        else:
            raise ValueError('anomalyType must be TRUE or MEAN')

    def timeToNearestAnomaly(self, anomaly: float, anomalyType: Anomaly) -> JulianDate:
        # anomaly - radians
        if self._elements.epoch is None:
            raise ValueError('epoch value in the elements field must be set to compute anomalies')
        if not isinstance(anomalyType, Anomaly):
            raise TypeError('anomalyType parameter must be an Anomaly type')

        meanMotion = _sma_to_mean_motion(self._elements.sma, self._body.mu)
        if anomalyType is TRUE:
            t0 = _mean_to_true_anomaly_newton(self._elements.meanAnomaly, self._elements.ecc)
            return _nearest_true_anomaly(meanMotion, self._elements.ecc, t0, self._elements.epoch, anomaly)
        elif anomalyType is MEAN:
            return _nearest_mean_anomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly)
        else:
            raise ValueError('anomalyType must be TRUE or MEAN')

    def getState(self, time: JulianDate) -> (EVector, EVector):
        trueAnomaly = self._get_true_anomaly_at(time)
        radius = _radius_at_anomaly(self._elements.sma, self._elements.ecc, trueAnomaly)
        flightAngle = _flight_angle_at_anomaly(self._elements.ecc, trueAnomaly)
        velocity = _velocity_at_anomaly(self._elements.sma, radius, self._body.mu)

        rotationMatrix = getEulerMatrix(ZXZ, EulerAngles(self._elements.raan, self._elements.inc,
                                                         self._elements.aop + trueAnomaly))
        position = rotateMatrixFrom(rotationMatrix, EVector.e1) * radius
        rotationMatrix @= getMatrix(Axis.Z_AXIS, -flightAngle)
        velocity = rotateMatrixFrom(rotationMatrix, EVector.e2) * velocity

        return position, velocity

    # todo: test this
    def impulse(self, time: JulianDate, prograde: float, normal: float, radial: float):
        # prograde, normal and radial-out are positive
        position, velocity = self.getState(time)
        # x - radial out, y - prograde, z - normal
        trueAnomaly = self._get_true_anomaly_at(time)
        flightAngle = _flight_angle_at_anomaly(self._elements.ecc, trueAnomaly)
        referenceFrame = ReferenceFrame(ZXZ, EulerAngles(self._elements.raan, self._elements.inc,
                                                         self._elements.aop + trueAnomaly - flightAngle))
        progradeVector = referenceFrame.RotateFrom(EVector.e2)
        radialVector = referenceFrame.RotateFrom(EVector.e1)
        normalVector = referenceFrame.RotateFrom(EVector.e3)

        newVelocity = velocity + progradeVector * prograde + radialVector * radial + normalVector * normal
        self._elements = Elements.fromState(position, newVelocity, self._elements.epoch)
        periapsis, apoapsis = self._get_apsides()
        self._periapsis = periapsis
        self._apoapsis = apoapsis

    def _get_true_anomaly_at(self, time: JulianDate):
        if self._elements.epoch is None:
            raise ValueError('epoch value in the elements field must be set to compute state')
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be a JulianDate type')
        meanMotion = _sma_to_mean_motion(self._elements.sma, self._body.mu)
        trueAnomaly_0 = _mean_to_true_anomaly(self._elements.meanAnomaly, self._elements.ecc)
        return _true_anomaly_at(meanMotion, self._elements.ecc, trueAnomaly_0, self._elements.epoch, time)

    def _get_apsides(self):
        periapsis = self._elements.sma * (1 - self._elements.ecc)
        apoapsis = self._elements.sma * (1 + self._elements.ecc)
        return periapsis, apoapsis


def _sma_to_mean_motion(semiMajorAxis: float, mu: float) -> float:
    # semiMajorAxis - km
    return sqrt(mu / (semiMajorAxis * semiMajorAxis * semiMajorAxis))


def _mean_motion_to_sma(meanMotion: float, mu: float) -> float:
    # meanMotion - rad / s
    oneThird = 1 / 3
    return (mu ** oneThird) / (meanMotion ** (oneThird * 2))


def _mean_to_eccentric_anomaly(meanAnomaly: float, eccentricity: float) -> float:
    Ej = meanAnomaly
    while True:
        numerator = Ej - (eccentricity * sin(Ej)) - meanAnomaly
        denominator = 1 - eccentricity * cos(Ej)
        Ej1 = Ej - (numerator / denominator)
        if abs(Ej1 - Ej) <= 1e-7:
            return Ej1
        else:
            Ej = Ej1


def _eccentric_to_true_anomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    beta = eccentricity / (1 + sqrt(1 - eccentricity * eccentricity))
    return eccentricAnomaly + 2 * atan2(beta * sin(eccentricAnomaly), 1 - beta * cos(eccentricAnomaly))


def _mean_to_true_anomaly_newton(meanAnomaly: float, eccentricity: float) -> float:
    eccAnom = _mean_to_eccentric_anomaly(meanAnomaly, eccentricity)
    return _eccentric_to_true_anomaly(eccAnom, eccentricity)


def _mean_to_true_anomaly_fast(meanAnomaly: float, eccentricity: float) -> float:
    return meanAnomaly + 2 * eccentricity * sin(meanAnomaly) \
           + 1.25 * eccentricity * eccentricity * sin(2 * meanAnomaly)


# todo: come up with an analytical rational for this number
_MAXIMUM_ECCENTRICITY = 0.17


def _mean_to_true_anomaly(meanAnomaly: float, eccentricity: float) -> float:
    if eccentricity < _MAXIMUM_ECCENTRICITY:
        return _mean_to_true_anomaly_fast(meanAnomaly, eccentricity)
    return _mean_to_true_anomaly_newton(meanAnomaly, eccentricity)


def _true_to_eccentric_anomaly(trueAnomaly: float, eccentricity: float) -> float:
    y = sqrt(1 - (eccentricity * eccentricity)) * sin(trueAnomaly)
    return atan3(y, cos(trueAnomaly) + eccentricity)


def _eccentric_to_mean_anomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    return eccentricAnomaly - eccentricity * sin(eccentricAnomaly)


def _true_to_mean_anomaly(trueAnomaly: float, eccentricity: float) -> float:
    eccentricAnomaly = _true_to_eccentric_anomaly(trueAnomaly, eccentricity)
    return _eccentric_to_mean_anomaly(eccentricAnomaly, eccentricity)


def _next_mean_anomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    # meanMotion - rev / day; m0, m1 - radians
    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion
    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m1 < m0:
        dm += TWOPI
    return time.future((dm / meanMotion) / TWOPI)


def _previous_mean_anomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    # meanMotion - rev / day; m0, m1 - radians
    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion
    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m0 < m1:
        dm -= TWOPI
    return time.future((dm / meanMotion) / TWOPI)


def _next_true_anomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate,
                       t1: float, time: JulianDate) -> JulianDate:
    # meanMotion - rev / day; m0, m1 - radians
    m0 = _true_to_mean_anomaly(t0, eccentricity)
    m1 = _true_to_mean_anomaly(t1, eccentricity)
    return _next_mean_anomaly(meanMotion, m0, epoch0, m1, time)


def _previous_true_anomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate,
                           t1: float, time: JulianDate) -> JulianDate:
    # meanMotion - rev / day; m0, m1 - radians
    m0 = _true_to_mean_anomaly(t0, eccentricity)
    m1 = _true_to_mean_anomaly(t1, eccentricity)
    return _previous_mean_anomaly(meanMotion, m0, epoch0, m1, time)


def _mean_anomaly_at(meanMotion: float, m0: float, epoch0: JulianDate, time: JulianDate) -> float:
    # meanMotion - rev / day; m0 radians
    dM = meanMotion * (time - epoch0) * 86400.0
    return (m0 + dM) % TWOPI


def _true_anomaly_at(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate, time: JulianDate) -> float:
    # meanMotion - rev / day; t0 - radians
    m0 = _true_to_mean_anomaly(t0, eccentricity)
    m1 = _mean_anomaly_at(meanMotion, m0, epoch0, time)
    return _mean_to_true_anomaly_newton(m1, eccentricity)


def _nearest_mean_anomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float) -> JulianDate:
    if m1 < m0:
        if m0 - m1 < pi:
            dma = m1 - m0
        else:
            dma = m1 + TWOPI - m0
    else:
        if m1 - m0 > pi:
            dma = m1 - TWOPI - m0
        else:
            dma = m1 - m0
    return epoch0.future(dma / (meanMotion * 86400))


def _nearest_true_anomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate, t1: float) \
        -> JulianDate:
    m0 = _true_to_mean_anomaly(t0, eccentricity)
    m1 = _true_to_mean_anomaly(t1, eccentricity)
    return _nearest_mean_anomaly(meanMotion, m0, epoch0, m1)


def _radius_at_anomaly(sma: float, eccentricity: float, trueAnomaly: float) -> float:
    numerator = sma * (1 - eccentricity * eccentricity)
    denominator = 1 + eccentricity * cos(trueAnomaly)
    return numerator / denominator


def _flight_angle_at_anomaly(eccentricity: float, trueAnomaly: float) -> float:
    numerator = eccentricity * sin(trueAnomaly)
    denominator = 1 + eccentricity * cos(trueAnomaly)
    return atan2(numerator, denominator)


def _velocity_at_anomaly(sma: float, radius: float, mu: float) -> float:
    temp = (2 / radius) - (1 / sma)
    return sqrt(mu * temp)
