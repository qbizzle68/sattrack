from math import radians, cos, sin, pi, acos

from pyevspace import EVector, norm, cross, dot
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import trueToMean
from sattrack.util.constants import TWOPI, EARTH_MU


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
        return  cls(radians(raan), radians(inclination), radians(argumentOfPeriapsis), eccentricity, semiMajorAxis,
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
        if not isinstance(time, JulianDate):
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

        rtn = cls(radians(raan), inc, radians(aop), ecc, sma, meanAnomaly, time)

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

    def setAnomaly(self, anomaly: float, epoch: JulianDate):
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


class Orbit:
    __slots__ = '_elements',

    def __init__(self, elements: Elements, ):

