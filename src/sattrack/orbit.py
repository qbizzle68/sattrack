import json
from abc import abstractmethod, ABC
from copy import deepcopy
from inspect import Parameter, signature
from math import radians, cos, sin, pi, sqrt, atan2, floor, degrees, asin
from typing import Callable

from pyevspace import Vector, ZXZ, Z_AXIS, Angles, getMatrixEuler, rotateMatrixFrom, getMatrixAxis, \
    ReferenceFrame
from sattrack.sgp4 import TwoLineElement, getState, elementsFromState

from sattrack._topocentric import _toTopocentricState
from sattrack.coordinates import GeoPosition, getSubPoint
from sattrack.sun import getSunPosition
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import siderealTime
from sattrack._orbit import _trueToMeanAnomaly, _trueToEccentricAnomaly, _eccentricToMeanAnomaly, \
    _smaToMeanMotion, _nearestTrueAnomaly, _nearestMeanAnomaly
from sattrack.util.constants import TWOPI, EARTH_MU, EARTH_POLAR_RADIUS, EARTH_EQUITORIAL_RADIUS, SUN_MU, SUN_RADIUS, \
    EARTH_SIDEREAL_PERIOD
from sattrack.util.conversions import atan3

_all__ = ('Elements', 'Body', 'SUN_BODY', 'EARTH_BODY', 'Orbitable', 'Orbit', 'Satellite', 'meanMotionToSma',
          'smaToMeanMotion', 'meanToTrueAnomaly', 'meanToEccentricAnomaly', 'eccentricToTrueAnomaly',
          'trueToMeanAnomaly', 'trueToEccentricAnomaly', 'eccentricToMeanAnomaly')


def _elementsFromTle(tle: TwoLineElement, time: JulianDate) -> (float, float, float, float, float, float):
    """Computes the classic orbital elements directly from values of a two-line element set. These values are less
    accurate instantaneously, however they are much more stable than the oscillating values computed from the SGP4
    module."""

    dt = time - tle.epoch
    inc = tle.inc
    # tle.meanMotion is radians / minute
    n0 = tle.meanMotion / TWOPI * 1440
    # tle.ndot is radians / minute^2
    ndot = tle.ndot / TWOPI * 1440 * 1440
    nddot = tle.nddot / TWOPI * (1440 ** 3)
    dM = (dt * (n0 + dt * (ndot + dt * nddot))) * TWOPI
    meanAnomaly = (tle.meanAnomaly + dM) % TWOPI
    ecc0 = tle.ecc
    n0dot = ndot * 2
    a0 = tle.sma
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
    raan = (tle.raan + radians(lanJ2Dot + lanMoon + lanSun) * dt) % TWOPI
    aop = (tle.aop + radians(aopJ2Dot + aopMoon + aopSun) * dt) % TWOPI

    return raan, inc, aop, ecc, sma, meanAnomaly


def _checkNumericType(value, paramName):
    """Checks the value is a numeric type, with paramName being used in the error message."""
    if not isinstance(value, (int, float)):
        raise TypeError(f'{paramName} parameter must be a numeric type', type(value))


def _checkJdType(value, paramName):
    """Checks the value is a JulianDate type, with paramName being used in the error message."""
    if not isinstance(value, JulianDate):
        raise TypeError(f'{paramName} parameter must be a JulianDate type', type(value))


class Elements:
    """An object containing the classical orbital elements of a satellite at a given time. All values of the constructor
    are in radians, for degrees use the fromDegrees() class method. A true anomaly can be set if known or will be
    computed otherwise. The epoch parameter is used to compute future or past anomalies."""

    __slots__ = '_raan', '_inc', '_aop', '_ecc', '_sma', '_meanAnomaly', '_trueAnomaly', '_epoch'

    def __init__(self, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                 semiMajorAxis: float, meanAnomaly: float, epoch: JulianDate, trueAnomaly: float = None):
        _checkNumericType(raan, 'raan')
        _checkNumericType(inclination, 'inclination')
        if not 0 <= inclination <= pi:
            raise ValueError('inclination parameter must be between 0 and pi', inclination)
        _checkNumericType(argumentOfPeriapsis, 'argumentOfPeriapsis')
        _checkNumericType(eccentricity, 'eccentricity')
        _checkNumericType(semiMajorAxis, 'semi-major axis')
        _checkNumericType(meanAnomaly, 'mean anomaly')
        _checkJdType(epoch, 'epoch')
        if trueAnomaly is not None:
            _checkNumericType(trueAnomaly, 'true anomaly')
            self._trueAnomaly = trueAnomaly % TWOPI
        else:
            self._trueAnomaly = _meanToTrueAnomaly(meanAnomaly, eccentricity)

        self._raan = raan % TWOPI
        self._inc = inclination
        self._aop = argumentOfPeriapsis % TWOPI
        self._ecc = eccentricity
        self._sma = semiMajorAxis
        self._meanAnomaly = meanAnomaly % TWOPI
        self._epoch = epoch

    @classmethod
    def fromDegrees(cls, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                    semiMajorAxis: float, meanAnomaly: float, epoch: JulianDate, trueAnomaly: float = None):
        """Create an Elements object with parameter units of degrees."""

        return cls(radians(raan), radians(inclination), radians(argumentOfPeriapsis), eccentricity, semiMajorAxis,
                   radians(meanAnomaly), epoch, radians(trueAnomaly))

    @classmethod
    def fromTle(cls, tle: TwoLineElement, epoch: JulianDate):
        """Create an Elements object from the components of a TLE. The epoch parameter is used to adjust the mean
        elements from the TLE for more accurate values. These elements are less accurate instantaneously, but they are
        much more stable than the oscillating elements computed from the SGP4 algorithm."""

        if not isinstance(tle, TwoLineElement):
            raise TypeError('tle parameter must be TwoLineElement type', type(tle))
        _checkJdType(epoch, 'epoch')

        elements = _elementsFromTle(tle, epoch)
        trueAnomaly = _meanToTrueAnomaly(elements[-1], elements[3])
        return cls(*elements, epoch, trueAnomaly)

    @classmethod
    def fromState(cls, position: Vector, velocity: Vector, epoch: JulianDate, MU: float = EARTH_MU):
        """Computes an Elements object from a set of state vectors. If using state vectors from the SGP4 module, these
        will be the oscillating elements. These are much more accurate for a given instance, but also vary more between
        two different times, which can make computations dependent on them more difficult."""

        elements = elementsFromState(position, velocity, MU)
        return cls(*elements[:-1], epoch, elements[-1])

    def __str__(self):
        """Creates a string representation of the elements."""
        header = ' elements |  raan   |   inc   |   aop   |   ecc    |   sma    | mean anom ' \
                 '| true anom |        epoch         '
        values = '  values  | {:^{w}.{p}} | {:^{w}.{p}} | {:^{w}.{p}} | {:^{w}.{p}f} | {:^{sw}.{sp}} | {:^{mw}.{mp}} ' \
                 '| {:^{mw}.{mp}} | {}'.format(degrees(self._raan), degrees(self._inc), degrees(self._aop), self._ecc,
                                               self._sma, degrees(self._meanAnomaly), degrees(self._trueAnomaly),
                                               self._epoch.date(), w=7, p=6, sw=8, sp=7, mw=9, mp=6)
        return f'{header}\n{values}'

    def __repr__(self):
        """Creates a string representation of the elements."""
        return f'Elements({self._raan}, {self._inc}, {self._aop}, {self._ecc}, {self._sma}, {self._meanAnomaly},' \
               f'{repr(self._epoch)}, {self._trueAnomaly})'

    def __reduce__(self):
        """Allows for copying and pickling Elements types."""
        return self.__class__, (self._raan, self._inc, self._aop, self._ecc, self._sma, self._meanAnomaly, self._epoch)

    @property
    def raan(self):
        """Returns right-ascension in radians."""
        return self._raan

    @raan.setter
    def raan(self, value):
        """Sets right-ascension in radians."""
        _checkNumericType(value, 'raan')
        self._raan = value % TWOPI

    @property
    def inc(self):
        """Returns inclination in radians."""
        return self._inc

    @inc.setter
    def inc(self, value):
        """Sets inclination in radians."""
        _checkNumericType(value, 'inclination')
        if not 0 <= value <= pi:
            raise ValueError('inclination parameter must be between 0 and pi', value)
        self._inc = value

    @property
    def aop(self):
        """Returns argument of periapsis in radians."""
        return self._aop

    @aop.setter
    def aop(self, value):
        """Sets argument of periapsis in radians."""
        _checkNumericType(value, 'argumentOfPeriapsis')
        self._aop = value % TWOPI

    @property
    def ecc(self):
        """Returns eccentricity of the orbit."""
        return self._ecc

    @ecc.setter
    def ecc(self, value):
        """Sets the eccentricity of the orbit."""
        _checkNumericType(value, 'eccentricity')
        self._ecc = value

    @property
    def sma(self):
        """Returns the semi-major axis in kilometers."""
        return self._sma

    @sma.setter
    def sma(self, value):
        """Sets the semi-major axis in kilometers."""
        _checkNumericType(value, 'semi-major axis')
        self._sma = value

    @property
    def meanAnomaly(self):
        """Returns the mean anomaly in radians."""
        return self._meanAnomaly

    @meanAnomaly.setter
    def meanAnomaly(self, value):
        """Sets the mean anomaly in radians without changing the associated epoch."""
        _checkNumericType(value, 'mean anomaly')
        self._meanAnomaly = value

    @property
    def epoch(self):
        """Returns the epoch as a JulianDate."""
        return self._epoch

    @epoch.setter
    def epoch(self, value: JulianDate):
        """Sets the epoch of the elements without altering the anomalies."""
        _checkJdType(value, 'value')
        self._epoch = value

    @property
    def trueAnomaly(self):
        """Returns the true anomaly in radians."""
        return self._trueAnomaly

    @trueAnomaly.setter
    def trueAnomaly(self, value):
        """Sets the true anomaly in radians without changing the associated epoch."""
        _checkNumericType(value, 'true anomaly')
        self._trueAnomaly = value

    def setMeanAnomaly(self, anomaly: float, epoch: JulianDate):
        """Sets the mean anomaly in radians and changes the associated epoch at the same time. The true anomaly is
        also updated, which means if the eccentricity must be updated, it should be done before calling this method so
        the correct true anomaly can be calculated."""

        _checkNumericType(anomaly, 'anomaly')
        _checkJdType(epoch, 'epoch')
        self.meanAnomaly = anomaly
        self._trueAnomaly = _trueToMeanAnomaly(self._meanAnomaly, self._ecc)
        self.epoch = epoch

    def setTrueAnomaly(self, anomaly: float, epoch: JulianDate):
        """Sets the true anomaly in radians and changes the associated epoch at the same time. The mean anomaly is
        also updated, which means if the eccentricity must be updated, it should be done before calling this method so
        the correct mean anomaly can be calculated."""

        _checkNumericType(anomaly, 'anomaly')
        _checkJdType(epoch, 'epoch')
        self.trueAnomaly = anomaly
        self._meanAnomaly = _meanToTrueAnomaly(self._trueAnomaly, self._ecc)
        self.epoch = epoch


def _checkCallable(func, errorName):
    """Makes sure an object is a callable with 1 parameter."""

    if isinstance(func, Callable):
        params = signature(func).parameters
        if len(params) != 1:
            ls = [p for p in params.values() if p.default is Parameter.empty and p.kind in
                  (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)]
            if len(ls) != 1:
                raise TypeError(f'{errorName} must be a Callable with 1 non-default positional parameter')
    else:
        raise TypeError(f'{errorName} must be a Callable type', type(func))


def _checkStrType(value, paramName):
    """Checks the value is a string type, with paramName being used in the error message."""
    if not isinstance(value, str):
        raise TypeError(f'{paramName} parameter must be a str type', type(value))


def _checkPositiveValue(value, paramName):
    """Checks the value is positive, with paramName being used in the error message."""
    if value < 0:
        raise ValueError(f'{paramName} must be a positive number', value)


def _checkNumericPositive(value, paramName):
    """Checks the value is a numeric type and positive, with paramName being used in the error message."""
    _checkNumericType(value, paramName)
    _checkPositiveValue(value, paramName)


class Body:
    __slots__ = '_name', '_MU', '_Re', '_Rp', '_flattening', '_revPeriod', '_orbit', '_parent', \
                '_offsetFunc', '_positionFunc'

    def __init__(self, name: str, MU: float, Re: float, revPeriod: float, offsetFunc: Callable, positionFunc: Callable,
                 *, Rp: float = 0, orbit: 'Orbit' = None, parent: 'Body' = None):
        """Creates an object representing a celestial body like a plant or moon. Any unknown or non-applicable
        parameters should be set to None."""

        _checkStrType(name, 'name')
        _checkNumericPositive(MU, 'MU')
        _checkNumericPositive(Re, 'Re')
        _checkNumericPositive(revPeriod, 'revPeriod')
        _checkCallable(offsetFunc, 'offsetFunc')
        _checkCallable(positionFunc, 'positionFunc')
        if orbit is not None:
            if not isinstance(orbit, Orbit):
                raise TypeError('orbit parameter must be an Orbit type', type(orbit))
            self._orbit = orbit
        if parent is not None:
            if not isinstance(parent, Body):
                raise TypeError('parent parameter must be a Body type', type(parent))
            self._parent = parent

        self._name = name
        self._MU = MU
        self._Re = Re
        if not Rp:
            _checkNumericPositive(Rp, 'Rp')
            self._flattening = (Re - Rp) / Re
            self._Rp = Rp
        else:
            self._flattening = 0
            self._Rp = Re
        self._revPeriod = revPeriod
        self._offsetFunc = offsetFunc
        self._positionFunc = positionFunc

    @property
    def name(self):
        """Returns the name of the celestial body."""
        return self._name

    @property
    def mu(self):
        """Returns the gravitational parameter of the celestial body."""
        return self._MU

    @property
    def equitorialRadius(self):
        """Returns the equitorial radius of the celestial body."""
        return self._Re

    @property
    def polarRadius(self):
        """Returns the polar radius of the celestial body."""
        return self._Rp

    @property
    def flattening(self):
        """Returns the flattening ratio of the celestial body (Re - Rp) / Re."""
        return self._flattening

    @property
    def period(self):
        """Returns the sidereal revolution period of the celestial body."""
        return self._revPeriod

    @property
    def orbit(self):
        """Returns the Orbit of the celestial object, or None if one wasn't set."""
        return self._orbit

    @property
    def parent(self):
        """Returns the parent body if it was set, otherwise None."""
        return self._parent

    def getOffsetAngle(self, time: JulianDate) -> float:
        """Returns the offset angle of the prime-meridian of the body from the celestial x-axis at a given time."""
        return self._offsetFunc(time)

    def getPosition(self, time: JulianDate) -> Vector:
        """Returns the position of the center of the celestial body relative to the Earth's center."""
        return self._positionFunc(time)

    def getAngularMomentum(self):
        """Returns the angular momentum of the rotation of the body. This assumes that the rotation is about the
        celestial Z-axis."""
        w = TWOPI / self._revPeriod
        return Vector(0, 0, w)


"""Global Body object representing the Sun."""
SUN_BODY = Body('Sun', SUN_MU, SUN_RADIUS, 0, lambda time: 0, getSunPosition)

"""Global Body object representing the Earth."""
EARTH_BODY = Body('Earth', EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_SIDEREAL_PERIOD, siderealTime,
                  lambda time: Vector((0, 0, 0)), Rp=EARTH_POLAR_RADIUS, parent=SUN_BODY)


class Orbitable(ABC):
    """Abstract object to represent an orbitable object around a parent body. This provides the similar attributes and
    methods of the Orbit and Satellite classes. An Orbit object is an Orbitable created from an Elements object, and
    a Satellite object is an Orbitable created with a TwoLineElement object."""
    __slots__ = '_name', '_body'

    def __init__(self, name: str, body: Body = EARTH_BODY):
        """Initialize the orbitable with the common requirements name and parent body, which defaults to EARTH_BODY."""
        if not isinstance(name, str):
            raise TypeError('name parameter must be a str type')
        if not isinstance(body, Body):
            raise TypeError('body parameter must be a Body type')
        self._name = name
        self._body = body

    @property
    def name(self):
        """Returns the name of the orbitable object."""
        return self._name

    @name.setter
    def name(self, value):
        """Sets the name of the orbitable object."""
        _checkStrType(value, 'value')
        self._name = value

    @property
    def body(self):
        """Returns the parent body of the orbitable object."""
        return self._body

    @body.setter
    def body(self, value):
        """Sets the parent body of the orbitable object."""
        if not isinstance(value, Body):
            raise TypeError('value parameter must be a Body type', type(value))

    @abstractmethod
    def anomalyAt(self, time: JulianDate, anomalyType: str = 'true') -> float:
        pass

    @abstractmethod
    def timeToNextAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        pass

    @abstractmethod
    def timeToPreviousAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        pass

    @abstractmethod
    def timeToNearestAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        pass

    @abstractmethod
    def getState(self, time: JulianDate) -> (Vector, Vector):
        pass

    def getTopocentricState(self, geo: GeoPosition, time: JulianDate) -> (Vector, Vector):
        """Computes the state vectors and rotates them to a topocentric reference frame."""

        if not isinstance(geo, GeoPosition):
            raise TypeError('geo parameter must be GeoPosition type', type(geo))
        if not isinstance(time, JulianDate):
            raise TypeError('time parameter must be JulianDate type', type(time))

        state = self.getState(time)
        return _toTopocentricState(*state, geo, time)

    @abstractmethod
    def getElements(self, time: JulianDate) -> Elements:
        pass

    @abstractmethod
    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        pass

    @abstractmethod
    def getPeriapsis(self, time: JulianDate = None) -> float:
        pass

    @abstractmethod
    def getApoapsis(self, time: JulianDate = None) -> float:
        pass

    def getDetails(self, time: JulianDate, geo: GeoPosition = None):
        """Computes satellite details at a given time and geo-position."""
        state = self.getState(time)
        subPoint = getSubPoint(state[0], time)
        height = state[0].mag() - subPoint.getRadius()
        pos = [i for i in state[0]]
        vel = [i for i in state[1]]

        if geo:
            topoStateRaw = self.getTopocentricState(geo, time)
            distance = topoStateRaw[0].mag()
            alt = degrees(asin(topoStateRaw[0][2] / topoStateRaw[0].mag()))
            az = degrees(atan3(topoStateRaw[0][1], -topoStateRaw[0][0]))
            topoPosition = [i for i in topoStateRaw[0]]
            topoVelocity = [i for i in topoStateRaw[1]]
        else:
            distance, alt, az = None, None, None, None
            topoPosition, topoVelocity = None, None

        details = {'state': {'position': pos, 'velocity': vel}, 'subpoint': {'latitude': subPoint.latitude,
                                                                             'longitude': subPoint.longitude},
                   'height': height, 'topoState': {'position': topoPosition, 'velocity': topoVelocity},
                   'distance': distance, 'altitude': alt, 'azimuth': az}
        return json.dumps(details)


def _radiusAtPeriapsis(semiMajorAxis: float, eccentricity: float) -> float:
    """Computes periapsis from semi-major axis in kilometers and eccentricity."""
    return semiMajorAxis * (1 - eccentricity)


def _radiusAtApoapsis(semiMajorAxis: float, eccentricity: float) -> float:
    """Computes apoapsis from semi-major axis in kilometers and eccentricity."""
    return semiMajorAxis * (1 + eccentricity)


def _meanToEccentricAnomaly(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to an eccentric anomaly in radians."""
    Ej = meanAnomaly
    while True:
        numerator = Ej - (eccentricity * sin(Ej)) - meanAnomaly
        denominator = 1 - eccentricity * cos(Ej)
        Ej1 = Ej - (numerator / denominator)
        if abs(Ej1 - Ej) <= 1e-7:
            break
        else:
            Ej = Ej1
    return Ej1


def _eccentricToTrueAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly in radians to a true anomaly in radians."""
    beta = eccentricity / (1 + sqrt(1 - eccentricity * eccentricity))
    return eccentricAnomaly + 2 * atan2(beta * sin(eccentricAnomaly), 1 - beta * cos(eccentricAnomaly))


def _meanToTrueAnomalyNewton(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to a true anomaly in radians via Newton's method to solve Kepler's
    equation."""
    eccAnom = _meanToEccentricAnomaly(meanAnomaly, eccentricity)
    return _eccentricToTrueAnomaly(eccAnom, eccentricity)


def _meanToTrueAnomalyFast(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to a true anomaly in radians via an approximation. The function is valid for
    small eccentricity as error is on the order of (ecc^3)."""
    return meanAnomaly + 2 * eccentricity * sin(meanAnomaly) + \
        1.25 * eccentricity * eccentricity * sin(2 * meanAnomaly)


# todo: come up with an analytical rational for this number
_MAXIMUM_ECCENTRICITY = 0.17


def _meanToTrueAnomaly(meanAnomaly: float, eccentricity: float) -> float:
    """Computes a mean anomaly in radians to a true anomaly in radians."""
    if eccentricity < _MAXIMUM_ECCENTRICITY:
        return _meanToTrueAnomalyFast(meanAnomaly, eccentricity)
    return _meanToTrueAnomalyNewton(meanAnomaly, eccentricity)


def _radiusAtAnomaly(sma: float, eccentricity: float, trueAnomaly: float) -> float:
    """Computes the radius of a satellite at a given position in its orbit. Semi-major axis is in kilometers and
    trueAnomaly in radians."""
    numerator = sma * (1 - eccentricity * eccentricity)
    denominator = 1 + eccentricity * cos(trueAnomaly)
    return numerator / denominator


def _flightAngleAtAnomaly(eccentricity: float, trueAnomaly: float) -> float:
    """Computes the flight angle between velocity and position vectors at a given true anomaly in radians."""
    numerator = eccentricity * sin(trueAnomaly)
    denominator = 1 + eccentricity * cos(trueAnomaly)
    return atan2(numerator, denominator)


def _velocityAtAnomaly(sma: float, radius: float, mu: float) -> float:
    """Compute the magnitude of the velocity at a certain position, where the position is specified by its radius
    in kilometers. Semi-major axis is also in kilometers."""
    temp = (2 / radius) - (1 / sma)
    return sqrt(mu * temp)


def _nextMeanAnomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    """Finds the time a satellite next achieves a mean anomaly based on the time of a previous anomaly and mean motion.
    Anomalies are in radians and mean motion is in revolutions / day."""
    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion
    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m1 < m0:
        dm += TWOPI
    return time.future((dm / meanMotion) / TWOPI)


def _previousMeanAnomaly(meanMotion: float, m0: float, epoch0: JulianDate, m1: float, time: JulianDate) -> JulianDate:
    """Finds the time a satellite previously achieved a mean anomaly based on the time of a previous anomaly and mean
    motion. Anomalies are in radians and mean motion is in revolutions / day."""
    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion
    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m0 < m1:
        dm -= TWOPI
    return time.future((dm / meanMotion) / TWOPI)


def _nextTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate,
                     t1: float, time: JulianDate) -> JulianDate:
    """Finds the time a satellite next achieves a true anomaly based on the time of a previous anomaly and mean motion.
    Anomalies are in radians and mean motion is in revolutions / day."""
    m0 = _trueToMeanAnomaly(t0, eccentricity)
    m1 = _trueToMeanAnomaly(t1, eccentricity)
    return _nextMeanAnomaly(meanMotion, m0, epoch0, m1, time)


def _previousTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate,
                         t1: float, time: JulianDate) -> JulianDate:
    """Finds the time a satellite previously achieved a true anomaly based on the time of a previous anomaly and mean
    motion. Anomalies are in radians and mean motion is in revolutions / day."""
    m0 = _trueToMeanAnomaly(t0, eccentricity)
    m1 = _trueToMeanAnomaly(t1, eccentricity)
    return _previousMeanAnomaly(meanMotion, m0, epoch0, m1, time)


def _meanAnomalyAt(meanMotion: float, m0: float, epoch0: JulianDate, time: JulianDate) -> float:
    """Computes the mean anomaly at a given time. Mean anomaly is in radians and mean motion is in revolutions / day."""
    dM = meanMotion * (time - epoch0) * 86400.0
    return (m0 + dM) % TWOPI


def _trueAnomalyAt(meanMotion: float, eccentricity: float, t0: float, epoch0: JulianDate, time: JulianDate) -> float:
    """Computes the true anomaly at a given time. True anomaly is in radians and mean motion is in revolutions / day."""
    m0 = _trueToMeanAnomaly(t0, eccentricity)
    m1 = _meanAnomalyAt(meanMotion, m0, epoch0, time)
    return _meanToTrueAnomaly(m1, eccentricity)


def _timeToCheck(anomaly, time, anomalyType):
    _checkStrType(anomalyType, 'anomalyType')
    _checkJdType(time, 'time')
    _checkNumericType(anomaly, 'anomaly')


class Orbit(Orbitable):
    """A class derived from the Orbitable abstract class that describes an orbitable object defined by a set of
    orbital elements around a parent body."""

    __slots__ = '_elements', '_periapsis', '_apoapsis'

    def __init__(self, elements: Elements, name: str = '', body: Body = EARTH_BODY):
        """Initializes an orbit from a set of orbital elements. Future and past positions can't be computed if the
        epoch attribute was not set when instantiating and elements parameter."""

        if not isinstance(elements, Elements):
            raise TypeError('elements parameter must be an Elements type', type(elements))
        _checkStrType(name, 'name')
        if not isinstance(body, Body):
            raise TypeError('body parameter must be a Body type', type(body))
        super().__init__(name, body)
        self._elements = elements
        self._periapsis = _radiusAtPeriapsis(self._elements.sma, self._elements.ecc)
        self._apoapsis = _radiusAtApoapsis(self._elements.sma, self._elements.ecc)

    def anomalyAt(self, time: JulianDate, anomalyType: str = 'true') -> float:
        """Finds the anomaly in radians at a given time. Anomaly type must be 'true' or 'mean'."""
        _checkJdType(time, 'time')
        _checkStrType(anomalyType, 'anomalyType')

        meanMotion = _smaToMeanMotion(self._elements.sma, self._body.mu)
        if anomalyType.lower() == 'mean':
            anomaly = _meanAnomalyAt(meanMotion, self._elements.meanAnomaly, self._elements.epoch, time)
        elif anomalyType.lower() == 'true':
            t0 = _meanToTrueAnomaly(self._elements.meanAnomaly, self._elements.ecc)
            anomaly = _trueAnomalyAt(meanMotion, self._elements.ecc, t0, self._elements.epoch, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return anomaly

    def timeToNextAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the next time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'."""
        _timeToCheck(anomaly, time, anomalyType)

        meanMotion = _smaToMeanMotion(self._elements.sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            rtn = _nextTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly, self._elements.epoch,
                                   anomaly, time)
        elif _anomalyType == 'false':
            rtn = _nextMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def timeToPreviousAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the previous time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'."""
        _timeToCheck(anomaly, time, anomalyType)

        meanMotion = _smaToMeanMotion(self._elements.sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            rtn = _previousTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly,
                                       self._elements.epoch, anomaly, time)
        elif _anomalyType == 'mean':
            rtn = _previousMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def timeToNearestAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the nearest time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'. The
        time parameter here is not used, and can be None."""
        if not isinstance(anomalyType, str):
            raise TypeError('anomalyType parameter must be a str type', type(anomalyType))
        _checkNumericType(anomaly, 'anomaly')

        meanMotion = _smaToMeanMotion(self._elements.sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            rtn = _nearestTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly, self._elements.epoch,
                                      anomaly)
        elif _anomalyType == 'mean':
            rtn = _nearestMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def getState(self, time: JulianDate) -> (Vector, Vector):
        """Returns the state vectors of the orbitable at a given time in kilometers and kilometers / second."""
        _checkJdType(time, 'time')

        # trueAnomaly = self._getTrueAnomalyAt(time)
        radius = _radiusAtAnomaly(self._elements.sma, self._elements.ecc, self._elements.trueAnomaly)
        flightAngle = _flightAngleAtAnomaly(self._elements.ecc, self._elements.trueAnomaly)
        velocity = _velocityAtAnomaly(self._elements.sma, radius, self._body.mu)

        angs = Angles(self._elements.raan, self._elements.inc, self._elements.aop + self._elements.trueAnomaly)
        rotationMatrix = getMatrixEuler(ZXZ, angs)
        position = rotateMatrixFrom(rotationMatrix, Vector.e1) * radius
        # since velocity is the flight angle away from 90 degrees from the position vector, rotating the first reference
        # frame by the flight angle and taking the y-axis gives the velocity direction
        rotationMatrix = rotationMatrix @ getMatrixAxis(Z_AXIS, -flightAngle)
        velocity = rotateMatrixFrom(rotationMatrix, Vector.e2) * velocity

        return position, velocity

    def getElements(self, time: JulianDate) -> Elements:
        """Returns the orbital elements of the orbitable adjusted by the given time."""
        _checkJdType(time, 'time')

        elements = deepcopy(self._elements)
        meanMotion = _smaToMeanMotion(self._elements.sma, self._body.mu)
        elements.meanAnomaly = _meanAnomalyAt(meanMotion, self._elements.meanAnomaly, self._elements.epoch, time)
        elements.epoch = time
        return elements

    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        """Computes the perifocal reference frame of the orbitable at a given time. The y-axis is the eccentric vector,
        which points towards periapsis, the z-axis points toward the angular momentum vector and the x-axis completes
        the right-handed coordinate system."""

        return ReferenceFrame(ZXZ, Angles(self._elements.raan, self._elements.inc, self._elements.aop))

    def getPeriapsis(self, time: JulianDate = None) -> float:
        """Returns the periapsis, the lowest altitude, of the orbit."""
        return self._periapsis

    def getApoapsis(self, time: JulianDate = None) -> float:
        """Returns the apoapsis, the highest altitude, of the orbit."""
        return self._apoapsis

    # todo: test this
    def impulse(self, time: JulianDate, prograde: float, normal: float, radial: float):
        """Adjusts the orbit of an orbitable based on an impulse taking place at a given instance. Prograde is positive
        in the direction of motion, normal is positive in the direction of the angular momentum vector and radial is
        positive away from the parent body. All delta-v values are in kilometers / second."""

        position, velocity = self.getState(time)
        # x - radial out, y - prograde, z - normal
        flightAngle = _flightAngleAtAnomaly(self._elements.ecc, self._elements.trueAnomaly)
        angs = Angles(self._elements.raan, self._elements.inc, self._elements.aop +
                      self._elements.trueAnomaly - flightAngle)
        referenceFrame = ReferenceFrame(ZXZ, angs)
        progradeVector = referenceFrame.rotateFrom(Vector.e2)
        radialVector = referenceFrame.rotateFrom(Vector.e1)
        normalVector = referenceFrame.rotateFrom(Vector.e3)

        newVelocity = velocity + progradeVector * prograde + radialVector * radial + normalVector * normal
        self._elements = Elements.fromState(position, newVelocity, self._elements.epoch)
        self._periapsis = _radiusAtPeriapsis(self._elements.sma, self._elements.ecc)
        self._apoapsis = _radiusAtApoapsis(self._elements.sma, self._elements.ecc)


class Satellite(Orbitable):
    """Derived from the Orbitable class, a Satellite object is an orbitable described by a TwoLineElement object around
    a parent body."""
    __slots__ = '_tle', '_propagator'

    def __init__(self, tle: TwoLineElement):
        """Initializes the satellite object with a TLE."""
        if not isinstance(tle, TwoLineElement):
            raise TypeError('tle parameter must be a TwoLineElement type')
        super().__init__(tle.name, EARTH_BODY)
        self._tle = tle

    @property
    def tle(self):
        return self._tle

    def anomalyAt(self, time: JulianDate, anomalyType: str = 'true') -> float:
        """Finds the anomaly in radians at a given time. Anomaly type must be 'true' or 'mean'."""
        _checkJdType(time, 'time')
        _checkStrType(anomalyType, 'anomalyType')

        _, _, _, eccentricity, _, meanAnomaly = _elementsFromTle(self._tle, time)
        if anomalyType == 'true':
            anomaly = _meanToTrueAnomaly(meanAnomaly, eccentricity)
        elif anomalyType == 'mean':
            anomaly = meanAnomaly
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return anomaly

    def timeToNextAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the next time the satellite achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""
        _timeToCheck(anomaly, time, anomalyType)

        _, _, _, eccentricity, sma, meanAnomaly = _elementsFromTle(self._tle, time)
        meanMotion = _smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            trueAnomaly = _meanToTrueAnomaly(meanAnomaly, eccentricity)
            rtn = _nextTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly, time)
        elif _anomalyType == 'mean':
            rtn = _nextMeanAnomaly(meanMotion, meanAnomaly, time, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def timeToPreviousAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the previous time the orbitable achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""
        _timeToCheck(anomaly, time, anomalyType)

        _, _, _, eccentricity, sma, meanAnomaly = _elementsFromTle(self._tle, time)
        meanMotion = _smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            trueAnomaly = _meanToTrueAnomaly(meanAnomaly, eccentricity)
            rtn = _previousTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly, time)
        elif _anomalyType == 'mean':
            rtn = _previousMeanAnomaly(meanMotion, meanAnomaly, time, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def timeToNearestAnomaly(self, anomaly: float, time: JulianDate, anomalyType: str = 'true') -> JulianDate:
        """Find the nearest time the orbitable achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""
        _timeToCheck(anomaly, time, anomalyType)

        _, _, _, eccentricity, sma, meanAnomaly = _elementsFromTle(self._tle, time)
        meanMotion = _smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        _anomalyType = anomalyType.lower()
        if _anomalyType == 'true':
            trueAnomaly = _meanToTrueAnomaly(meanAnomaly, eccentricity)
            rtn = _nearestTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly)
        elif _anomalyType == 'mean':
            rtn = _nearestMeanAnomaly(meanMotion, meanAnomaly, time, anomaly)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

        return rtn

    def getState(self, time: JulianDate) -> (Vector, Vector):
        """Computes the position and velocity state vectors in kilometers and kilometers / second from the SGP4
        model."""
        _checkJdType(time, 'time')

        return getState(self._tle, time)

    def getElements(self, time: JulianDate) -> Elements:
        """Returns a set of stable orbital elements at a specified time. For oscillating elements see
        Elements.fromState()."""
        return Elements.fromTle(self._tle, time)

    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        """Generates a pyevspace.ReferenceFrame object representing the perifocal reference frame of the satellite's
        orbit."""
        #  todo: compute the vectors from eccentric and angular momentum vectors ?
        raan, inc, aop, *_ = _elementsFromTle(self._tle, time)
        return ReferenceFrame(ZXZ, Angles(raan, inc, aop))

    def getPeriapsis(self, time: JulianDate = None) -> float:
        """Computes the perigee (the smallest altitude) of the satellites orbit relative to the Earth's equitorial
        radius."""
        _checkJdType(time, 'time')
        return self._tle.perigee

    def getApoapsis(self, time: JulianDate = None) -> float:
        """Computes the apogee (the highest altitude) of the satellites orbit relative to the Earth's equitorial
        radius."""
        _checkJdType(time, 'time')
        return self._tle.apogee


def _meanMotionToSma(meanMotion: float, mu: float) -> float:
    """Converts a mean motion in radians / second to semi-major axis in kilometers."""
    oneThird = 1 / 3
    return (mu ** oneThird) / (meanMotion ** (oneThird * 2))


def meanMotionToSma(meanMotion: float, mu: float) -> float:
    """Converts a mean motion in radians / second to a semi-major axis in kilometers."""
    _checkNumericType(meanMotion, 'meanMotion')
    _checkNumericType(mu, 'mu')
    return _meanMotionToSma(meanMotion, mu)


def smaToMeanMotion(semiMajorAxis: float, mu: float) -> float:
    """Converts a semi-major axis in kilometers to a mean motion in radians / second."""
    _checkNumericType(semiMajorAxis, 'semiMajorAxis')
    _checkNumericType(mu, 'mu')
    return _smaToMeanMotion(semiMajorAxis, mu)


def meanToTrueAnomaly(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to a true anomaly in radians."""
    _checkNumericType(meanAnomaly, 'meanAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _meanToTrueAnomaly(meanAnomaly, eccentricity)


def meanToEccentricAnomaly(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to an eccentric anomaly in radians."""
    _checkNumericType(meanAnomaly, 'meanAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _meanToEccentricAnomaly(meanAnomaly, eccentricity)


def eccentricToTrueAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly in radians to a true anomaly in radians."""
    _checkNumericType(eccentricAnomaly, 'eccentricAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _eccentricToTrueAnomaly(eccentricAnomaly, eccentricity)


def trueToMeanAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly in radians to a mean anomaly in radians."""
    _checkNumericType(trueAnomaly, 'trueAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _trueToMeanAnomaly(trueAnomaly, eccentricity)


def trueToEccentricAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly in radians to an eccentric anomaly in radians."""
    _checkNumericType(trueAnomaly, 'trueAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _trueToEccentricAnomaly(trueAnomaly, eccentricity)


def eccentricToMeanAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly in radians to a mean anomaly in radians."""
    _checkNumericType(eccentricAnomaly, 'eccentricAnomaly')
    _checkNumericType(eccentricity, 'eccentricity')
    return _eccentricToMeanAnomaly(eccentricAnomaly, eccentricity)

# meanAnomalyAt
# trueAnomalyAt
# computeMeanAnomaly
# computeTrueAnomaly
