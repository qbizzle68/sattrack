from enum import Enum, IntEnum
from math import sqrt, radians, cos, sin, pi, atan2

from _pyevspace import Vector, getMatrixEuler, ZXZ, Angles, rotateMatrixFrom, dot

from sattrack import TWOPI, earthOffsetAngle, SIDEREAL_PER_SOLAR
from sattrack.coordinates import GeoPosition
from sattrack.orbit import Orbitable, Satellite
from sattrack.spacetime.juliandate import JulianDate
from sattrack.tle import TwoLineElement

tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23163.82406868  .00013295  00000+0  23456-3 0  9999
2 25544  51.6428 350.1551 0005377  78.9264  55.4332 15.50710892401113''')
iss = Satellite(tle)

if __debug__ is True:
    debug__all__ = ['_getZj', '_getZjPrime', '_getZjPPrime', '_getZjPPPrime', '_getZiValue', '_getZiPrime',
                    '_getZiPPrime', '_getZiPPPrime', '_getEllipseVectors', '_getConstants', '_angleDifference',
                    '_signOf'] + ['tle', 'iss']
else:
    debug__all__ = []

__all__ = ['getGeoValues', 'getEllipseVectors', 'getConstants', 'FunctionSpec', 'Point', 'Extrema', 'Boundary',
           'ALMOST_ONE', 'DOMAIN_ONE', 'Domain', 'EXTREMA_ZERO_EPSILON', 'DUPLICATE_ZERO_EPSILON', 'OrbitPathDirection',
           'ZeroIntersection', 'ZeroDisappeared', 'EPSILON', 'INITIAL_REDUCTION_COUNT', 'NEWTON_EPSILON', 'NEWTON_GAP',
           'ZeroFunction', 'OccurrenceDirection', 'SMALL', 'CHECK_MINIMUM', 'TIME_DIFFERENCE', 'SPECIAL_EPSILON',
           'ZeroCountChange', 'OrbitPath'] + debug__all__


# this is for expediting development, not intended to be in the public api
# todo: remove this
def getGeoValues(sat, geo, jd):
    angle = radians(geo.latitude)
    zk = sin(angle)
    rho = cos(angle)
    ellipseVectors = _getEllipseVectors(sat, jd)
    k = _getConstants(*ellipseVectors, geo, jd)

    return zk, rho, k


def _getZj(zi: float, rho: float, sign: int) -> float:
    assert sign == 1 or sign == -1, f'sign expected to be +/- 1, was {sign}'

    return sign * sqrt(rho * rho - zi * zi)


def _getZjPrime(zi: float, zj: float) -> float:
    return -zi / zj


def _getZjPPrime(zi: float, zj: float, zjP: float = None) -> float:
    if zjP is None:
        zjP = _getZjPrime(zi, zj)

    return (zjP * zi - zj) / (zj * zj)


def _getZjPPPrime(zi: float, zj: float, zjP: float = None, zjPP: float = None) -> float:
    if zjP is None:
        zjP = _getZjPrime(zi, zj)
    if zjPP is None:
        zjPP = _getZjPPrime(zi, zj, zjP)

    zj2 = zj * zj
    return (zj2 * (zjP + zjPP * zi - zjP) - 2 * zj * zjP * (zjP * zi - zj)) / (zj2 * zj2)


def _getZiValue(zi: float, zj: float, zk: float, k: list[float]) -> float:
    term1 = k[0] * zi * zi
    term2 = k[1] * zj * zj
    term3 = k[2] * zk * zk
    term4 = k[3] * zi * zj
    term5 = k[4] * zi * zk
    term6 = k[5] * zj * zk
    term7 = k[6] * zi
    term8 = k[7] * zj
    term9 = k[8] * zk

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 - k[9]


def _getZiPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj

    term1 = 2 * (k[0] * zi + k[1] * zj * zjP)
    term2 = k[3] * (zi * zjP + zj)
    term3 = k[4] * zk
    term4 = k[5] * zk * zjP

    return term1 + term2 + term3 + term4 + k[6] + k[7] * zjP


def _getZiPPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None, zjPP: float = None) -> float:
    if zjP is None:
        zjP = _getZjPrime(zi, zj)
    if zjPP is None:
        zjPP = _getZjPPrime(zi, zj, zjP)

    term1 = 2 * (k[0] + (k[1] * (zj * zjPP + zjP * zjP)))
    term2 = k[3] * (zi * zjPP + 2 * zjP)
    term3 = k[5] * zk * zjPP
    term4 = k[7] * zjPP

    return term1 + term2 + term3 + term4


def _getZiPPPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None, zjPP: float = None, zjPPP: float =
                  None) -> float:
    if zjP is None:
        zjP = _getZjPrime(zi, zj)
    if zjPP is None:
        zjPP = _getZjPPrime(zi, zj, zjP)
    if zjPPP is None:
        zjPPP = _getZjPPPrime(zi, zj, zjP, zjPP)

    term1 = (2 * k[1]) * (zj * zjPPP + zjP * zjPP + 2 * zjP * zjPP)
    term2 = k[3] * (zi * zjPPP + 3 * zjPP)
    term3 = k[5] * zk * zjPPP
    term4 = k[7] * zjPPP

    return term1 + term2 + term3 + term4


# Used to select the function from a FunctionSpec enumeration value.
# This methodology requires the signature of these methods to be equivalent.
_FUNCTION_ARRAY = [_getZiValue, _getZiPrime, _getZiPPrime, _getZiPPPrime]


def _getEllipseVectors(sat: Orbitable, jd: JulianDate) -> (Vector, Vector, Vector):
    elements = sat.getElements(jd)
    matrix = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    aMag = elements.sma
    bMag = aMag * sqrt(1 - elements.ecc * elements.ecc)
    cMag = aMag * elements.ecc

    a = rotateMatrixFrom(matrix, Vector(aMag, 0, 0))
    b = rotateMatrixFrom(matrix, Vector(0, bMag, 0))
    c = rotateMatrixFrom(matrix, Vector(-cMag, 0, 0))

    return a, b, c


def getEllipseVectors(sat: Orbitable, jd: JulianDate) -> (Vector, Vector, Vector):
    if not isinstance(sat, Orbitable):
        raise TypeError(f'sat must be an Orbitable type, not {type(sat)}')
    if not isinstance(jd, JulianDate):
        raise TypeError(f'jd must be a JulianDate type, not {type(jd)}')

    return _getEllipseVectors(sat, jd)


def _getConstants(a: Vector, b: Vector, c: Vector, geo: GeoPosition, jd: JulianDate) -> list[float]:
    gamma = geo.getPositionVector(jd)
    angle = radians(geo.latitude) - geo.getGeocentricLatitude()
    d = gamma.mag() * cos(angle)
    d2 = 2 * d

    k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
    k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
    k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
    k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
    k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
    k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
    k7 = d2 * c[0]
    k8 = d2 * c[1]
    k9 = d2 * c[2]
    k10 = d * d

    return [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]


def getConstants(a: Vector, b: Vector, c: Vector, geo: GeoPosition, jd: JulianDate) -> list[float]:
    if not isinstance(a, Vector):
        raise TypeError(f'a must be a Vector type, not {type(a)}')
    if not isinstance(b, Vector):
        raise TypeError(f'ba must be a Vector type, not {type(b)}')
    if not isinstance(c, Vector):
        raise TypeError(f'c must be a Vector type, not {type(c)}')
    if not isinstance(geo, GeoPosition):
        raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
    if not isinstance(jd, JulianDate):
        raise TypeError(f'jd must be a JulianDate type, not {type(jd)}')

    return _getConstants(a, b, c, geo, jd)


def _checkFloat(value, name):
    if not isinstance(value, float):
        raise TypeError(f'{name} must be a float type, not {type(value)}')


class FunctionSpec(IntEnum):
    ZERO_FUNCTION = 0
    FIRST_DERIVATIVE = 1
    SECOND_DERIVATIVE = 2
    THIRD_DERIVATIVE = 3


class Point:
    """The Point class represents a point on the geo/altitude plot.

    The x-coordinate of the point represents the zi component of the geo-position.
    The y-coordinate of the point represents the value of the zero-function.

    The point may be on any of the zero-function or its time-dependent derivatives."""

    __slots__ = '_zi', '_value', '_spec'

    def __init__(self, zi: float, value: float, spec: FunctionSpec = FunctionSpec.ZERO_FUNCTION):
        _checkFloat(zi, 'zi')
        _checkFloat(value, 'value')
        if not isinstance(spec, FunctionSpec):
            raise TypeError(f'spec must be a FunctionSpec enumeration, not {type(spec)}')

        self._zi = zi
        self._value = value
        self._spec = spec

    @classmethod
    def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, k: list[float], spec: FunctionSpec):
        # def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, k: list[float], method=_getZiValue):
        """Create a Point instance from a ZeroFunction object.

        The value of the point is computed by the function represented by the spec argument. The parameters
        to the function are retrieved or computed using the zi, func, sign and k parameters. This allows
        to avoid needing to compute values only for instantiating a Point type, which makes for cleaner
        looking code."""
        assert sign == 1 or sign == -1, f'sign expected to be +/- 1, was {sign}'

        zj = _getZj(zi, func.rho, sign)
        # value = method(zi, zj, func.zk, k)
        value = _FUNCTION_ARRAY[spec](zi, zj, func.zk, k)

        rtn = object.__new__(cls)
        rtn.__init__(zi, value, spec)

        return rtn

    def __str__(self) -> str:
        return str((self._zi, self._value, self._spec))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}{self.__str__()}'

    @property
    def zi(self):
        return self._zi

    @property
    def value(self):
        return self._value

    @property
    def spec(self):
        return self._spec


class Extrema(Point):
    """A more specific point, which represents an extrema on the geo/altitude plot.

    Like the Point class, the extrema may be a minimum or maximum of the zero-function or
    any of its immediate time-dependent derivatives."""

    pass


class Boundary:
    """The Boundary class represents a section of the geo/altitude plot, bounded on either side
    by two points, represented by Point objects.

    The boundary can represent any arbitrary bound within a valid domain, however the intention
    of the Boundary class is to represent a boundary between two extrema, defining the search area
    for a zero value.

    The two points are internally stored as two separate pairs of points: one a left and right bound,
    the other as an upper and lower bound. This is an implementation detail that is unimportant for the
    public API, however it should be known the Boundary class handles determining both types of bounds
    from a single pair of Points objects. This allows for more efficient comparisons while implementing
    a bisection procedure using the Boundary object."""

    __slots__ = '_upper', '_lower', '_left', '_right', '_middle'

    def __init__(self, point1: Point, point2: Point):

        if not isinstance(point1, Point):
            raise TypeError(f'point1 must be a Point type, not {type(point1)}')
        if not isinstance(point2, Point):
            raise TypeError(f'point2 must be a Point type, not {type(point2)}')

        if point1.value > point2.value:
            self._upper = point1
            self._lower = point2
        else:
            self._upper = point2
            self._lower = point1

        if point1.zi < point2.zi:
            self._left = point1
            self._right = point2
        else:
            self._left = point2
            self._right = point1

        self._middle = (point1.zi + point2.zi) / 2

    def __str__(self) -> str:
        return str((self._upper, self._lower))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}{self.__str__()}'

    def __contains__(self, item: float):
        _checkFloat(item, 'item')

        return self._left.zi <= item <= self._right.zi

    @property
    def left(self):
        return self._left

    @property
    def right(self):
        return self._right

    @property
    def upper(self):
        return self._upper

    @property
    def lower(self):
        return self._lower

    @property
    def middle(self):
        return self._middle


ALMOST_ONE = 0.9999999
# Used while computing extreme values in valid domains. This value must be extremely close to
# one to avoid missing a significant value between the computed extreme and the true domain limit,
# which must be avoided to also avoid singularities or imaginary numbers (complex solutions are of
# no interest to us here).
DOMAIN_ONE = 0.999999999999999


class Domain(Boundary):
    """The Domain class is derived from the Boundary class, mostly for explicitly declaring
    the boundary represents the domain as well.

    The main difference between the Domain and Boundary class form which the former is derived from,
    is the fromZeroFunction class method, which allows a Domain instance to be created from a
    ZeroFunction object."""

    def __init__(self, point1: Point, point2: Point):
        super().__init__(point1, point2)

    @classmethod
    def fromZeroFunction(cls, func: 'ZeroFunction', sign: int, k: list[float],
                         spec: FunctionSpec = FunctionSpec.ZERO_FUNCTION):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = func.rho * DOMAIN_ONE
        method = _FUNCTION_ARRAY[spec]
        leftValue = method(-limit, _getZj(-limit, func.rho, sign), func.zk, k)
        rightValue = method(limit, _getZj(limit, func.rho, sign), func.zk, k)
        left = Point(-limit, leftValue, spec)
        right = Point(limit, rightValue, spec)

        rtn = object.__new__(cls)
        rtn.__init__(left, right)

        return rtn


EXTREMA_ZERO_EPSILON = 1e-3
DUPLICATE_ZERO_EPSILON = 1e-3


class OrbitPathDirection(IntEnum):
    ASCENDING = 1
    SPECIAL = 0
    DESCENDING = -1

    @classmethod
    def fromValue(cls, primeValue: float, value: float = None):
        if not isinstance(primeValue, (int, float)):
            raise TypeError(f'value must be a numeric type, not {type(primeValue)}')

        if value is not None:
            if abs(value) < EXTREMA_ZERO_EPSILON:
                return cls.SPECIAL

        if primeValue > 0:
            return cls.ASCENDING
        else:
            return cls.DESCENDING


class ZeroIntersection:
    __slots__ = '_zi', '_sign', '_direction', '_bound'

    def __init__(self, zi: float, sign: int, direction: OrbitPathDirection, bound: Boundary):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        _checkFloat(zi, 'zi')
        if not isinstance(direction, OrbitPathDirection):
            raise TypeError(f'direction must be an enum value from OrbitPathDirection, not {type(direction)}')

        self._zi = zi
        self._sign = sign
        self._direction = direction
        self._bound = bound

    def __str__(self) -> str:
        return str((self._zi, self._sign, self._direction, self._bound))

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}{self.__str__()}'

    def __eq__(self, other):
        if isinstance(other, float):
            return abs(self._zi - other) < DUPLICATE_ZERO_EPSILON
        elif isinstance(other, ZeroIntersection):
            return abs(self._zi - other._zi) < DUPLICATE_ZERO_EPSILON and self._sign == other._sign
        return NotImplemented

    @property
    def zi(self):
        return self._zi

    @property
    def sign(self):
        return self._sign

    @zi.setter
    def zi(self, value: float):
        _checkFloat(value, 'value')

        self._zi = value

    @sign.setter
    def sign(self, value: int):
        if value != 1 and value != -1:
            raise TypeError(f'value must be +/- 1, not {value}')

        self._sign = value

    @property
    def direction(self):
        # so this attribute can be used instantiating another intersection
        # return float(self._direction)
        return self._direction

    # we shouldn't need this ?
    @direction.setter
    def direction(self, value: int):
        assert value in OrbitPathDirection, f'value must be an OrbitPathDirection, not {type(value)}'

        self._direction = value

    @property
    def boundary(self):
        return self._bound

    @boundary.setter
    def boundary(self, value):
        if not isinstance(value, Boundary):
            raise TypeError(f'value must be a {Boundary.__module__}.{Boundary.__name__}, not {type(value)}')

        self._bound = value


class ZeroDisappeared(Exception):
    def __init__(self, message, prev: JulianDate, time: JulianDate):
        super().__init__(message)
        self._prev = prev
        self._time = time

    @property
    def prev(self):
        return self._prev

    @property
    def time(self):
        return self._time


EPSILON = 1e-7
INITIAL_REDUCTION_COUNT = 4
NEWTON_EPSILON = 10
NEWTON_GAP = 1e-7


def _hasSameSign(a: (float, int), b: (float, int)) -> bool:
    """Returns True if a and b have the same sign, false otherwise.

    The function returns NotImplemented if either a or b is zero."""

    if a < 0:
        if b < 0:
            return True
        elif b > 0:
            return False
    elif a > 0:
        if b > 0:
            return True
        elif b < 0:
            return False
    return NotImplemented


class Intersections:
    """Helper class that contains all zeros for a particular function. The class also orders the zeros
    in order of occurrence of a geo-position vector throughout the Earth's rotation."""

    __slots__ = '_zeros', '_count', '_ordered'

    def __init__(self, zeros: list[ZeroIntersection], zeta: Vector):
        count = len(zeros)
        assert count in (0, 2, 4), f'expected an iterable with 0, 2, or 4 zeros, not {count}'

        if not isinstance(zeta, Vector):
            raise TypeError(f'zeta must be a {Vector.__module__}.{Vector.__name__} type, not {type(zeta)}')

        self._zeros = zeros
        self._count = len(zeros)
        self._ordered = self._orderZeros(zeros, zeta[0], zeta[1])

    @staticmethod
    def _orderZeros(zeros: list[ZeroIntersection], zi: float, zj: float) -> list[ZeroIntersection]:
        """Order the zeros based on the i and j components of the geo-position vector, so the order of the
        zeros are the order in which the geo-position encounters them throughout the next revolution."""

        (zPos := [z for z in zeros if z.sign == 1]).sort(key=lambda o: o.zi, reverse=True)
        (zNeg := [z for z in zeros if z.sign == -1]).sort(key=lambda o: o.zi)

        relativeIndex = 0
        if zj > 0 or (zj == 0 and zi > 0):
            for z in zPos:
                if z.zi > zi:
                    relativeIndex += 1
                else:
                    break
            return zPos[relativeIndex:] + zNeg + zPos[:relativeIndex]
        else:
            for z in zNeg:
                if z.zi < zi:
                    relativeIndex += 1
                else:
                    break
            return zNeg[relativeIndex:] + zPos + zNeg[:relativeIndex]

    def __len__(self):
        return self._count

    def __getitem__(self, item):
        return self._zeros[item]

    def get(self, enum: int):
        assert isinstance(enum, int), f'enum must be an integer type, not {type(enum)}'
        assert 0 <= enum < self._count, f'enum must be in [0-{self._count}), not {enum}'

        return self._ordered[enum]


class ZeroFunction:
    """The ZeroFunction class emulates the zero function we are trying to find the zeros of. This class
    defines the methods needed to find the zeros at any given time."""

    __slots__ = '_zk', '_rho', '_zeta'

    def __init__(self, geo: GeoPosition):  # , jd: JulianDate):
        # if not isinstance(geo, GeoPosition):
        #     raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
        assert isinstance(geo, GeoPosition), f'geo must be a {GeoPosition.__module__}.{GeoPosition.__name__}, not ' \
                                             f'{type(geo)}'
        # assert isinstance(jd, JulianDate), f'jd must be a {JulianDate.__module__}.{JulianDate.__name__},
        # not {type(jd)}'

        lat = radians(geo.latitude)
        self._rho = cos(lat)
        self._zk = sin(lat)

    @classmethod
    def fromValues(cls, zk: float, rho: float):
        # _checkFloat(zk, 'zk')
        # _checkFloat(rho, 'rho')

        rtn = object.__new__(cls)
        rtn._zk = zk
        rtn._rho = rho

        return rtn

    def _halfBoundary(self, bound: Boundary, sign: int, k: list[float], spec: FunctionSpec) -> Boundary:
        """Half the boundary represented by bound.

        The new bounds have the middle point as the new upper or lower bounds, whichever ensures that
        the upper and lower properties have opposite signs."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        point = Point.fromZi(bound.middle, self, sign, k, spec)

        if point.value >= 0:
            rtn = Boundary(point, bound.lower)
        else:
            rtn = Boundary(point, bound.upper)

        return rtn

    def _getZeroBifurcate(self, sign: int, k: list[float], spec: FunctionSpec) -> Extrema:
        # def _getZeroBifurcate(self, sign: int, k: list[float], function, derivative) -> Extrema:
        """Finds the zero of a function specified by spec, by continually bifurcating a boundary
        until the left and right bounds are satisfactorily close.

        The method is only guaranteed to work iff only a single zero exists. The boundary taken is the
        domain of the function, where each tail have opposite signs. If more than one zero exists, it's
        possible that during a single iteration the upper and lower bounds take the same sign, which leaves
        the upper and lower definition ambiguous.

        While it is possible to find a zero if more than one exists, it's not guaranteed it will be found.
        Additionally, if more than one zero exists, this methodology on its own is not an ideal way of
        finding all zeros. Therefore, this method should only be called when exactly one zero exists, and
        you need a foolproof way of finding it every time."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'
        assert spec.value > FunctionSpec.ZERO_FUNCTION, f'spec cannot be {FunctionSpec.ZERO_FUNCTION}'

        method = FunctionSpec(spec - 1)

        limit = self._rho * DOMAIN_ONE
        left = Point.fromZi(-limit, self, sign, k, spec)
        right = Point.fromZi(limit, self, sign, k, spec)
        bound = Boundary(left, right)

        while abs(bound.left.zi - bound.right.zi) > EPSILON:
            bound = self._halfBoundary(bound, sign, k, spec)

        return Extrema.fromZi(bound.middle, self, sign, k, method)

    def _computeZero(self, bound: Boundary, sign: int, k: list[float], spec: FunctionSpec) -> float:
        """Computes the zero of a specified function, using a boundary as a safety for ensuring a domain
        error does not occur, a specific zero is found, not any arbitrary one, and to ensure the procedure
        converges on an answer, avoiding an infinite loop by jumping around due to an unfortunate function shape.

        The sign, k and spec parameters define the function and its constants. bound is used to ensure the
        zero we are interested in is found."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'
        assert spec.value < FunctionSpec.THIRD_DERIVATIVE, f'spec cannot be FunctionSpec.THIRD_DERIVATIVE'

        # Reduce the initial boundary by halving it a set number of times. The end boundary is 1 / 2^N times
        # the original, so it only takes a few iterations to drastically reduce the search boundary. The most
        # important reason for this is the closer to the zero we start, the steeper the gradient, and the fewer
        # the number of the pitfalls of Newton's method we will run into.
        # for i in range(INITIAL_REDUCTION_COUNT):
        #     bound = self._halfBoundary(bound, sign, k, spec)
        #
        # zi0 = -2
        # zi1 = bound.middle
        # zj = _getZj(zi1, self._rho, sign)
        # function = _FUNCTION_ARRAY[spec]
        # derivative = _FUNCTION_ARRAY[spec + 1]
        #
        # # It is difficult to find a NEWTON_EPSILON that satisfies the desired accuracy, without being
        # # too restrictive for some instances. So we choose a reasonably small epsilon value, but also
        # # ensure successive zi1 values are also very close together.
        # while (value := abs(function(zi1, zj, self._zk, k))) > NEWTON_EPSILON or abs(zi1 - zi0) > NEWTON_GAP:
        #     if zi1 not in bound:
        #         bound = self._halfBoundary(bound, sign, k, spec)
        #         zi1 = bound.middle
        #     else:
        #         zi0 = zi1
        #         zj = _getZj(zi0, self._rho, sign)
        #         zi1 = zi0 - value / derivative(zi0, zj, self._zk, k)
        #
        # return zi1
        while abs(bound.left.zi - bound.right.zi) > NEWTON_GAP:
            bound = self._halfBoundary(bound, sign, k, spec)

        return bound.middle

    def _getEvenExtrema(self, sign: int, k: list[float]) -> list[Extrema]:
        """Compute the extrema points for an even number of extrema.

        The function assumes the only number of event extrema is two. The case of four
        extrema is presumed to not occur, and therefore no differentiation between the two needs
        to occur.

        The procedure is to utilize the fact that the second derivative has a single zero, with tails of
        different signs. A simple iterative bifurcation methodology will guarantee the zero will be found.
        Once found, the point is used as an extrema of the first derivative, with the accompanying tails
        of the first derivative, to generate two boundaries that surround the extrema of the zero-function.
        The zero-function extrema and tail points are then used as boundaries to check for, and find when
        they exist, the zeros of the zero-function."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get the zero from the 2nd derivative and make an extrema value
        firstExtrema = self._getZeroBifurcate(sign, k, FunctionSpec.SECOND_DERIVATIVE)

        # use the tails and the extrema to make left and right bounds of 1st derivative
        domain = Domain.fromZeroFunction(self, sign, k, FunctionSpec.FIRST_DERIVATIVE)
        leftBound = Boundary(domain.left, firstExtrema)
        rightBound = Boundary(firstExtrema, domain.right)

        # use the bounds to find the two zeros, and the extrema they represent
        leftExtremaZi = self._computeZero(leftBound, sign, k, FunctionSpec.FIRST_DERIVATIVE)
        rightExtremaZi = self._computeZero(rightBound, sign, k, FunctionSpec.FIRST_DERIVATIVE)
        leftExtrema = Extrema.fromZi(leftExtremaZi, self, sign, k, FunctionSpec.ZERO_FUNCTION)
        rightExtrema = Extrema.fromZi(rightExtremaZi, self, sign, k, FunctionSpec.ZERO_FUNCTION)

        # todo: should we return a tuple for more secure value? (i.e. cant be changed)
        return [leftExtrema, rightExtrema]

    def _getOddExtrema(self, sign: int, k: list[float]) -> list[Extrema]:
        """Compute the extrema points for an odd number of extrema.

        An extrema count is assumed to not be larger than 3, meaning the function computes
        a single or three extrema.

        The procedure utilizes the existence or absence of a pair of zeros in the second derivative.
        This can be guaranteed to be found, since the third derivative has exactly one zero and tails
        of opposite signs, an iterative bifurcation method can find this extreme of the second derivative.
        The sign of the extreme determines how many extrema of the zero-function there are. If both tails
        and extrema have the same sign, three zero-function extrema exists. If the extreme differs in sign
        to the tails only a single zero-function extrema exists.

        If a single extrema exists the same procedure can be used to bifurcate the tails of the first derivative
        to find the extrema point of the zero-function. If three extrema exist, the extreme of the second derivative,
        along with the tails, create two boundaries for the zeros of the second derivative. These bounds are
        used to find the extrema of the first derivative, which are combined with their tails to create boundaries
        for the extrema of the zero function.

        Once the boundaries of the zero-function are found, they are used to check for and find any zeros of
        the zero-function."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get zero from the 3rd derivative as an extrema value for the 2nd derivative
        secondExtrema = self._getZeroBifurcate(sign, k, FunctionSpec.THIRD_DERIVATIVE)

        # compare the signs and the 2nd derivative tails and the extrema
        domainDD = Domain.fromZeroFunction(self, sign, k, FunctionSpec.SECOND_DERIVATIVE)
        assert _hasSameSign(domainDD.left.value, domainDD.right.value) is True, \
            f'the function (sign: {sign}) expected to have an odd number of extrema'
        domainD = Domain.fromZeroFunction(self, sign, k, FunctionSpec.FIRST_DERIVATIVE)

        # if they are not the same:
        if _hasSameSign(domainDD.left.value, secondExtrema.value):
            # find the zero of the 1st derivative as an extrema value for the zero function
            zi = self._computeZero(Boundary(domainD.left, domainD.right), sign, k, FunctionSpec.FIRST_DERIVATIVE)
            return [Extrema.fromZi(zi, self, sign, k, FunctionSpec.ZERO_FUNCTION)]
        # if they are different:
        else:
            # find the zeros to the 2nd derivative as extrema for the 1st derivative
            leftDDBound = Boundary(domainDD.left, secondExtrema)
            rightDDBound = Boundary(secondExtrema, domainDD.right)
            leftDExtremaZi = self._computeZero(leftDDBound, sign, k, FunctionSpec.SECOND_DERIVATIVE)
            rightDExtremaZi = self._computeZero(rightDDBound, sign, k, FunctionSpec.SECOND_DERIVATIVE)
            leftDExtrema = Extrema.fromZi(leftDExtremaZi, self, sign, k, FunctionSpec.FIRST_DERIVATIVE)

            rightDExtrema = Extrema.fromZi(rightDExtremaZi, self, sign, k, FunctionSpec.FIRST_DERIVATIVE)

            # use tails and extrema of the 1st derivative to check for zeros
            bounds = [Boundary(domainD.left, leftDExtrema), Boundary(leftDExtrema, rightDExtrema),
                      Boundary(rightDExtrema, domainD.right)]
            # existingZeros = [bound for bound in bounds if bound.left.value * bound.right.value < 0]
            existingZeros = [bound for bound in bounds if not _hasSameSign(bound.left.value, bound.right.value)]

            # use the tails and valid extrema of the 1st derivative to get zeros as extrema values for the zero function
            extremaZis = [self._computeZero(bound, sign, k, FunctionSpec.FIRST_DERIVATIVE) for bound in existingZeros]
            return [Extrema.fromZi(zi, self, sign, k, FunctionSpec.ZERO_FUNCTION) for zi in extremaZis]

    def _computeExtrema(self, sign: int, k: list[float]) -> list[Extrema]:
        """Compute the extrema locations of the zero function."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # Determine the number of extrema in the zero function based on the derivative of the tails.
        # If they have the same sign, there are an even number of extrema, odd if they differ in sign.
        domainD = Domain.fromZeroFunction(self, sign, k, FunctionSpec.FIRST_DERIVATIVE)

        if _hasSameSign(domainD.left.value, domainD.right.value):
            return self._getEvenExtrema(sign, k)
        else:
            return self._getOddExtrema(sign, k)

    @staticmethod
    def _oppositeSign(value):
        if value > 0:
            return -1
        else:
            return 1

    def _findZeros(self, sign: int, k: list[float]) -> list[ZeroIntersection]:
        """Compute the zeros of the function depending on the sign and time, expressed as the
        constants found in k."""

        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # getExtrema of the zero function
        extrema = self._computeExtrema(sign, k)

        # get tails as extrema from a Domain
        domain = Domain.fromZeroFunction(self, sign, k, FunctionSpec.ZERO_FUNCTION)

        # get bounds for every tail/extrema
        if not extrema:
            bounds = [domain]
        else:
            bounds = [Boundary(domain.left, extrema[0])] + [Boundary(extrema[i], extrema[i + 1]) for i in
                                                            range(len(extrema) - 1)] + [
                         Boundary(extrema[-1], domain.right)]

        # find where the actual zeros exist
        existingZeros = [bound for bound in bounds if not _hasSameSign(bound.left.value, bound.right.value)]

        # find and return the zeros
        zeros = [self._computeZero(bound, sign, k, FunctionSpec.ZERO_FUNCTION) for bound in existingZeros]
        rtn = []
        for zi, bound in zip(zeros, existingZeros):
            zj = _getZj(zi, self._rho, sign)
            value = _getZiValue(zi, zj, self._zk, k)
            primeValue = _getZiPrime(zi, zj, self._zk, k)
            direction = OrbitPathDirection.fromValue(primeValue * sign * -1, value=value)
            rtn.append(ZeroIntersection(zi, sign, direction, bound))
        return rtn
        # aux = [_getZiPrime(z, _getZj(z, self._rho, sign), self._zk, k) for z in zeros]
        # return [ZeroIntersection(zi, sign, OrbitPathDirection.fromValue(ziP * sign * -1))
        #         for zi, ziP in zip(zeros, aux)]

    def _orderZeros(self, zeros: list[ZeroIntersection]) -> list[ZeroIntersection]:
        """Order the zeros based on the i and j components of the geo-position vector, so the order of the
        zeros are the order in which the geo-position encounters them throughout the next revolution."""

        zi = self._zeta[0]
        zj = self._zeta[1]
        (zPos := [z for z in zeros if z.sign == 1]).sort(key=lambda o: o.zi, reverse=True)
        (zNeg := [z for z in zeros if z.sign == -1]).sort(key=lambda o: o.zi)

        relativeIndex = 0
        if zj > 0 or (zj == 0 and zi > 0):
            for z in zPos:
                if z.zi > zi:
                    relativeIndex += 1
                else:
                    break
            return zPos[relativeIndex:] + zNeg + zPos[:relativeIndex]
        else:
            for z in zNeg:
                if z.zi < zi:
                    relativeIndex += 1
                else:
                    break
            return zNeg[relativeIndex:] + zPos + zNeg[:relativeIndex]

    # todo: we can make this argument optional in type, JulianDate will compute the k values, or directly
    #   pass k. not sure how we'll end up calling this
    def computeZeros(self, k: list[float]) -> Intersections:
        """Compute the zeros for the function from the given constant values.

        The zeros are ordered based on the occurrence during the next rotation of Earth. If the orbit path
        is currently above the horizon, the last zeros is moved to the first position."""

        # Zeros that occur in both positive and negative functions (have the same zi) refer to different
        # geo-positions (have opposite zj values) so they are not duplicates.
        return self._findZeros(1, k) + self._findZeros(-1, k)

    @property
    def rho(self):
        return self._rho

    @property
    def zk(self):
        return self._zk

    def plot(self, k: list[float]):
        if __debug__:
            import matplotlib.pyplot as plt

            limit = self._rho * DOMAIN_ONE
            plotCount = 10000
            m = int(limit * plotCount)
            x = [-limit] + [i / plotCount for i in range(-m, m)] + [limit]
            f, fPrime, fPPrime, fPPPrime = [], [], [], []
            fN, fPrimeN, fPPrimeN, fPPPrimeN = [], [], [], []

            for zi in x:
                zj = _getZj(zi, self._rho, 1)
                f.append(_getZiValue(zi, zj, self._zk, k))
                fPrime.append(_getZiPrime(zi, zj, self._zk, k))
                fPPrime.append(_getZiPPrime(zi, zj, self._zk, k))
                fPPPrime.append(_getZiPPPrime(zi, zj, self._zk, k))
                zj = _getZj(zi, self._rho, -1)
                fN.append(_getZiValue(zi, zj, self._zk, k))
                fPrimeN.append(_getZiPrime(zi, zj, self._zk, k))
                fPPrimeN.append(_getZiPPrime(zi, zj, self._zk, k))
                fPPPrimeN.append(_getZiPPPrime(zi, zj, self._zk, k))

            figure = plt.figure()
            ax = figure.add_subplot(111)
            ax.plot(x, f, '-', c='blue')
            ax.plot(x, fPrime, '-', c='dodgerblue')
            ax.plot(x, fPPrime, '-', c='powderblue')
            ax.plot(x, fPPPrime, '-', c='lightcyan')
            ax.plot(x, fN, '-', c='red')
            ax.plot(x, fPrimeN, '-', c='coral')
            ax.plot(x, fPPrimeN, '-', c='pink')
            ax.plot(x, fPPPrimeN, '-', c='mistyrose')

            maxVal = max(max(f), max(fN)) * 1.2
            minVal = min(min(f), min(fN)) * 1.2
            ax.set_ylim((minVal, maxVal))
            ax.grid()

            plt.show()
        else:
            raise NotImplemented('matplotlib is not a production dependency')


def _angleDifference(angle: float) -> float:
    """Converts an angle difference into a valid range of [-pi-pi].

    Many angle additions or subtractions may take place before validating the range of the output,
    therefore there are no restrictions to the domain for this function."""

    validAngle = angle % TWOPI
    if validAngle > pi:
        return validAngle - TWOPI
    return validAngle


def _signOf(value: float) -> int:
    """Returns 1 if value is greater than or equal to zero, false otherwise."""

    if value >= 0:
        return 1
    return -1


class OccurrenceDirection(Enum):
    """Enumeration to describe a previous or next occurrence of a cyclic event."""
    PREVIOUS_OCCURRENCE = 0
    NEXT_OCCURRENCE = 1


SMALL = 1e-7
CHECK_MINIMUM = 0.9999
TIME_DIFFERENCE = 1e-5  # less than 1/10 of a second
SPECIAL_EPSILON = 1e-5


class ZeroCountChange(Exception):
    def __init__(self, msg, sat, geo, jd):
        super().__init__(msg)
        self._sat = sat
        self._geo = geo
        self._jd = jd
        self._time = None
        self._prev = None

    @property
    def sat(self):
        return self._sat

    @property
    def geo(self):
        return self._geo

    @property
    def jd(self):
        return self._jd

    @property
    def prev(self):
        return self._prev

    @prev.setter
    def prev(self, value):
        if not isinstance(value, JulianDate):
            raise TypeError(f'value must be JulianDate type, not {type(value)}')
        self._prev = value

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        if not isinstance(value, JulianDate):
            raise TypeError(f'value must be JulianDate type, not {type(value)}')
        self._time = value


class OrbitPath:
    __slots__ = '_sat', '_geo', '_jd', '_zk', '_rho', '_zeta', '_initCount', '_above', '_func'

    def __init__(self, sat: Orbitable, geo: GeoPosition):
        if not isinstance(sat, Orbitable):
            raise TypeError(f'sat must be an Orbitable type, not {type(sat)}')
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
        # if not isinstance(jd, JulianDate):
        #     raise TypeError(f'jd must be a JulianDate type, not {type(geo)}')

        self._sat = sat
        self._geo = geo

        angle = radians(geo.latitude)
        self._zk = sin(angle)
        self._rho = cos(angle)
        self._func = ZeroFunction(geo)

        self._jd = None
        self._zeta = None
        self._initCount = None
        self._above = None

    def _orderZeros(self, zeros):

        zi = self._zeta[0]
        zj = self._zeta[1]
        (zPos := [z for z in zeros if z.sign == 1]).sort(key=lambda o: o.zi, reverse=True)
        (zNeg := [z for z in zeros if z.sign == -1]).sort(key=lambda o: o.zi)

        relativeIndex = 0
        if zj > 0 or (zj == 0 and zi > 0):
            for z in zPos:
                if z.zi > zi:
                    relativeIndex += 1
                else:
                    break
            orderedZeros = zPos[relativeIndex:] + zNeg + zPos[:relativeIndex]
        else:
            for z in zNeg:
                if z.zi < zi:
                    relativeIndex += 1
                else:
                    break
            orderedZeros = zNeg[relativeIndex:] + zPos + zNeg[:relativeIndex]

        if orderedZeros[0].direction == OrbitPathDirection.ASCENDING:
            assert orderedZeros[-1].direction == OrbitPathDirection.DESCENDING
            orderedZeros[:-1].insert(0, orderedZeros[-1])

        return orderedZeros

    def _getGeneralZeros(self, jd: JulianDate) -> list[ZeroIntersection]:
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)

        zeros = self._func.computeZeros(k)
        return self._orderZeros(zeros)

    def _checkTime(self, jd: JulianDate) -> float:
        a, b, c = _getEllipseVectors(self._sat, jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        numerator = dot(zeta, gamma - c)
        aDotZ = dot(zeta, a)
        bDotZ = dot(zeta, b)
        denominator = sqrt(aDotZ * aDotZ + bDotZ * bDotZ)

        return numerator / denominator

    # todo: do we just pass the ZeroIntersection into this?
    def _getTimeTo(self, zi: float, zj: float, jd: JulianDate, direction: OccurrenceDirection):
        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = _angleDifference(lng - radians(self._geo.longitude))
        if direction == OccurrenceDirection.NEXT_OCCURRENCE and dl < 0:
            dl += TWOPI
        elif direction == OccurrenceDirection.PREVIOUS_OCCURRENCE and dl > 0:
            dl -= TWOPI

        # We need to use sidereal day length as one revolution, not solar day length.
        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    def _computeSpecificZero(self, k: list[float], enum: int):
        assert enum in (1, 2, 3, 4), f'enum must be an integer from 1 to 4'

        # func = ZeroFunction(self._geo, self._jd)
        zeros = self._func.computeZeros(k)
        if len(zeros) != self._initCount:
            raise ZeroCountChange(f'zero count changed from {self._initCount} to {len(zeros)}',
                                  self._sat, self._geo, self._jd)

        # order zeros relative to the original zi and zj
        orderedZeros = self._orderZeros(zeros)
        # make sure the zeros aren't shifted by their procession
        # if self._above is True and orderedZeros[0].direction != OrbitPathDirection.DESCENDING or \
        #         self._above is False and orderedZeros[0].direction != OrbitPathDirection.ASCENDING:
        #     # assume the first zero when we started refining zeros has moved before the start position, so
        #     # the zeros are offset to the left by one.
        #     return orderedZeros[enum - 2]
        # else:
        #     return orderedZeros[enum - 1]
        return orderedZeros[enum - 1]

    def _switchSides(self, z: ZeroIntersection, direction: OccurrenceDirection, time: JulianDate, k) -> float:
        DT = 1 / 86400
        zj = _getZj(z.zi, self._rho, z.sign)
        moveDirection = z.sign * _signOf(_getZiPrime(z.zi, zj, self._zk, k)) * -1

        while self._checkTime(time) > 1:
            # fixme:
            #   it's possible this moves too far past both zeros. If this happens we don't really need to
            #   care about this section, but functionally we need to worry about it to avoid infinite loops.
            time = time.future(moveDirection * DT)

        return time

    def _refineZeroSpecial(self, zero: ZeroIntersection, enum: int, direction: OccurrenceDirection, e: ZeroDisappeared
                           ) -> ZeroIntersection:
        pass

    def _computeSpecialTime(self, zero: ZeroIntersection, enum: int, previousTime: JulianDate, time: JulianDate) -> JulianDate:
        middleTime = previousTime.future((time - previousTime) / 2)
        # todo: make a copy method for this
        z = ZeroIntersection(zero.zi, zero.sign, zero.direction, zero.boundary)

        if z.direction == OrbitPathDirection.ASCENDING:
            otherEnum = enum + 1
        else:
            otherEnum = enum - 1

        # middleTime = previousTime.future((time - previousTime) / 2)
        ellipseVectors = _getEllipseVectors(self._sat, previousTime)
        k = _getConstants(*ellipseVectors, self._geo, previousTime)
        otherZ = self._computeSpecificZero(k, otherEnum)

        # todo: computing the first derivative zeros (for the value maximas) in the same way we compute
        #   normal zeros (ordering and using enumerations to reference them) will be faster probably. for
        #   now we need both zeros so might as well use them to find the special point

        while abs(z.zi - otherZ.zi) > SPECIAL_EPSILON:
            middleTime = previousTime.future((time - previousTime) / 2)
            ellipseVectors = _getEllipseVectors(self._sat, middleTime)
            k = _getConstants(*ellipseVectors, self._geo, middleTime)
            try:
                z = self._computeSpecificZero(k, enum)
                otherZ = self._computeSpecificZero(k, otherEnum)
                previousTime = middleTime
            except ZeroCountChange:
                time = middleTime

        # # the best way to do this is to work towards when the value at the maxima is sufficiently small (this will work
        # # fine, just probably slower).
        # if z.direction == OrbitPathDirection.ASCENDING:
        #     # the next thing to do is to change the procedure so that if the orbit path is above (path._above is True),
        #     # make sure the first zero in the general zeros list is always ascending. This will need to be done in the
        #     # ZeroFunction class, and at all interfaces of that class in this class. Then we can safely call
        #     # (enum + 1) for an enum value if we know the current enum is an ascending zero.
        #
        # # while
        # # while z.direction != OrbitPathDirection.SPECIAL:
        #     middleTime = previousTime.future((time - previousTime) / 2)
        #     ellipseVectors = _getEllipseVectors(self._sat, middleTime)
        #     k = _getConstants(*ellipseVectors, self._geo, middleTime)
        #     try:
        #         z = self._computeSpecificZero(k, enum)
        #         previousTime = middleTime
        #     except ZeroCountChange:
        #         time = middleTime

        # middleTime might be just after the zero disappears, meaning the returned time will still raise Exceptions
        # if the time is used to find the zeros again. Need to crawl towards valid times.
        rtn = middleTime
        try:
            ellipseVectors = _getEllipseVectors(self._sat, rtn)
            k = _getConstants(*ellipseVectors, self._geo, rtn)
            _ = self._computeSpecificZero(k, enum)
            needToRefine = False
        except:
            needToRefine = True

        while needToRefine:
            # Move about 0.01 seconds at a time
            rtn = rtn.future(-1e-7)

        return rtn

    def _refineZero(self, zero: ZeroIntersection, enum: int, direction: OccurrenceDirection, start: JulianDate) -> ZeroIntersection:
        assert enum in (1, 2, 3, 4), f'enum must be an integer from 1 to 4'

        time = start
        ellipseVectors = _getEllipseVectors(self._sat, time)
        k = _getConstants(*ellipseVectors, self._geo, time)
        # todo: make a copy method for this
        z = ZeroIntersection(zero.zi, zero.sign, zero.direction, zero.boundary)
        previousTime = time.future(-1)

        # Check must be close to 1 but also less than to avoid domain errors on inverse trig functions.
        #   To avoid being too strict on the check value and inducing an infinite loop, we also check the
        #   difference in successive times computed. The more liberal CHECK_MINIMUM value is backed up by
        #   ensuring the change in time is necessarily small.
        CHECK_EPSILON = 1 - CHECK_MINIMUM
        while abs((check := self._checkTime(time)) - 1) > CHECK_EPSILON or abs(time - previousTime) > TIME_DIFFERENCE:
            ellipseVectors = _getEllipseVectors(self._sat, time)
            k = _getConstants(*ellipseVectors, self._geo, time)
            try:
                z = self._computeSpecificZero(k, enum)
            except ZeroCountChange as e:
                e.prev = previousTime
                e.time = time
                raise e
            # We need to be weary of cases where we jump past the disappearance of a zero.
            # try:
            #     z = self._computeSpecificZero(k, enum)
            # except ZeroDisappeared as e:
            #     # fixme: do we just return from here?
            #     z, time = self._refineZeroSpecial(z, enum, direction, e)
            #     break

            # Update the geo-position whose zenith vector corresponds to zi, and update the time to that point.
            zj = _getZj(z.zi, self._rho, z.sign)
            previousTime = time
            time = self._getTimeTo(z.zi, zj, self._jd, direction)

        if check > 1:
            return self._switchSides(z, direction, time, k)
        else:
            return time

    def _computeTimes(self, zeros: list[ZeroIntersection], directions: list[OccurrenceDirection]) -> list[JulianDate]:
        times = []
        for i, (z, d) in enumerate(zip(zeros, directions)):
            jd = self._jd if i == 0 else times[-1]

            try:
                time = self._refineZero(z, i + 1, d, jd)
                times.append(time)
            except ZeroCountChange as e:
                disappearTime = self._computeSpecialTime(z, i + 1, e.prev, e.time)
                if z.direction == OrbitPathDirection.ASCENDING:
                    # if directions[i] is PREVIOUS (i.e. path._above is True) then z will exist and we find it
                    if directions[i] == OccurrenceDirection.PREVIOUS_OCCURRENCE:
                        # can this actually happen? raising this exception to solve if it's ever encountered
                        raise Exception('special case occurred')
                    # if directions[i] is NEXT it needs to be checked if the zero has a solution before its
                    #   pair disappears
                    else:
                        # if zenith vector at disappearTime is 'before'
                        pass
                pass

        return times

    def computeIntersectionTimes(self, jd: JulianDate):
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a JulianDate type, not {type(jd)}')
        # set internal state dependent on the time
        self._jd = jd
        self._zeta = self._geo.getZenithVector(jd)
        self._func = ZeroFunction(self._geo)

        # set internal state dependent on initial position relative to orbit path
        zeros = self._getGeneralZeros(jd)
        self._initCount = len(zeros)
        # self._above = zeros[0].direction == OrbitPathDirection.DESCENDING
        self._above = self._checkTime(jd) < 1
        directions = [OccurrenceDirection.NEXT_OCCURRENCE] * self._initCount
        if self._above:
            # directions[-1] = OccurrenceDirection.PREVIOUS_OCCURRENCE
            directions[0] = OccurrenceDirection.PREVIOUS_OCCURRENCE

        times = self._computeTimes(zeros, directions)
        return sorted(times, key=lambda o: o.value)
