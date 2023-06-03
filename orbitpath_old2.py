from math import sqrt, cos, radians, pi, atan2

from _pyevspace import getMatrixEuler, ZXZ, Angles, rotateMatrixFrom, Vector, dot

from sattrack.coordinates import GeoPosition
from sattrack.orbit import Orbitable
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.constants import TWOPI, SIDEREAL_PER_SOLAR

from sattrack.api import *
tle = TwoLineElement('''ISS (ZARYA)
1 25544U 98067A   23145.82312521  .00013516  00000+0  24287-3 0  9990
2 25544  51.6416  79.3943 0005436  19.5209  14.6217 15.50172027398320''')
iss = Satellite(tle)
tle = TwoLineElement('''EROS C3
1 54880U 22179A   23146.46287103  .00008319  00000+0  27339-3 0  9990
2 54880 139.3580   2.7579 0012797 254.3011 105.6391 15.29433100 22500''')
eros = Satellite(tle)

jd = now()

from random import seed, uniform
seed()
testingGeo = 0
def testGeo(sat):
    global testingGeo
    geo = GeoPosition(uniform(-90, 90), uniform(-180, 180))
    testingGeo = geo
    eVecs = _getEllipseVectors(sat, jd)
    k = _getConstants(*eVecs, geo, jd)
    zk = geo.getZenithVector(jd)[2]
    rho = cos(radians(geo.latitude))
    func = ZeroFunction(zk, rho, k)
    times = func.computeGeneralZeros(1, k) + func.computeGeneralZeros(-1, k)
    return geo, times

errors = []
results = []
def longTest(sat, count):
    global results
    results = []
    for i in range(count):
        try:
            geo, times = testGeo(sat)
            results.append((geo, times))
        except Exception as e:
            errors.append((testingGeo, e))


def getFunc(sat, geo):
    eVecs = _getEllipseVectors(sat, jd)
    k = _getConstants(*eVecs, geo, jd)
    zk = geo.getZenithVector(jd)[2]
    rho = cos(radians(geo.latitude))
    func = ZeroFunction(zk, rho, k)
    return func

def _getZj(zi: float, rho: float, sign: int) -> float:
    assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

    return sign * sqrt(rho * rho - zi * zi)


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


def _getZiPrime(zi: float, zj: float, zk: float, k: list[float]) -> float:
    zjP = -zi / zj

    term1 = 2 * (k[0] * zi + k[1] * zj * zjP)
    term2 = k[3] * (zi * zjP + zj)
    term3 = k[4] * zk
    term4 = k[5] * zk * zjP

    return term1 + term2 + term3 + term4 + k[6] + k[7] * zjP


def _getZiPPrime(zi: float, zj: float, zk: float, k: list[float]) -> float:
    zjP = -zi / zj
    zjPP = (zjP * zi - zj) / (zj * zj)

    term1 = 2 * (k[0] + (k[1] * (zj * zjPP + zjP * zjP)))
    term2 = k[3] * (zi * zjPP + 2 * zjP)
    term3 = k[5] * zk * zjPP
    term4 = k[7] * zjPP

    return term1 + term2 + term3 + term4


def _getZiPPPrime(zi: float, zj: float, zk: float, k: list[float]) -> float:
    zjP = -zi / zj
    zjPP = (zjP * zi - zj) / (zj * zj)
    zjPPP = ((zj*zj) * (zjP + zjPP * zi - zjP) - 2 * zj * zjP * (zjP * zi - zj)) / (zj*zj*zj*zj)
    return (2 * k[1]) * (zj * zjPPP + zjP * zjPP + 2 * zjP * zjPP) + k[3] * (zi * zjPPP + 3 * zjPP) + k[
        5] * zk * zjPPP + k[7] * zjPPP


def _getEllipseVectors(sat: Orbitable, jd: JulianDate) -> list[Vector]:
    elements = sat.getElements(jd)
    matrix = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    aMag = elements.sma
    bMag = aMag * sqrt(1 - elements.ecc*elements.ecc)
    cMag = aMag * elements.ecc

    a = rotateMatrixFrom(matrix, Vector(aMag, 0, 0))
    b = rotateMatrixFrom(matrix, Vector(0, bMag, 0))
    c = rotateMatrixFrom(matrix, Vector(-cMag, 0, 0))

    return a, b, c


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


class Point:
    __slots__ = '_zi', '_value'

    def __init__(self, zi: float, value: float):
        if not isinstance(zi, float):
            raise TypeError(f'zi must be a float type, not {type(zi)}')
        if not isinstance(value, float):
            raise TypeError(f'value must be a float type, not {type(value)}')

        self._zi = zi
        self._value = value

    @classmethod
    def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zj = _getZj(zi, func.rho, sign)
        value = _getZiValue(zi, zj, func.zk, k)

        rtn = object.__new__(cls)
        rtn.__init__(zi, value)

        return rtn

    def __str__(self) -> str:
        return str((self._zi, self._value))

    def __repr__(self) -> str:
        return f'Point({self.__str__()})'

    @property
    def zi(self):
        return self._zi

    @property
    def value(self):
        return self._value


EXTREMA_TAIL = 0
EXTREMA_MINMAX = 1


class Extrema(Point):
    __slots__ = '_type'

    def __init__(self, zi: float, value: float, xType: int):
        super().__init__(zi, value)

        if not isinstance(xType, int):
            raise TypeError(f'xType must be an int type, not {type(xType)}')

        self._type = xType

    @classmethod
    def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zj = _getZj(zi, func.rho, sign)
        value = _getZiValue(zi, zj, func.zk, k)

        rtn = object.__new__(cls)
        rtn.__init__(zi, value, )

        return rtn

    def __str__(self) -> str:
        return str((self._zi, self._value, self._type))

    def __repr__(self) -> str:
        return f'Extrema({self.__str__()})'

    @property
    def type(self):
        return self._type


BOUND_EXTREMA = 2
BOUND_OTHER = 3


class Bound:
    __slots__ = '_upper', '_lower', '_left', '_right'   # , '_upperClass', '_lowerClass'

    def __init__(self, point1: Point, point2: Point):
        # this doesn't need to be checked with the Domain.__init__() method
        # assert point1.value * point2.value < 0, 'bound values are both positive or both negative'

        if point1.value > 0:
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

        # if isinstance(self._upper, Extrema):
        #     self._upperClass = BOUND_EXTREMA
        # else:
        #     if isinstance(self._upper, Point):
        #         self._upperClass = BOUND_OTHER
        #     else: #elif not isinstance(left, Point):
        #         raise TypeError(f'left must be Point type, not {type(self._upper)}')
        # if isinstance(self._lower, Extrema):
        #     self._lowerClass = BOUND_EXTREMA
        # else:
        #     if isinstance(self._lower, Point):
        #         self._lowerClass = BOUND_OTHER
        #     else: #if not isinstance(right, Point):
        #         raise TypeError(f'right must be Point type, not {type(self._lower)}')
        if not isinstance(point1, Point):
            raise TypeError(f'right must be Point type, not {type(point1)}')
        if not isinstance(point2, Point):
            raise TypeError(f'right must be Point type, not {type(point2)}')

    def __str__(self) -> str:
        return str((self._upper, self._lower))

    def __repr__(self) -> str:
        return f'Bound({self.__str__()})'

    def __contains__(self, item: float):
        if not isinstance(item, float):
            raise TypeError(f'item must be float type, not {type(item)}')

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

    # @property
    # def leftClass(self):
    #     return self._leftClass
    #
    # @property
    # def rightClass(self):
    #     return self._rightClass


class Domain(Bound):

    def __init__(self, left: Point, right: Point):
        super().__init__(left, right)

    @classmethod
    def fromFunction(cls, func: 'ZeroFunction', sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = ALMOST_ONE * func.rho
        left = Point.fromZi(-limit, func, sign, k)
        right = Point.fromZi(limit, func, sign, k)

        rtn = object.__new__(cls)
        rtn.__init__(left, right)

        return rtn


ALMOST_ONE = 0.9999999
EPSILON = 1e-7
NEWTON_DERIVATIVE_EPSILON = 1e-3
NEWTON_EPSILON = 0.1
NEWTON_GAP = 1e-7
EXTREMA_DERIVATIVE_EPSILON = 0.1
ZERO_FUNCTION = 0
DERIVATIVE_FUNCTION = 1
INITIAL_REDUCTION_COUNT = 4


def _checkFloat(value, name: str):
    if not isinstance(value, float):
        raise TypeError(f'{name} must be a float type, not {type(value)}')


def _checkArray(array, name):
    if not isinstance(array, (list, tuple)):
        raise TypeError(f'{name} must be a list or tuple, not {type(array)}')
    length = len(array)
    if length != 10:
        raise ValueError(f'{name} must have exactly 10 elements, {length} found')


class ZeroFunction:
    __slots__ = '_zk', '_rho', '_k'

    def __init__(self, zk: float, rho: float, k: list[float]):

        _checkFloat(zk, 'zk')
        _checkFloat(rho, 'rho')
        _checkArray(k, 'k')

        self._zk = zk
        self._rho = rho
        self._k = k

    # todo: this might not be used for the derivative, should remove function arg if that's the case
    def _halfBound(self, bound: Bound, sign: int, k: list[float], function) -> (Bound, float):
    # def _halfBound(self, bound: Bound, sign: int, k: list[float], function: int):
    #     assert function == ZERO_FUNCTION or function == DERIVATIVE_FUNCTION, 'function must be ZERO_FUNCTION or ' \
    #                                                                          'DERIVATIVE_FUNCTION '

        middle = (bound.upper.zi + bound.lower.zi) / 2
        value = function(middle, _getZj(middle, self._rho, sign), self._zk, k)
        # if function == ZERO_FUNCTION:
        #     value = _getZiValue(middle, _getZj(middle, self._rho, sign), self._zk, k)
        # else:
        #     value = _getZiPrime(middle, _getZj(middle, self._rho, sign), self._zk, k)

        point = Point(middle, value)
        if value >= 0:
            rtn = Bound(point, bound.lower)
        else:
            rtn = Bound(point, bound.upper)

        newZi = (rtn.upper.zi + rtn.lower.zi) / 2

        return rtn, newZi

    def _getThirdZero(self, sign: int, k: list[float]) -> Extrema:
        """compute the zeros of the third derivative"""
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = self._rho * ALMOST_ONE
        left = Point(-limit, _getZiPPPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k))
        right = Point(limit, _getZiPPPrime(limit, _getZj(limit, self._rho, sign), self._zk, k))
        bound = Bound(left, right)

        zi = 0
        while abs(bound.left.zi - bound.right.zi) > EPSILON:
            bound, zi = self._halfBound(bound, sign, k, _getZiPPPrime)

        return Extrema(zi, _getZiPPrime(zi, _getZj(zi, self._rho, sign), self._zk, k), EXTREMA_MINMAX)

    def _getSecondZero(self, sign: int, k: list[float]) -> Extrema:
        """compute the zeros of the second derivative"""
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'
        
        limit = self._rho * ALMOST_ONE
        left = Point(-limit, _getZiPPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k))
        right = Point(limit, _getZiPPrime(limit, _getZj(limit, self._rho, sign), self._zk, k))
        bound = Bound(left, right)

        zi = 0
        while abs(bound.left.zi - bound.right.zi) > EPSILON:
            bound, zi = self._halfBound(bound, sign, k, _getZiPPrime)

        return Extrema(zi, _getZiPrime(zi, _getZj(zi, self._rho, sign), self._zk, k), EXTREMA_MINMAX)

    def _getOddExtremaCount(self, extrema: Extrema, sign: int, k: list[float]) -> int:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = self._rho * ALMOST_ONE
        leftValue = _getZiPPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k)
        if __debug__:
            rightValue = _getZiPPrime(limit, _getZj(limit, self._rho, sign), self._zk, k)
        assert leftValue * rightValue > 1

        if leftValue * extrema.value > 0:
            return 1
        else:
            return 3

    def _reduceBound(self, bound: Bound, sign: int, k: list[float], function, derivative):
        """reduces a bound to a single value using newton's method and checking the next iteration stays within the
        valid bound, bifurcating the bound each time to ensure the algorithm converges
        """
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zi1 = 0
        for i in range(INITIAL_REDUCTION_COUNT):
            bound, zi1 = self._halfBound(bound, sign, k, function)

        zi0 = -2
        zj = _getZj(zi1, self._rho, sign)
        while abs(function(zi1, zj, self._zk, k)) > NEWTON_EPSILON and abs(zi1 - zi0) > NEWTON_GAP:
            if zi1 not in bound:
                bound, zi1 = self._halfBound(bound, sign, k, function)
            else:
                zi0 = zi1
                zj = _getZj(zi0, self._rho, sign)
                zi1 = zi0 - function(zi0, zj, self._zk, k) / derivative(zi0, zj, self._zk, k)

        return zi1

    def _getEvenExtrema(self, sign: int, k: list[float]):
        """get the extrema for a zero function whose first derivative tails have the same sign"""
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get the zero from the 2nd derivative and make an extrema value
        limit = self._rho * ALMOST_ONE
        firstExtrema = self._getSecondZero(sign, k)

        # use the tails and the extrema to make left and right bounds of 1st derivative
        leftSecond = Point(-limit, _getZiPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k))
        rightSecond = Point(limit, _getZiPrime(limit, _getZj(limit, self._rho, sign), self._zk, k))
        leftBound = Bound(leftSecond, firstExtrema)
        rightBound = Bound(firstExtrema, rightSecond)

        # use the bounds to find the two zeros, and the extrema they represent
        leftExtremaZi = self._reduceBound(leftBound, sign, k, _getZiPrime, _getZiPPrime)
        rightExtremaZi = self._reduceBound(rightBound, sign, k, _getZiPrime, _getZiPPrime)
        leftExtremaValue = _getZiValue(leftExtremaZi, _getZj(leftExtremaZi, self._rho, sign), self._zk, k)
        rightExtremaValue = _getZiValue(rightExtremaZi, _getZj(rightExtremaZi, self._rho, sign), self._zk, k)

        return Extrema(leftExtremaZi, leftExtremaValue, EXTREMA_MINMAX), Extrema(rightExtremaZi, rightExtremaValue, EXTREMA_MINMAX)


    def _getOddExtrema(self, sign: int, k: list[float]):
        """get the extrema for a zero function whose first derivative tails have opposite signs"""
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get zero from the 3rd derivative as an extrema value for the 2nd derivative
        secondExtrema = self._getThirdZero(sign, k)

        # compare the signs and the 2nd derivative tails and the extrema
        limit = self._rho * ALMOST_ONE
        leftDDTail = Extrema(-limit, _getZiPPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k), EXTREMA_TAIL)
        rightDDTail = Extrema(limit, _getZiPPrime(limit, _getZj(limit, self._rho, sign), self._zk, k), EXTREMA_TAIL)
        leftDTail = Extrema(-limit, _getZiPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k), EXTREMA_TAIL)
        rightDTail = Extrema(limit, _getZiPrime(limit, _getZj(limit, self._rho, sign), self._zk, k), EXTREMA_TAIL)
        assert leftDDTail.value * rightDDTail.value > 0
        # if they are the same:
        if leftDDTail.value * secondExtrema.value > 0:
            # find the zero of the 1st derivative as an extrema value for the zero function
            zi = self._reduceBound(Bound(leftDTail, rightDTail), sign, k, _getZiPrime, _getZiPPrime)
            return Extrema(zi, _getZiValue(zi, _getZj(zi, self._rho, sign), self._zk, k), EXTREMA_MINMAX)
        # if they are difference:
        else:
            # find the zeros to the 2nd derivative as extrema for the 1st derivative
            leftDDBound = Bound(leftDDTail, secondExtrema)
            rightDDBound = Bound(secondExtrema, rightDDTail)
            leftDExtremaZi = self._reduceBound(leftDDBound, sign, k, _getZiPPrime, _getZiPPPrime)
            rightDExtremaZi = self._reduceBound(rightDDBound, sign, k, _getZiPPrime, _getZiPPPrime)
            leftDExtrema = Extrema.fromZi(leftDExtremaZi, self, sign, k)
            rightDExtrema = Extrema.fromZi(rightDExtremaZi, self, sign, k)
            # use the tails and extrema of the 1st derivative to get zeros as extrema values for the zero function
            leftDBound = Bound(leftDTail, leftDExtrema)
            middleDBound = Bound(leftDExtrema, rightDExtrema)
            rightDBound = Bound(rightDExtrema, rightDTail)
            leftExtrema = self._reduceBound(leftDBound, sign, k, _getZiPrime, _getZiPPrime)
            middleExtrema = self._reduceBound(middleDBound, sign, k, _getZiPrime, _getZiPPrime)
            rightExtrema = self._reduceBound(rightDBound, sign, k, _getZiPrime, _getZiPPrime)
            return leftExtrema, middleExtrema, rightExtrema

    def _computeExtrema(self, sign: int, k: list[float]) -> list[Extrema]:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = self._rho * ALMOST_ONE
        leftRate = _getZiPrime(-limit, _getZj(-limit, self._rho, sign), self._zk, k)
        rightRate = _getZiPrime(limit, _getZj(limit, self._rho, sign), self._zk, k)

        if leftRate * rightRate > 0:
            return self._getEvenExtrema(sign, k)
        else:
            return self._getOddExtrema(sign, k)

    START_LEFT = -1
    START_ZERO = 0
    START_RIGHT = 1

    def _getPrimeStartingZi(self, start: int, sign: int, k: list[float]):
        assert start in (self.START_LEFT, self.START_ZERO, self.START_RIGHT), 'start must be START_LEFT, START_ZERO ' \
                                                                              'or START_RIGHT '
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        if start == self.START_ZERO:
            return 0
        else:
            zi1 = start * self._rho * ALMOST_ONE
            zj = _getZj(zi1, self._rho, sign)
            ziP = _getZiPrime(zi1, zj, self._zk, k)
            ziPP = _getZiPPrime(zi1, zj, self._zk, k)

            if start == self.START_LEFT:
                while (ziP > 0 and ziPP > 0) or (ziP < 0 and ziPP < 0):
                    zi1 *= 0.99
                    zj = _getZj(zi1, self._rho, sign)
                    ziP = _getZiPrime(zi1, zj, self._zk, k)
                    ziPP = _getZiPPrime(zi1, zj, self._zk, k)
            elif start == self.START_RIGHT:
                while (ziP > 0 and ziPP < 0) or (ziP < 0 and ziPP > 0):
                    zi1 *= 0.99
                    zj = _getZj(zi1, self._rho, sign)
                    ziP = _getZiPrime(zi1, zj, self._zk, k)
                    ziPP = _getZiPPrime(zi1, zj, self._zk, k)

        return zi1

    NO_EXTREMA_FOUND = 2

    def _getPrimeZero(self, start: int, sign: int, k: list[float]):
        assert start in (self.START_LEFT, self.START_ZERO, self.START_RIGHT), 'start must be START_LEFT, START_ZERO ' \
                                                                              'or START_RIGHT '
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zi0 = 2
        zi1 = self._getPrimeStartingZi(start, sign, k)
        zj = _getZj(zi1, self._rho, sign)
        # while abs(_getZiPrime(zi1, zj, self._zk, k)) > NEWTON_DERIVATIVE_EPSILON :
        while abs(_getZiPrime(zi1, zj, self._zk, k)) > NEWTON_EPSILON and abs(zi1 - zi0) > NEWTON_GAP:
            zi0 = zi1
            try:
                # don't call _getZj because of the assert statement, want to catch the error here
                zj = sign * sqrt(self._rho*self._rho - zi0*zi0)
            except ValueError:
                # todo: reraise a custom exception here?
                return self.NO_EXTREMA_FOUND
            zi1 = zi0 - _getZiPrime(zi0, zj, self._zk, k) / _getZiPPrime(zi0, zj, self._zk, k)

        return zi1

    def _computeExtremaValues(self, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        primeZeros = [self._getPrimeZero(s, sign, k) for s in (self.START_LEFT, self.START_ZERO, self.START_RIGHT)]
        args = [(z, _getZj(z, self._rho, sign)) for z in primeZeros if z != self.NO_EXTREMA_FOUND]
        extrema = [Extrema(z, _getZiValue(z, zj, self._zk, k), EXTREMA_MINMAX) for z, zj in args]

        return extrema

    def _getTails(self, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        domain = Domain.fromFunction(self, sign, k)

        zjRight = _getZj(domain.right.zi, self._rho, sign)
        if abs(_getZiPrime(domain.right.zi, zjRight, self._zk, k)) < EXTREMA_DERIVATIVE_EPSILON:
            rightType = EXTREMA_MINMAX
        else:
            rightType = EXTREMA_TAIL

        zjLeft = _getZj(domain.left.zi, self._rho, sign)
        if abs(_getZiPrime(domain.left.zi, zjLeft, self._zk, k)) < EXTREMA_DERIVATIVE_EPSILON:
            leftType = EXTREMA_MINMAX
        else:
            leftType = EXTREMA_TAIL

        rightTail = Extrema(domain.right.zi, domain.right.value, rightType)
        leftTail = Extrema(domain.left.zi, domain.left.value, leftType)

        return rightTail, leftTail

    def _getStartingZi(self, bound: Bound) -> float:
        pass

    def _getZiSafely(self, bound: Bound, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zi1 = 0
        for i in range(INITIAL_REDUCTION_COUNT):
            bound, zi1 = self._halfBound(bound, sign, k, ZERO_FUNCTION)

        zi0 = 2
        zj = _getZj(zi1, self._rho, sign)
        while abs(_getZiValue(zi1, zj, self._zk, k)) > NEWTON_EPSILON and abs(zi0 - zi1) > NEWTON_GAP:
            if zi1 not in bound:
                bound, zi1 = self._halfBound(bound, sign, k, ZERO_FUNCTION)
            else:
                zj = _getZj(zi1, self._rho, sign)
                zi0 = zi1
                zi1 = zi0 - _getZiValue(zi0, zj, self._zk, k) / _getZiPrime(zi0, zj, self._zk, k)

        return zi1

    def _computeBounds(self, sign: int, k: list[float]) -> list[Bound]:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        extrema = list(self._getTails(sign, k)) + self._computeExtremaValues(sign, k)
        sortedExtrema = sorted(extrema, key=lambda o: o.zi)

        bounds = []
        for i, x, in enumerate(sortedExtrema[:-1]):
            if x.value * sortedExtrema[i + 1].value < 0:
                bound = Bound(x, sortedExtrema[i + 1])
                bounds.append(bound)

        return bounds

    def computeGeneralZeros(self, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        bounds = self._computeBounds(sign, k)

        intersections = [self._getZiSafely(bound, sign, k) for bound in bounds]

        return intersections

    def computeSpecificZero(self, zi: float, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        bounds = self._computeBounds(sign, k)
        bound = None
        for b in bounds:
            if zi in b:
                bound = b

        # assert bound != None, 'zi not in any bounds'
        # if the zero switched signs, the bound it was in will no longer be returned by _computeBounds,
        #   since the value of each extrema will have the same sign. here's where we need to raise a
        #   specific exception to tell the calling method the zi should be found in the function with
        #   the current k array but with flipped sign
        if bound is None:
            raise ValueError('zi disappeared, check the other sign')

        # another solution coule be call _getZiSafely with sign arg as -sign, but that doesn't tell the
        # calling function next
        return self._getZiSafely(bound, sign, k)

    def plot(self, sign, k):
        if __debug__:
            import matplotlib.pyplot as plt

            limit = self._rho * ALMOST_ONE
            plotCount = 1000
            m = int(limit * plotCount)
            x = [-limit] + [i / plotCount for i in range(-m, m)] + [limit]
            f, fPrime, fPPrime, fPPPrime = [], [], [], []

            for zi in x:
                zj = _getZj(zi, self._rho, sign)
                f.append(_getZiValue(zi, zj, self._zk, k))
                fPrime.append(_getZiPrime(zi, zj, self._zk, k))
                fPPrime.append(_getZiPPrime(zi, zj, self._zk, k))
                fPPPrime.append(_getZiPPPrime(zi, zj, self._zk, k))

            figure = plt.figure()
            ax = figure.add_subplot(111)
            ax.plot(x, f, '-', c='blue')
            ax.plot(x, fPrime, '-', c='dodgerblue')
            ax.plot(x, fPPrime, '-', c='powderblue')
            ax.plot(x, fPPPrime, '-', c='lightcyan')

            maxVal = max(f) * 1.2
            minVal = min(f) * 1.2
            ax.set_ylim((minVal, maxVal))
            ax.grid()

            plt.show()
        else:
            raise NotImplemented('matplotlib is not a production dependency')

    @property
    def rho(self):
        return self._rho

    @property
    def zk(self):
        return self._zk


SMALL = 1e-7


class OrbitPath:
    __slots__ = '_sat', '_geo', '_jd', '_zk', '_rho'

    def __init__(self, sat: Orbitable, geo: GeoPosition, jd: JulianDate):
        if not isinstance(sat, OrbitPath):
            raise TypeError(f'sat must be Orbitable type, not {type(sat)}')
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be GeoPosition type, not {type(geo)}')
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be JulianDate type, not{type(jd)}')

        self._sat = sat
        self._geo = geo
        self._jd = jd

        self._zk = geo.getZenithVector(jd)
        self._rho = cos(radians(geo.latitude))

    def _getGeneralZeros(self, jd):
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)

        func = ZeroFunction(self._zk, self._rho, k)
        zeros = []
