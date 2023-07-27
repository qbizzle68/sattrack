from copy import copy
from math import sqrt, radians, cos, sin, pi, atan2

from _pyevspace import Vector, getMatrixEuler, ZXZ, Angles, rotateMatrixFrom, dot

from sattrack.coordinates import GeoPosition
from sattrack.orbit import Orbitable
from sattrack.spacetime.juliandate import JulianDate

if __debug__ is True:
    debug__all__ = ['_getZj', '_getZjPrime', '_getZjPPrime', '_getZjPPPrime', '_getZiValue', '_getZiPrime',
                    '_getZiPPrime', '_getZiPPPrime', '_getEllipseVectors', '_getConstants', '_angleDifference',
                    '_signOf']
else:
    debug__all__ = []

__all__ = ['getGeoValues', 'getEllipseVectors', 'getConstants', 'Point', 'Extrema', 'Boundary', 'ALMOST_ONE',
           'DOMAIN_ONE', 'Domain', 'DUPLICATE_ZERO_EPSILON', 'ZeroIntersection', 'ZeroDisappeared', 'EPSILON',
           'INITIAL_REDUCTION_COUNT', 'NEWTON_EPSILON', 'NEWTON_GAP', 'ZeroFunction', 'PREVIOUS_OCCURRENCE',
           'NEXT_OCCURRENCE', 'SMALL', 'CHECK_MINIMUM', 'TIME_DIFFERENCE', 'OrbitPath'] + debug__all__

from sattrack.util.constants import SIDEREAL_PER_SOLAR

# for developing
from sattrack import *

tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23147.70061491  .00011858  00000+0  21377-3 0  9990
2 25544  51.6408  70.0897 0005510  26.3638  52.3806 15.50214517398618''')
iss = Satellite(tle)
jd = now()
geo = GeoPosition(38, -98)
__all__ += ['iss', 'jd', 'geo']


# disappearing zero happens with !!!!!!!!!!!!!!!!!!!!1
# tle = TwoLineElement('''ISS (ZARYA)
# 1 25544U 98067A   23147.70061491  .00011858  00000+0  21377-3 0  9990
# 2 25544  51.6408  70.0897 0005510  26.3638  52.3806 15.50214517398618''')
# jd = JulianDate(6, 7, 2023, 10, 3, 6.379000000000815, -5.0)
# geo = GeoPosition(-50.5146337828902, -28.027297912220664, 0)
# jd2 = jd.future(0.1)

def getGeoValues(sat, geo, jd):
    zk = sin(radians(geo.latitude))
    rho = cos(radians(geo.latitude))
    eVecs = _getEllipseVectors(sat, jd)
    k = _getConstants(*eVecs, geo, jd)

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


class Point:
    __slots__ = '_zi', '_value'

    def __init__(self, zi: float, value: float):
        _checkFloat(zi, 'zi')
        _checkFloat(value, 'value')

        self._zi = zi
        self._value = value

    @classmethod
    def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, k: list[float], method=_getZiValue):
        assert sign == 1 or sign == -1, f'sign expected to be +/- 1, was {sign}'

        zj = _getZj(zi, func.rho, sign)
        value = method(zi, zj, func.zk, k)
        # value = _getZiValue(zi, zj, func.zk, k)

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


class Extrema(Point):

    def __repr__(self) -> str:
        return f'Extrema({self.__str__()})'


class Boundary:
    __slots__ = '_upper', '_lower', '_left', '_right'

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

    def __str__(self) -> str:
        return str((self._upper, self._lower))

    def __repr__(self) -> str:
        return f'Bound({self.__str__()})'

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


ALMOST_ONE = 0.9999999
DOMAIN_ONE = 0.999999999999999


class Domain(Boundary):

    def __init__(self, point1: Point, point2: Point):
        super().__init__(point1, point2)

    @classmethod
    def fromZeroFunction(cls, func: 'ZeroFunction', sign: int, k: list[float], method=_getZiValue):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = func.rho * DOMAIN_ONE
        leftValue = method(-limit, _getZj(-limit, func.rho, sign), func.zk, k)
        rightValue = method(limit, _getZj(limit, func.rho, sign), func.zk, k)
        left = Point(-limit, leftValue)
        right = Point(limit, rightValue)

        rtn = object.__new__(cls)
        rtn.__init__(left, right)

        return rtn


DUPLICATE_ZERO_EPSILON = 1e-3
ZERO_ASCENDING = 1
ZERO_DESCENDING = -1


class ZeroIntersection:
    __slots__ = '_zi', '_sign', '_direction'

    def __init__(self, zi: float, sign: int, ziP: float):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        _checkFloat(zi, 'zi')
        _checkFloat(ziP, 'ziP')

        if ziP > 0:
            self._direction = ZERO_ASCENDING
        else:
            self._direction = ZERO_DESCENDING

        self._zi = zi
        self._sign = sign

    def __str__(self) -> str:
        return str((self._zi, self._sign))

    def __repr__(self) -> str:
        return f'ZeroIntersection({self.__str__()})'

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
        return float(self._direction)

    @direction.setter
    def direction(self, value: int):
        assert value == ZERO_ASCENDING or value == ZERO_DESCENDING, f'value must be ZERO_ASCENDING or ' \
                                                                    f'ZERO_DESCENDING, not {value} '
        self._direction = value


class ZeroDisappeared(Exception):
    def __init__(self, message):
        super().__init__(message)


EPSILON = 1e-7
INITIAL_REDUCTION_COUNT = 4
NEWTON_EPSILON = 0.1
NEWTON_GAP = 1e-7


class ZeroFunction:
    __slots__ = '_zk', '_rho', '_geo', '_zeros'

    # def __new__(cls, geo: GeoPosition, jd: JulianDate):
    #     if not isinstance(geo, GeoPosition):
    #         raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
    #     if not isinstance(jd, JulianDate):
    #         raise TypeError(f'jd must be a JulianDate type, not {type(jd)}')
    #
    #     rtn = object.__new__(cls)
    #     lat = radians(geo.latitude)
    #     rtn._zk = sin(lat)
    #     rtn._rho = cos(lat)
    #     rtn._geo = geo
    #     rtn._zeros = []
    #
    # def __init__(self, geo: GeoPosition, jd: JulianDate):
    #     if not isinstance(geo, GeoPosition):
    #         raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
    #
    #     lat = radians(geo.latitude)
    #     self._zk = sin(lat)
    #     self._rho = cos(lat)
    #     self._geo = geo
    #     self._zeros = []

    def __init__(self, geo: GeoPosition):
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')

        lat = radians(geo.latitude)
        self._zk = sin(lat)
        self._rho = cos(lat)

    @classmethod
    def fromValues(cls, zk: float, rho: float):
        _checkFloat(zk, 'zk')
        _checkFloat(rho, 'rho')

        rtn = object.__new__(cls)
        rtn._zk = zk
        rtn._rho = rho
        # rtn._geo = None
        # rtn._zeros = []

        return rtn

    def _halfBoundary(self, bound: Boundary, sign: int, k: list[float], function) -> (Boundary, float):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        middle = (bound.upper.zi + bound.lower.zi) / 2
        point = Point.fromZi(middle, self, sign, k, function)

        if point.value >= 0:
            rtn = Boundary(point, bound.lower)
        else:
            rtn = Boundary(point, bound.upper)

        newZi = (rtn.upper.zi + rtn.lower.zi) / 2

        return rtn, newZi

    def _getZeroBifurcate(self, sign: int, k: list[float], function, derivative) -> Extrema:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        limit = self._rho * DOMAIN_ONE
        left = Point.fromZi(-limit, self, sign, k, derivative)
        right = Point.fromZi(limit, self, sign, k, derivative)
        bound = Boundary(left, right)

        zi = 0
        while abs(bound.left.zi - bound.right.zi) > EPSILON:
            bound, zi = self._halfBoundary(bound, sign, k, derivative)

        return Extrema.fromZi(zi, self, sign, k, function)

    def _reduceBoundary(self, bound: Boundary, sign: int, k: list[float], function, derivative):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        zi1 = 0
        for i in range(INITIAL_REDUCTION_COUNT):
            bound, zi1 = self._halfBoundary(bound, sign, k, function)

        zi0 = -2
        zj = _getZj(zi1, self._rho, sign)
        while abs(function(zi1, zj, self._zk, k)) > NEWTON_EPSILON and abs(zi1 - zi0) > NEWTON_GAP:
            if zi1 not in bound:
                bound, zi1 = self._halfBoundary(bound, sign, k, function)
            else:
                zi0 = zi1
                zj = _getZj(zi0, self._rho, sign)
                zi1 = zi0 - function(zi0, zj, self._zk, k) / derivative(zi0, zj, self._zk, k)

        return zi1

    def _getEvenExtrema(self, sign: int, k: list[float]):
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get the zero from the 2nd derivative and make an extrema value
        firstExtrema = self._getZeroBifurcate(sign, k, _getZiPrime, _getZiPPrime)

        # use the tails and the extrema to make left and right bounds of 1st derivative
        domain = Domain.fromZeroFunction(self, sign, k, _getZiPrime)
        leftBound = Boundary(domain.left, firstExtrema)
        rightBound = Boundary(firstExtrema, domain.right)

        # use the bounds to find the two zeros, and the extrema they represent
        leftExtremaZi = self._reduceBoundary(leftBound, sign, k, _getZiPrime, _getZiPPrime)
        rightExtremaZi = self._reduceBoundary(rightBound, sign, k, _getZiPrime, _getZiPPrime)
        leftExtrema = Extrema.fromZi(leftExtremaZi, self, sign, k, _getZiValue)
        rightExtrema = Extrema.fromZi(rightExtremaZi, self, sign, k, _getZiValue)

        return [leftExtrema, rightExtrema]

    def _getOddExtrema(self, sign: int, k: list[float]) -> tuple[Extrema]:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # get zero from the 3rd derivative as an extrema value for the 2nd derivative
        secondExtrema = self._getZeroBifurcate(sign, k, _getZiPPrime, _getZiPPPrime)

        # compare the signs and the 2nd derivative tails and the extrema
        domainDD = Domain.fromZeroFunction(self, sign, k, _getZiPPrime)
        assert domainDD.left.value * domainDD.right.value > 0
        domainD = Domain.fromZeroFunction(self, sign, k, _getZiPrime)

        # if they are the same:
        if domainDD.left.value * secondExtrema.value > 0:
            # find the zero of the 1st derivative as an extrema value for the zero function
            zi = self._reduceBoundary(Boundary(domainD.left, domainD.right), sign, k, _getZiPrime, _getZiPPrime)
            return [Extrema.fromZi(zi, self, sign, k, _getZiValue)]
        # if they are different:
        else:
            # find the zeros to the 2nd derivative as extrema for the 1st derivative
            leftDDBound = Boundary(domainDD.left, secondExtrema)
            rightDDBound = Boundary(secondExtrema, domainDD.right)
            leftDExtremaZi = self._reduceBoundary(leftDDBound, sign, k, _getZiPPrime, _getZiPPPrime)
            rightDExtremaZi = self._reduceBoundary(rightDDBound, sign, k, _getZiPPrime, _getZiPPPrime)
            leftDExtrema = Extrema.fromZi(leftDExtremaZi, self, sign, k, _getZiPrime)
            rightDExtrema = Extrema.fromZi(rightDExtremaZi, self, sign, k, _getZiPrime)

            # use tails and extrema of the 1st derivative to check for zeros
            bounds = [Boundary(domainD.left, leftDExtrema), Boundary(leftDExtrema, rightDExtrema),
                      Boundary(rightDExtrema, domainD.right)]
            existingZeros = [bound for bound in bounds if bound.left.value * bound.right.value < 0]

            # use the tails and valid extrema of the 1st derivative to get zeros as extrema values for the zero function
            extremaZis = [self._reduceBoundary(bound, sign, k, _getZiPrime, _getZiPPrime) for bound in existingZeros]
            return [Extrema.fromZi(zi, self, sign, k, _getZiValue) for zi in extremaZis]
            # rtn = tuple(Extrema.fromZi(zi, self, sign, k, _getZiValue) for zi in extremaZis)
            # if len(rtn) == 1:
            #     return rtn[0]
            # return rtn

    def _computeExtrema(self, sign: int, k: list[float]) -> list[Extrema]:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        domainD = Domain.fromZeroFunction(self, sign, k, _getZiPrime)

        if domainD.left.value * domainD.right.value > 0:
            return self._getEvenExtrema(sign, k)
        else:
            return self._getOddExtrema(sign, k)

    def _computeZeros(self, sign: int, k: list[float]) -> list[ZeroIntersection]:
        assert sign == 1 or sign == -1, f'plus expected to be +/- 1, was {sign}'

        # getExtrema of the zero function
        extrema = self._computeExtrema(sign, k)

        # get tails as extrema from a Domain
        domain = Domain.fromZeroFunction(self, sign, k, _getZiValue)

        # get bounds for every tail/extrema
        if not extrema:
            bounds = [domain]
        else:
            bounds = [Boundary(domain.left, extrema[0])] + [Boundary(extrema[i], extrema[i + 1]) for i in
                                                            range(len(extrema) - 1)] + [
                         Boundary(extrema[-1], domain.right)]

        # find where the actual zeros exist
        existingZeros = [bound for bound in bounds if bound.left.value * bound.right.value < 0]

        # find and return the zeros
        zeros = [self._reduceBoundary(bound, sign, k, _getZiValue, _getZiPrime) for bound in existingZeros]
        aux = [_getZiPrime(z, _getZj(z, self._rho, sign), self._zk, k) for z in zeros]
        return [ZeroIntersection(zi, sign, ziP) for zi, ziP in zip(zeros, aux)]

    def computeGeneralZeros(self, k: list[float]) -> list[ZeroIntersection]:
        # Zeros that occur in both positive and negative functions refer to different geo-positions
        #   (they have difference zjs' (one negative and one positive) so they are not duplicates).
        return self._computeZeros(1, k) + self._computeZeros(-1, k)

    # this method will have difficulties recognizing the same zero as time increased between k arrays
    def computeSpecificZero(self, zero: ZeroIntersection, k: list[float]) -> ZeroIntersection:
        if not isinstance(zero, ZeroIntersection):
            raise TypeError(f'zero must be a ZeroIntersection type, not {type(zero)}')

        zeros = self._computeZeros(zero.sign, k)
        zCount = len(zeros)
        if zCount == 0:
            raise ZeroDisappeared(f'zero: {zero.zi} not found in function with sign: {zero.sign}')
        elif zCount == 1:
            z = zeros[0]
        else:
            dz = [abs(z.zi - zero.zi) for z in zeros]
            z = zeros[dz.index(min(dz))]

        # the closest zero might not be the same one
        # allowable difference should be dependent on the size of the domain
        ALLOWABLE_DIFFERENCE = self._rho * .1
        if abs(z.zi - zero.zi) < ALLOWABLE_DIFFERENCE:
            return z
        raise ZeroDisappeared(f'zero: {zero.zi} not found in function with sign: {zero.sign}')

    @property
    def rho(self) -> float:
        return self._rho

    @property
    def zk(self) -> float:
        return self._zk

    @property
    def zeros(self):
        return self._zeros

    def plot(self, k: list[float]):
        if __debug__:
            import matplotlib.pyplot as plt

            limit = self._rho * DOMAIN_ONE
            plotCount = 1000
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
    # The nearest absolute distance to any two angles should be less than pi, so if the measured difference
    #   is < -pi, add TWOPI to get the value between [0-pi]. Same reasoning for the difference > pi.

    # An angle can be subtracted and added many times without range checking, so the angle needs to be
    #   'modded' by TWOPI (or compute how many multiples of TWOPI to add/subtract by hand).
    validAngle = angle % TWOPI
    if validAngle > pi:
        return validAngle - TWOPI
    return validAngle


def _signOf(value: float) -> int:
    if value >= 0:
        return 1
    return -1


PREVIOUS_OCCURRENCE = 0
NEXT_OCCURRENCE = 1
SMALL = 1e-7
CHECK_MINIMUM = 0.9999
TIME_DIFFERENCE = 1e-5  # less than 1/10 of a second


class OrbitPath:
    __slots__ = '_sat', '_geo', '_jd', '_zk', '_rho'

    def __init__(self, sat: Orbitable, geo: GeoPosition, jd: JulianDate):
        if not isinstance(sat, Orbitable):
            raise TypeError(f'sat must be an Orbitable type, not {type(sat)}')
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a JulianDate type, not {type(geo)}')

        self._sat = sat
        self._geo = geo
        self._jd = jd

        angle = radians(geo.latitude)
        self._zk = sin(angle)
        self._rho = cos(angle)

    def _getGeneralZeros(self, jd: JulianDate) -> list[ZeroIntersection]:
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)

        func = ZeroFunction.fromValues(self._zk, self._rho)

        # return [ZeroIntersection(zi, 1) for zi in func.computeGeneralZeros(k)]
        return func.computeGeneralZeros(k)

    def _checkTime(self, jd: JulianDate) -> float:
        a, b, c = _getEllipseVectors(self._sat, jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        numerator = dot(zeta, gamma - c)
        aDotZ = dot(zeta, a)
        bDotZ = dot(zeta, b)
        denominator = sqrt(aDotZ*aDotZ + bDotZ*bDotZ)

        return numerator / denominator

    def _orderZero(self, zeros: list[ZeroIntersection], dls: list[float], zCount: int, jd: JulianDate) -> (
                   list[ZeroIntersection], list[bool]):
        assert zCount == 2 or zCount == 4, f'expected zCount of 2 or 4, not {zCount}'

        # Sort zeros and angle differences based on the angle differences
        zeroSorted = [z for _, z in sorted(zip(dls, zeros))]
        dlSorted = sorted(dls)

        # We want to order the zeros to find the next time to the intersecting geo-positions, with
        #   the exception being we want to find the previous zero time if the orbit path is currently
        #   above the geo-position. We only need to sort the zeros accordingly, the sign of the dl is
        #   not significant outside this method, as the direction array specifies which direction to
        #   move in time while searching for intersecting times.
        # todo: should this be in [0, 1] or -[1, 1] ?
        up = -1 < self._checkTime(jd) < 1
        # Find the first positive dl then determine which indices to append to the positive values.
        positiveIndices = [i for i, dl, in enumerate(dlSorted) if dl > 0]
        if positiveIndices:
            firstPositive = positiveIndices[0]
        else:
            # The cases with no positive dls and all positive dls are treated the same, so set firstPositive
            #   to 0 if positiveIndices is empty (see explanation below).
            firstPositive = 0

        # Here is the breakdown for creating the updated zeroSorted arrays, note: fp = firstPositive,
        #   zs = zeroSorted and dls = dlSorted for brevity
        #
        #   dls: [1, 2, 3, 4] [-1, 2, 3, 4] [-1,-2, 3, 4] [-1,-2,-3, 4] [-1,-2,-3,-4]
        #   zs:  [a, b, c, d] [ a, b, c, d] [ a, b, c, d] [ a, b, c, d] [ a, b, c, d]
        #
        #   zs should become the following when up is true
        #   zs:  [d, a, b, c] [ a, b, c, d] [ b, c, d, a] [ c, d, a, b] [ d, a, b, c]
        #   and zs should become the following when up is false
        #   zs:  [a, b, c, d] [ b, c, d, a] [ c, d, a, b] [ d, a, b, c] [ a, b, c, d]
        #
        #   Therefore the correct ordering of zs should be zs[fp-1:length] + zs[0:fp-1] when up is True
        #   and zs[fp:length] + zs[0:fp] when up is False. Notice that the first and last dls values in
        #   the above outline result in the same output, hence why firstPositive is 0 if positiveIndices
        #   is empty. The same index pattern is valid when the length of the zero array is 2.
        #   The direction for each zero should always be the next occurrence, except when up is True, then
        #   the first zero's previous occurrence should be found.
        if up:
            updatedZeros = zeroSorted[firstPositive-1:zCount] + zeroSorted[0:firstPositive-1]
            direction = [PREVIOUS_OCCURRENCE] + [NEXT_OCCURRENCE] * (zCount - 1)
        else:
            updatedZeros = zeroSorted[firstPositive:zCount] + zeroSorted[0:firstPositive]
            direction = [NEXT_OCCURRENCE] * zCount

        return updatedZeros, direction

    def _getTimeTo(self, zi: float, zj: float, jd: JulianDate, direction: int):
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = _angleDifference(lng - radians(self._geo.longitude))
        if direction == NEXT_OCCURRENCE and dl < 0:
            dl += TWOPI
        elif direction == PREVIOUS_OCCURRENCE and dl > 0:
            dl -= TWOPI

        # We need to use sidereal day length as one revolution, not solar day length.
        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    # Check the fraction but as if the satellite ellipse vectors didn't change.
    def _checkTimeStatic(self, satTime: JulianDate, checkTime: JulianDate):
        a, b, c = _getEllipseVectors(self._sat, satTime)

        zeta = self._geo.getZenithVector(checkTime)
        aDotZ = dot(a, zeta)
        bDotZ = dot(b, zeta)
        numerator = dot(zeta, self._geo.getPositionVector(checkTime) - c)
        denominator = sqrt(aDotZ*aDotZ + bDotZ*bDotZ)

        return numerator / denominator

    def _getMinimumValid(self, z: ZeroIntersection, k: list[float], jd, time, direction) -> (float, float):
        limit = self._rho * 0.999999999999999
        # The nudging factor should be proportional to rho.
        nudgeFactor = self._rho * 0.05
        zi = z.zi
        zj = _getZj(zi, self._rho, z.sign)
        derivativeSign = _signOf(_getZiPrime(zi, zj, self._zk, k))
        derivativeSign2 = derivativeSign
        previousTime = time
        previousZi = zi
        func = ZeroFunction.fromValues(self._zk, self._rho)

        while self._checkTime(time) > 1:
            if derivativeSign2 != derivativeSign:
                nudgeFactor /= 2
                time = previousTime
                zi = previousZi
            else:
                previousTime = time
                previousZi = zi

            tryZi = min(max(zi + derivativeSign * nudgeFactor, -limit), limit)
            zj = _getZj(zi, self._rho, z.sign)
            time = self._getTimeTo(zi, zj, jd, direction)
            ellipseVectors = _getEllipseVectors(self._sat, time)
            k = _getConstants(*ellipseVectors, self._geo, time)
            try:
                zi = func.computeSpecificZero(tryZi, k)
            except ZeroDisappeared:
                z.sign = -z.sign
                # Known issue: moving too far forwards/backwards moves the zero so much it's not
                #   recognized in the computeSpecificZero method. This could cause issues in the future.
                z = func.computeSpecificZero(z, k)

        while _getZiValue(zi, zj, self._zk, k) < 0:
            tmpZi = min(max(zi + derivativeSign * nudgeFactor, -limit), limit)
            zj = _getZj(tmpZi, self._rho, z.sign)
            if _getZiValue(tmpZi, zj, self._zk, k) < 0:
                if _signOf(_getZiPrime(tmpZi, zj, self._zk, k)) != derivativeSign:
                    # We've sandwiched the positive portion of the function, bisect from here.
                    pass

                temperZi = tmpZi
                tmpZj = _getZj(temperZi, self._rho, z.sign)
                while _getZiValue(temperZi, tmpZj, self._zk, k) < 0:
                    temperZi = (temperZi + tmpZi) / 2
                    tmpZj = _getZj(temperZi, self._rho, z.sign)

                zi = temperZi
            else:
                zi = tmpZi

                tmpZj = _getZj(tmpZi, self._rho, z.sign)
                temperZi = tmpZi
                while _signOf(_getZiPrime(temperZi, tmpZj, self._zk, k)) != derivativeSign:
                    temperZi = (temperZi + tmpZi) / 2
                    tmpZj = _getZj(temperZi, self._rho, z.sign)


            newNudge = nudgeFactor
            while _signOf(_getZiPrime(tmpZi, zj, self._zk, k)) != derivativeSign:
                tmpZi = (tmpZi + zi) / 2




        while (check := self._checkTimeStatic(time, newTime)) > 1:
            tmpZi = min(max(zi + derivativeSign * nudgeFactor, -limit), limit)
            zj = _getZj(tmpZi, self._rho, z.sign)
            newTime2 = self._getTimeTo(tmpZi, zj, jd, direction)
            if self._checkTimeStatic(time, newTime2) > 1:
                derivativeSign2 = _signOf(_getZiPrime(tmpZi, zj, self._zk, k))
                if derivativeSign2 != derivativeSign:
                    nudgeFactor /= 2
                else:
                    zi = tmpZi
                    newTime = newTime2
            else:
                zi = tmpZi
                newTime = newTime2

        # while (check := self._checkTime(time)) > 1:
        #     tmpZi = min(max(zi + derivativeSign * nudgeFactor, -limit), limit)
        #     zj = _getZj(tmpZi, self._rho, z.sign)
        #     time2 = self._getTimeTo(tmpZi, zj, jd, direction)
        #     if self._checkTime(time2) > 1:
        #         derivativeSign2 = _signOf(_getZiPrime(tmpZi, zj, self._zk, k))
        #         if derivativeSign2 != derivativeSign:
        #             nudgeFactor /= 2
        #         else:
        #             zi = tmpZi
        #             time = time2
        #             # derivativeSign = derivativeSign2
        #     else:
        #         zi = tmpZi
        #         time = time2

        return check, zi

    def _getMinimumValidOld(self, z: ZeroIntersection, k: list[float], jd, time, direction) -> (float, float):
        limit = self._rho * 0.999999999999999
        # The nudging factor should be proportional to rho.
        nudgeFactor = self._rho * 0.001
        zi = z.zi
        # tmpZi = zi

        # Want the value of zi to be > 0, so move in the direction of the derivative.
        zj = _getZj(zi, self._rho, z.sign)
        derivativeSign = _signOf(_getZiPrime(zi, zj, self._zk, k))

        # Nudge zi, then check zj at that time. If zi is still > 1, check derivative of that zi, if it's
        #   the opposite of where we started, half the nudgeFactor and try again. If the derivative has
        #   the same sign, nudge again and repeat.
        while (check := self._checkTime(time)) > 1:
            tmpZi = min(max(zi + derivativeSign * nudgeFactor, -limit), limit)
            # tmpZi = zi + derivativeSign * nudgeFactor
            zj = _getZj(tmpZi, self._rho, z.sign)
            time2 = self._getTimeTo(tmpZi, zj, jd, direction)
            if self._checkTime(time2) > 1:
                derivativeSign2 = _signOf(_getZiPrime(tmpZi, zj, self._zk, k))
                if derivativeSign2 != derivativeSign:
                    nudgeFactor /= 2
                else:
                    zi = tmpZi
                    time = time2
            else:
                zi = tmpZi
                time = time2

        return check, zi

    def _switchSides(self, z: ZeroIntersection, direction: int, time: JulianDate, k) -> float:
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        DT = 1. / 86400.
        zj = _getZj(z.zi, self._rho, z.sign)
        moveDirection = z.sign * _signOf(_getZiPrime(z.zi, zj, self._zk, k)) * -1

        while self._checkTime(time) > 1:
            # fixme:
            #   it's possible this moves too far past both zeros. If this happens we don't really need to
            #   care about this section, but functionally we need to worry about it to avoid infinite loops.
            time = time.future(moveDirection * DT)

        return time

    def _refineZero(self, zero: ZeroIntersection, jd: JulianDate, direction: int):
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        time = jd
        z = ZeroIntersection(zero.zi, zero.sign, zero.direction)
        func = ZeroFunction.fromValues(self._zk, self._rho)
        previousTime = time.future(-1)

        # Check must be close to 1 but also less than to avoid domain errors on inverse trig functions.
        #   To avoid being too strict on the check value and inducing an infinite loop, we also check the
        #   difference in successive times computed. The more liberal CHECK_MINIMUM value is backed up by
        #   ensuring the change in time is necessarily small.
        CHECK_EPSILON = 1 - CHECK_MINIMUM
        while abs((check := self._checkTime(time)) - 1) > CHECK_EPSILON or abs(time - previousTime) > TIME_DIFFERENCE:
        # while abs((check := self._checkTime(time))) > 1 or abs(time - previousTime) > TIME_DIFFERENCE:
            ellipseVectors = _getEllipseVectors(self._sat, time)
            k = _getConstants(*ellipseVectors, self._geo, time)
            # Update zi to the computed zero for the ZeroFunction with updated k values, which the
            #   previous zi as out initial guess for the zero finding procedure.
            # It is possible a zero may move from the positive to the negative function or vise versa.
            #   This is signaled by a ZeroDisappeared exception and should be handled here.
            try:
                z = func.computeSpecificZero(z, k)
            except ZeroDisappeared:
                z.sign = -z.sign
                # Known issue: moving too far forwards/backwards moves the zero so much it's not
                #   recognized in the computeSpecificZero method. This could cause issues in the future.
                z = func.computeSpecificZero(z, k)

            # Update the geo-position whose zenith vector corresponds to zi, and update the time to that point.
            zj = _getZj(z.zi, self._rho, z.sign)
            previousTime = time
            time = self._getTimeTo(z.zi, zj, jd, direction)

        if check > 1:
            return self._switchSides(z, direction, time, k)
        else:
            return time

    # def _refineZero_old(self, zero: ZeroIntersection, jd: JulianDate, direction: int):
    #     assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
    #                                                                              f'NEXT_OCCURRENCE or ' \
    #                                                                              f'PREVIOUS_OCCURRENCE, not{direction}'
    #     # These are named maximum, but they are the maximum which check is still less than 1.
    #     maximumValidCheck = -1e10
    #     maximumValidZi = zero.zi
    #     # firstIteration = True
    # 
    #     time = jd
    #     z = ZeroIntersection(zero.zi, zero.sign)
    #     zj = _getZj(z.zi, self._rho, zero.sign)
    #     ellipseVectors = _getEllipseVectors(self._sat, time)
    #     k = _getConstants(*ellipseVectors, self._geo, time)
    #     func = ZeroFunction.fromValues(self._zk, self._rho)
    #     previousTime = time.future(-1)
    # 
    #     # Check must be close to 1 but also less than to avoid domain errors on inverse trig functions.
    #     #   To avoid being too strict on the check value and inducing an infinite loop, we also check the
    #     #   difference in successive times computed. The more liberal CHECK_MINIMUM value is backed up by
    #     #   ensuring the change in time is necessarily small.
    #     while not ((CHECK_MINIMUM < (check := self._checkTime(time)) <= 1) and abs(time - previousTime) < TIME_DIFFERENCE):
    #         # If check is close enough but still greater than 1, we need to 'nudge' it in the proper direction.
    #         if 1 < check < 1.001:
    #             if maximumValidZi != zero.zi:
    #             # if not firstIteration:
    #                 z.zi = (z.zi + maximumValidZi) / 2
    #             else:
    #                 maximumValidCheck, maximumValidZi = self._getMinimumValid(z, k, jd, time, direction)
    #                 z.zi = (z.zi + maximumValidZi) / 2
    #         # Otherwise, just keep iterating towards the intersection geo-position.
    #         else:
    #             if maximumValidCheck < check < 1:
    #                 maximumValidCheck = check
    #                 maximumValidZi = z.zi
    #             ellipseVectors = _getEllipseVectors(self._sat, time)
    #             k = _getConstants(*ellipseVectors, self._geo, time)
    #             # Update zi to the computed zero for the ZeroFunction with updated k values, which the
    #             #   previous zi as out initial guess for the zero finding procedure.
    #             # It is possible a zero may move from the positive to the negative function or vise versa.
    #             #   This is signaled by a ZeroDisappeared exception and should be handled here.
    #             try:
    #                 z = func.computeSpecificZero(z, k)
    #             except ZeroDisappeared:
    #                 z.sign = -z.sign
    #                 # Known issue: moving too far forwards/backwards moves the zero so much it's not
    #                 #   recognized in the computeSpecificZero method. This could cause issues in the future.
    #                 z = func.computeSpecificZero(z, k)
    # 
    #         # Update the geo-position whose zenith vector corresponds to zi, and update the time to that point.
    #         zj = _getZj(z.zi, self._rho, z.sign)
    #         previousTime = time
    #         time = self._getTimeTo(z.zi, zj, jd, direction)
    # 
    #     return time

    def computeIntersectionTimes(self, jd: JulianDate):
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a JulianDate type, not {type(geo)}')

        # Get zeros and determine if any intersections occur.
        zeros = self._getGeneralZeros(jd)
        zCount = len(zeros)
        if zCount == 0:
            # The path is either always or never visible, so no intersection times exist.
            return []
        assert zCount == 2 or zCount == 4, f'expected to find 0, 2 or 4 intersections, not {zCount}'

        # Find the difference in longitude from each intersection geo-position relative to the current geo-position.
        # todo: need to add a property to GeoPosition class that returns the values in radians
        geoLng = radians(self._geo.longitude)
        lngs = [atan2(_getZj(z.zi, self._rho, z.sign), z.zi) - earthOffsetAngle(jd) for z in zeros]
        dls = [_angleDifference(lng - geoLng) for lng in lngs]

        # Sort the order and which occurrence of the zeros to find.
        zeroSorted, direction = self._orderZero(zeros, dls, zCount, jd)
        times = [self._refineZero(z, jd, d) for z, d in zip(zeroSorted, direction)]
        return times


class OrbitPath2:
    __slots__ = '_sat', '_geo', '_jd', '_zk', '_rho'

    def __init__(self, sat: Orbitable, geo: GeoPosition, jd: JulianDate):
        if not isinstance(sat, Orbitable):
            raise TypeError(f'sat must be an Orbitable type, not {type(sat)}')
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be a GeoPosition type, not {type(geo)}')
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a JulianDate type, not {type(geo)}')

        self._sat = sat
        self._geo = geo
        self._jd = jd

        angle = radians(geo.latitude)
        self._zk = sin(angle)
        self._rho = cos(angle)

    def _getGeneralZeros(self, jd: JulianDate) -> list[ZeroIntersection]:
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)

        func = ZeroFunction.fromValues(self._zk, self._rho)

        # return [ZeroIntersection(zi, 1) for zi in func.computeGeneralZeros(k)]
        return func.computeGeneralZeros(k)

    def _checkTime(self, jd: JulianDate) -> float:
        a, b, c = _getEllipseVectors(self._sat, jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        numerator = dot(zeta, gamma - c)
        aDotZ = dot(zeta, a)
        bDotZ = dot(zeta, b)
        denominator = sqrt(aDotZ * aDotZ + bDotZ * bDotZ)

        return numerator / denominator

    def _orderZero(self, zeros: list[ZeroIntersection], dls: list[float], zCount: int, jd: JulianDate) -> (
            list[ZeroIntersection], list[bool]):
        assert zCount == 2 or zCount == 4, f'expected zCount of 2 or 4, not {zCount}'

        # Sort zeros and angle differences based on the angle differences
        zeroSorted = [z for _, z in sorted(zip(dls, zeros))]
        dlSorted = sorted(dls)

        # We want to order the zeros to find the next time to the intersecting geo-positions, with
        #   the exception being we want to find the previous zero time if the orbit path is currently
        #   above the geo-position. We only need to sort the zeros accordingly, the sign of the dl is
        #   not significant outside this method, as the direction array specifies which direction to
        #   move in time while searching for intersecting times.
        # todo: should this be in [0, 1] or -[1, 1] ?
        up = -1 < self._checkTime(jd) < 1
        # Find the first positive dl then determine which indices to append to the positive values.
        positiveIndices = [i for i, dl, in enumerate(dlSorted) if dl > 0]
        if positiveIndices:
            firstPositive = positiveIndices[0]
        else:
            # The cases with no positive dls and all positive dls are treated the same, so set firstPositive
            #   to 0 if positiveIndices is empty (see explanation below).
            firstPositive = 0

        # Here is the breakdown for creating the updated zeroSorted arrays, note: fp = firstPositive,
        #   zs = zeroSorted and dls = dlSorted for brevity
        #
        #   dls: [1, 2, 3, 4] [-1, 2, 3, 4] [-1,-2, 3, 4] [-1,-2,-3, 4] [-1,-2,-3,-4]
        #   zs:  [a, b, c, d] [ a, b, c, d] [ a, b, c, d] [ a, b, c, d] [ a, b, c, d]
        #
        #   zs should become the following when up is true
        #   zs:  [d, a, b, c] [ a, b, c, d] [ b, c, d, a] [ c, d, a, b] [ d, a, b, c]
        #   and zs should become the following when up is false
        #   zs:  [a, b, c, d] [ b, c, d, a] [ c, d, a, b] [ d, a, b, c] [ a, b, c, d]
        #
        #   Therefore the correct ordering of zs should be zs[fp-1:length] + zs[0:fp-1] when up is True
        #   and zs[fp:length] + zs[0:fp] when up is False. Notice that the first and last dls values in
        #   the above outline result in the same output, hence why firstPositive is 0 if positiveIndices
        #   is empty. The same index pattern is valid when the length of the zero array is 2.
        #   The direction for each zero should always be the next occurrence, except when up is True, then
        #   the first zero's previous occurrence should be found.
        if up:
            updatedZeros = zeroSorted[firstPositive - 1:zCount] + zeroSorted[0:firstPositive - 1]
            direction = [PREVIOUS_OCCURRENCE] + [NEXT_OCCURRENCE] * (zCount - 1)
        else:
            updatedZeros = zeroSorted[firstPositive:zCount] + zeroSorted[0:firstPositive]
            direction = [NEXT_OCCURRENCE] * zCount

        return updatedZeros, direction

    def _getTimeTo(self, zi: float, zj: float, jd: JulianDate, direction: int):
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = _angleDifference(lng - radians(self._geo.longitude))
        if direction == NEXT_OCCURRENCE and dl < 0:
            dl += TWOPI
        elif direction == PREVIOUS_OCCURRENCE and dl > 0:
            dl -= TWOPI

        # We need to use sidereal day length as one revolution, not solar day length.
        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    def _switchSides(self, z: ZeroIntersection, direction: int, time: JulianDate, k) -> float:
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        DT = 1. / 86400.
        zj = _getZj(z.zi, self._rho, z.sign)
        moveDirection = z.sign * _signOf(_getZiPrime(z.zi, zj, self._zk, k)) * -1

        while self._checkTime(time) > 1:
            # fixme:
            #   it's possible this moves too far past both zeros. If this happens we don't really need to
            #   care about this section, but functionally we need to worry about it to avoid infinite loops.
            time = time.future(moveDirection * DT)

        return time

    def _refineZero(self, zero: ZeroIntersection, jd: JulianDate, direction: int):
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        time = jd
        z = copy(zero)
        func = ZeroFunction.fromValues(self._zk, self._rho)
        previousTime = time.future(-1)

        CHECK_EPSILON = 1 - CHECK_MINIMUM
        while abs((check := self._checkTime(time)) - 1) > CHECK_EPSILON or abs(time - previousTime) > TIME_DIFFERENCE:
            ellipseVectors = _getEllipseVectors(self._sat, time)
            k = _getConstants(*ellipseVectors, self._geo, time)

            # get the specified zero

    def _refineZero(self, zero: ZeroIntersection, jd: JulianDate, direction: int):
        assert direction == NEXT_OCCURRENCE or direction == PREVIOUS_OCCURRENCE, f'direction expected to be ' \
                                                                                 f'NEXT_OCCURRENCE or ' \
                                                                                 f'PREVIOUS_OCCURRENCE, not{direction}'
        time = jd
        z = ZeroIntersection(zero.zi, zero.sign, zero.direction)
        func = ZeroFunction.fromValues(self._zk, self._rho)
        previousTime = time.future(-1)

        # Check must be close to 1 but also less than to avoid domain errors on inverse trig functions.
        #   To avoid being too strict on the check value and inducing an infinite loop, we also check the
        #   difference in successive times computed. The more liberal CHECK_MINIMUM value is backed up by
        #   ensuring the change in time is necessarily small.
        CHECK_EPSILON = 1 - CHECK_MINIMUM
        while abs((check := self._checkTime(time)) - 1) > CHECK_EPSILON or abs(time - previousTime) > TIME_DIFFERENCE:
            # while abs((check := self._checkTime(time))) > 1 or abs(time - previousTime) > TIME_DIFFERENCE:
            ellipseVectors = _getEllipseVectors(self._sat, time)
            k = _getConstants(*ellipseVectors, self._geo, time)
            # Update zi to the computed zero for the ZeroFunction with updated k values, which the
            #   previous zi as out initial guess for the zero finding procedure.
            # It is possible a zero may move from the positive to the negative function or vise versa.
            #   This is signaled by a ZeroDisappeared exception and should be handled here.
            try:
                z = func.computeSpecificZero(z, k)
            except ZeroDisappeared:
                z.sign = -z.sign
                # Known issue: moving too far forwards/backwards moves the zero so much it's not
                #   recognized in the computeSpecificZero method. This could cause issues in the future.
                z = func.computeSpecificZero(z, k)

            # Update the geo-position whose zenith vector corresponds to zi, and update the time to that point.
            zj = _getZj(z.zi, self._rho, z.sign)
            previousTime = time
            time = self._getTimeTo(z.zi, zj, jd, direction)

        if check > 1:
            return self._switchSides(z, direction, time, k)
        else:
            return time

    def computeIntersectionTimes(self, jd: JulianDate):
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a JulianDate type, not {type(geo)}')

        # Get zeros and determine if any intersections occur.
        zeros = self._getGeneralZeros(jd)
        zCount = len(zeros)
        if zCount == 0:
            # The path is either always or never visible, so no intersection times exist.
            return []
        assert zCount == 2 or zCount == 4, f'expected to find 0, 2 or 4 intersections, not {zCount}'

        # Find the difference in longitude from each intersection geo-position relative to the current geo-position.
        # todo: need to add a property to GeoPosition class that returns the values in radians
        geoLng = radians(self._geo.longitude)
        lngs = [atan2(_getZj(z.zi, self._rho, z.sign), z.zi) - earthOffsetAngle(jd) for z in zeros]
        dls = [_angleDifference(lng - geoLng) for lng in lngs]

        # Sort the order and which occurrence of the zeros to find.
        zeroSorted, direction = self._orderZero(zeros, dls, zCount, jd)
        times = [self._refineZero(z, jd, d) for z, d in zip(zeroSorted, direction)]
        return times
