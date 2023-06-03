import sys
from math import sqrt, cos, radians, pi, atan2

from _pyevspace import getMatrixEuler, ZXZ, Angles, rotateMatrixFrom, Vector, dot

from sattrack.coordinates import GeoPosition
from sattrack.orbit import Orbitable
from sattrack.spacetime.juliandate import JulianDate
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.util.constants import TWOPI, SIDEREAL_PER_SOLAR

#   DELETE THIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from sattrack.api import *
tle = TwoLineElement('''ISS (ZARYA)
1 25544U 98067A   23143.90152745  .00015170  00000+0  27201-3 0  9994
2 25544  51.6418  88.9189 0005432  13.3637  90.1823 15.50123922398022''')
iss = Satellite(tle)
jd = now()

from random import seed, uniform
seed()
testingGeo = 0
def testGeo():
    global testingGeo
    geo = GeoPosition(uniform(-90, 90), uniform(-180, 180))
    testingGeo = geo
    path = OrbitPath(iss, geo, jd)
    times = path.computeIntersectionTimes(jd)
    return geo, times


def plotZeroFunction(sat: Orbitable, geo: GeoPosition, jd: JulianDate):
    if 'matplotlib' in sys.modules:
        if not isinstance(sat, Orbitable):
            raise TypeError(f'sat must be a {Orbitable.__module__ + Orbitable.__name__} type, not {type(sat)}')
        if not isinstance(geo, GeoPosition):
            raise TypeError(f'geo must be a {GeoPosition.__module__ + GeoPosition.__name__} type, not {type(geo)}')
        if not isinstance(jd, JulianDate):
            raise TypeError(f'jd must be a {JulianDate.__module__ + JulianDate.__name__} type, not {type(jd)}')

        a, b, c = _getEllipseVectors(sat, jd)
        k = _getConstants(a, b, c, geo, jd)
        zk = geo.getZenithVector(jd)[2]
        rho = cos(radians(geo.latitude))

        m = int(rho * 1000)
        x = [i / 1000 for i in range(-m, m)]
        fPos, fNeg, fPrimePos, fPrimeNeg, fPPPos, fPPNeg = [], [], [], [], [], []

        for zi in x:
            zj = _getZj(zi, rho, 1)

            fPos.append(_getZiValue(zi, zj, zk, k))
            fNeg.append(_getZiValue(zi, -zj, zk, k))
            fPrimePos.append(_getZiPrime(zi, zj, zk, k))
            fPrimeNeg.append(_getZiPrime(zi, -zj, zk, k))
            fPPPos.append(_getZiPPrime(zi, zj, zk, k))
            fPPNeg.append(_getZiPPrime(zi, -zj, zk, k))

        figure = matplotlib.pyplot.figure()
        ax = figure.add_subplot(111)
        ax.plot(x, fPos, '-', c='blue')
        ax.plot(x, fPrimePos, '-', c='dodgerblue')
        ax.plot(x, fPPPos, '-', c='powderblue')
        ax.plot(x, fNeg, '-', c='red')
        ax.plot(x, fPrimeNeg, '-', c='coral')
        ax.plot(x, fPPNeg, '-', c='pink')

        maxVal = max(max(fPos), max(fNeg)) * 1.2
        minVal = min(min(fPos), min(fNeg)) * 1.2
        ax.set_ylim((minVal, maxVal))
        ax.grid()

        matplotlib.pyplot.show()
    else:
        raise ImportError('matplotlib not found in sys.modules')


def _getZj(zi: float, rho: float, plus: int) -> float:
    assert plus == 1 or plus == -1, f'plus expected to be +/- 1, was {plus}'
    # assert abs(zi) < rho, f'zi must be less than rho; zi: {zi}, rho: {rho}'

    return plus * sqrt(rho*rho - zi*zi)


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


EXTREMA_TAIL = 0
EXTREMA_MINMAX = 1


class Extrema:
    __slots__ = '_zi', '_value', '_ziPP', '_type'

    def __init__(self, zi: float, value: float, ziPP: float, xType: int):
        assert xType == EXTREMA_TAIL or xType == EXTREMA_MINMAX

        self._zi = zi
        self._value = value
        self._ziPP = ziPP
        self._type = xType

    def __str__(self) -> str:
        return str((self._zi, self._value, self._ziPP, self._type))

    def __repr__(self) -> str:
        return f'Extrema({str(self)})'

    @property
    def zi(self) -> float:
        return self._zi

    @property
    def value(self) -> float:
        return self._value

    @property
    def ziPP(self) -> float:
        return self._ziPP

    @property
    def type(self) -> int:
        return self._type


class Bound:
    __slots__ = '_lhs', '_rhs'

    def __init__(self, lhs: Extrema, rhs: Extrema):
        if not isinstance(lhs, Extrema):
            raise TypeError(f'lhs must be an Extrema type, not {type(lhs)}')
        if not isinstance(rhs, Extrema):
            raise TypeError(f'rhs must be an Extrema type, not {type(rhs)}')

        self._lhs = lhs
        self._rhs = rhs

    def __str__(self):
        return str((self._lhs, self._rhs))

    def __repr__(self):
        return f'Bound({self._lhs}, {self._rhs})'

    def __contains__(self, item):
        if not isinstance(item, float):
            raise TypeError(f'item must be a float type, not {type(item)}')

        return self._lhs.zi <= item <= self._rhs.zi

    @property
    def lhs(self):
        return self._lhs

    @property
    def rhs(self):
        return self._rhs


ALMOST_ONE = 0.9999999
EPSILON = 1e-7
NEWTON_EPSILON = 1e-7
EXTREMA_DERIVATIVE_EPSILON = 0.1


def _getZiNewton(zi: float, zk: float, rho: float, k: list[float], plus: int, epsilon: float = NEWTON_EPSILON):
    zi0 = -2
    zi1 = zi

    while abs(zi0 - zi1) > epsilon:
        zi0 = zi1
        zj = _getZj(zi0, rho, plus)
        zi1 = zi0 - _getZiValue(zi0, zj, zk, k) / _getZiPrime(zi0, zj, zk, k)

    return zi1


def _centerGuess(upper: float, lower: float, zi: float, zk: float, rho: float, k: list[float], plus: int) -> float:
    # lowerValue = _getZiValue(lower, _getZj(lower, rho, plus), zk, k)
    # upperValue = _getZiValue(upper, _getZj(upper, rho, plus), zk, k)
    ziValue = _getZiValue(zi, _getZj(zi, rho, plus), zk, k)

    if ziValue > 0:
        middle = (zi + lower) / 2
        return zi, middle, lower
    else:
        middle = (upper + zi) / 2
        return upper, middle, zi


def _getStartingZi(bound: Bound) -> float:
    if bound.lhs.type == EXTREMA_TAIL:
        return bound.lhs.zi
    elif bound.rhs.type == EXTREMA_TAIL:
        return bound.rhs.zi
    else:
        return (bound.lhs.zi + bound.rhs.zi) / 2


def _halfBounds(upper, lower, zi, zk, rho, k, plus):
    value = _getZiValue(zi, _getZj(zi, rho, plus), zk, k)
    if value >= 0:
        upper = zi
    else:
        lower = zi
    newZi = (upper + lower) / 2

    return upper, newZi, lower


def _getZiSafely(bound: Bound, zk: float, rho: float, k: list[float], plus: int, epsilon: float = NEWTON_EPSILON):
    # start by bisecting the bounds several times to severely reduce the changes of a domain error
    zi1 = _getStartingZi(bound)
    value = _getZiValue(bound.lhs.zi, _getZj(bound.lhs.zi, rho, plus), zk, k)
    if value > 0:
        upper = bound.lhs.zi
        lower = bound.rhs.zi
    else:
        upper = bound.rhs.zi
        lower = bound.lhs.zi
    INITIAL_REDUCTION_COUNT = 4
    for i in range(INITIAL_REDUCTION_COUNT):
        upper, zi1, lower = _halfBounds(upper, lower, zi1, zk, rho, k, plus)

    # To find the zero, we iterate over Newton's method, checking that the most recently calculated
    #   zi1 is in the bounds of the zero. This will keep the zi guess from being outside it's domain,
    #   and will also keep the guesses moving towards the actual zero even if it is still in a valid
    #   domain but not converging very quickly or is stagnant. This will protect against the domain
    #   problem and with the initial reduction above, should prevent slow convergence or even stagnation
    #   resulting in an infinite loop.
    zi0 = -2
    while abs(zi0 - zi1) > epsilon:
        # we don't know if upper is > or < lower, so need to check both cases, only 1 need be True
        if not (upper >= zi1 >= lower or lower >= zi1 >= upper):
            upper, zi1, lower = _halfBounds(upper, lower, zi0, zk, rho, k, plus)
        else:
            zj = _getZj(zi1, rho, plus)
            zi0 = zi1
            zi1 = zi0 - _getZiValue(zi0, zj, zk, k) / _getZiPrime(zi0, zj, zk, k)

    return zi1


def _getZiSafely_old(bound: Bound, zk: float, rho: float, k: list[float], plus: int, epsilon: float = NEWTON_EPSILON):
    zi0 = -2
    zi1 = _getStartingZi(bound)
    value = _getZiValue(bound.lhs.zi, _getZj(bound.lhs.zi, rho, plus), zk, k)
    if value > 0:
        upper = bound.lhs.zi
        lower = bound.rhs.zi
    else:
        upper = bound.rhs.zi
        lower = bound.lhs.zi
    centered = False

    while abs(zi0 - zi1) > epsilon or centered:
        centered = False
        try:
            zj = _getZj(zi1, rho, plus)
            zi0 = zi1
            zi1 = zi0 - _getZiValue(zi0, zj, zk, k) / _getZiPrime(zi0, zj, zk, k)
        except ValueError as e:
            if str(e) != 'math domain error':
                raise e
            else:
                # if an iteration of Newton's method pushed zi out of the valid domain,
                #   use a bifurcation method to improve zi guess until Newton's method converges

                upper, zi1, lower = _centerGuess(upper, lower, zi0, zk, rho, k, plus)
                # If the new zi1 results in the next zi1 being out of bounds, the old zi0 will produce
                #   the same zi1, and we'll end up in an infinite loop. Need to update zi0 as well, but
                #   this will make the epsilon check fail when it shouldn't, so need an extre boolean to
                #   keep that from happening.
                zi0 = zi1
                centered = True

    return zi1


def _checkFloat(value, name):
    if not isinstance(value, float):
        raise TypeError(f'{name} must be a float type, not {type(value)}')


def _checkArray(array, name):
    if not isinstance(array, (list, tuple)):
        raise TypeError(f'{name} must be a list or tuple, not {type(array)}')
    length = len(array)
    if length != 10:
        raise ValueError(f'{name} must have exactly 10 elements, {length} found')


class ZeroFunction:
    __slots__ = '_zk', '_rho', '_k', '_plus'

    def __init__(self, zk: float, rho: float, k: list[float], plus: int):
        assert plus == 1 or plus == -1, f'plus expected to be +/- 1, was {plus}'

        _checkFloat(zk, 'zk')
        _checkFloat(rho, 'rho')
        _checkArray(k, 'k')
        if not isinstance(plus, int):
            raise TypeError(f'plus must be an int type, not {type(plus)}')
        if plus != 1 and plus != -1:
            raise ValueError(f'plus must be equal to 1 or -1; plus: {plus}')

        self._zk = zk
        self._rho = rho
        self._k = k
        self._plus = plus

    NO_EXTREMA_FOUND = -1

    def _computeExtremaValues(self, start: int) -> tuple[float]:
        assert start == 1 or start == 0 or start == -1,\
            f'expected start value to be 1, 0 or -1; start: {start}'

        if start == 1:
            zi1 = self._rho * ALMOST_ONE
            zj = _getZj(zi1, self._rho, self._plus)
            ziP = _getZiPrime(zi1, zj, self._zk, self._k)
            ziPP = _getZiPPrime(zi1, zj, self._zk, self._k)

            # reduce initial guess until it's either positive and increasing or negative and decreasing
            while (ziP > 0 and ziPP < 0) or (ziP < 0 and ziPP > 0):
                zi1 *= 0.99
                ziP = _getZiPrime(zi1, zj, self._zk, self._k)
                ziPP = _getZiPPrime(zi1, zj, self._zk, self._k)
        elif start == -1:
            zi1 = -self._rho * ALMOST_ONE
            zj = _getZj(zi1, self._rho, self._plus)
            ziP = _getZiPrime(zi1, zj, self._zk, self._k)
            ziPP = _getZiPPrime(zi1, zj, self._zk, self._k)

            # reduce initial guess until it's either positive and decreasing or negative and increasing
            while (ziP > 0 and ziPP > 0) or (ziP < 0 and ziPP < 0):
                zi1 *= 0.99
                ziP = _getZiPrime(zi1, zj, self._zk, self._k)
                ziPP = _getZiPPrime(zi1, zj, self._zk, self._k)
        else:
            zi1 = 0

        zi0 = -2
        zj = _getZj(zi1, self._rho, self._plus)
        # FIXME: successive iterations might still be less than EPSILON, should check the value of zi1,
        #   or try bifurcation method
        # while abs(zi0 - zi1) > EPSILON:
        while abs(_getZiPrime(zi1, zj, self._zk, self._k)) > 0.1:
            zi0 = zi1
            try:
                # don't call _getZj because of the assert statement
                zj = self._plus * sqrt(self._rho*self._rho - zi0*zi0)
            except ValueError:
                # todo: reraise a custom exception here?
                return self.NO_EXTREMA_FOUND
            zi1 = zi0 - _getZiPrime(zi0, zj, self._zk, self._k) / _getZiPPrime(zi0, zj, self._zk, self._k)

        return zi1, _getZiValue(zi1, zj, self._zk, self._k), _getZiPPrime(zi1, zj, self._zk, self._k)

    def _getTails(self) -> (Extrema, Extrema):
        ziRight = self._rho * ALMOST_ONE
        zjRight = _getZj(ziRight, self._rho, self._plus)
        rightValue = _getZiValue(ziRight, zjRight, self._zk, self._k)
        rightPPrime = _getZiPPrime(ziRight, zjRight, self._zk, self._k)
        # the tail might also be an extrema, need to check
        if abs(_getZiPrime(ziRight, zjRight, self._zk, self._k)) < EXTREMA_DERIVATIVE_EPSILON:
            rightType = EXTREMA_MINMAX
        else:
            rightType = EXTREMA_TAIL
        rightTail = Extrema(ziRight, rightValue, rightPPrime, rightType)

        ziLeft = -self._rho * ALMOST_ONE
        zjLeft = _getZj(ziLeft, self._rho, self._plus)
        leftValue = _getZiValue(ziLeft, zjLeft, self._zk, self._k)
        leftPPrime = _getZiPPrime(ziLeft, zjLeft, self._zk, self._k)
        # the tail might also be an extrema, need to check
        if abs(_getZiPrime(ziLeft, zjLeft, self._zk, self._k)) < EXTREMA_DERIVATIVE_EPSILON:
            leftType = EXTREMA_MINMAX
        else:
            leftType = EXTREMA_TAIL
        leftTail = Extrema(ziLeft, leftValue, leftPPrime, leftType)

        return rightTail, leftTail

    def _checkZi(self, zi):
        _checkFloat(zi, 'zi')
        if abs(zi) > self._rho:
            raise ValueError(f'zi must be less than rho; zi: {zi}, rho: {self._rho}')

    @staticmethod
    def _checkZj(zj):
        _checkFloat(zj, 'zj')
        # todo: check zj value? unmake static if we do

    def getValue(self, zi: float, zj: float = None) -> float:
        self._checkZi(zi)

        if zj is not None:
            self._checkZj(zj)
        else:
            zj = _getZj(zi, self._rho, self._plus)

        return _getZiValue(zi, zj, self._zk, self._k)

    def getPrime(self, zi: float, zj: float = None) -> float:
        self._checkZi(zi)

        if zj is not None:
            self._checkZj(zj)
        else:
            zj = _getZj(zi, self._rho, self._plus)

        return _getZiPrime(zi, zj, self._zk, self._k)

    def getPPrime(self, zi, zj=None):
        self._checkZi(zi)

        if zj is not None:
            self._checkZj(zj)
        else:
            zj = _getZj(zi, self._rho, self._plus)

        return _getZiPPrime(zi, zj, self._zk, self._k)

    def _centerGuess(self, lhs: [float, int], rhs: [float, int]) -> float:
        lhsValue = _getZiValue(lhs[0], _getZj(lhs[0], self._rho, self._plus), self._zk, self._k)
        if lhsValue > 0:
            upper = lhs[0]
            lower = rhs[0]
        else:
            upper = rhs[0]
            lower = lhs[0]

        rhsValue = _getZiValue(rhs[0], _getZj(rhs[0], self._rho, self._plus), self._zk, self._k)
        # todo: find an analytic rational for this number (or change it)
        epsilon = min(abs(lhsValue), abs(rhsValue)) * 0.001
        middle = (lhs[0] + rhs[0]) / 2
        value = _getZiValue(middle, _getZj(middle, self._rho, self._plus), self._zk, self._k)
        while abs(value) > epsilon:
            if value > 0:
                upper = middle
            else:
                lower = middle
            middle = (upper + lower) / 2
            value = _getZiValue(middle, _getZj(middle, self._rho, self._plus), self._zk, self._k)

        return middle

    def computeGeneralZeros(self):
        rightTail, leftTail = self._getTails()
        xArgs = [self._computeExtremaValues(b) for b in (1, 0, -1)]
        extrema = [Extrema(*args, EXTREMA_MINMAX) for args in xArgs if args != self.NO_EXTREMA_FOUND]
        sortedExtrema = sorted([rightTail, leftTail] + extrema, key=lambda o: o.zi)

        bounds = []
        for i, x in enumerate(sortedExtrema[:-1]):
            if x.value * sortedExtrema[i+1].value < 0:
                bound = Bound(x, sortedExtrema[i+1])
                # bound = (x.zi, x.type), (sortedExtrema[i+1].zi, sortedExtrema[i+1].type)
                bounds.append(bound)

        # intersections = []
        # for lhs, rhs in bounds:
        #     # todo: do we bisect the tails too for precaution?
        #     if lhs[1] == EXTREMA_TAIL:
        #         zi = lhs[0]
        #     elif rhs[1] == EXTREMA_TAIL:
        #         zi = rhs[0]
        #     else:
        #         zi = self._centerGuess(lhs, rhs)
        #
        #     intersection = _getZiNewton(zi, self._zk, self._rho, self._k, self._plus)
        #     intersections.append(intersection)
        intersections = [_getZiSafely(bound, self._zk, self._rho, self._k, self._plus) for bound in bounds]

        return intersections

    def computeSpecificZero(self, zi: float, k: list[float]) -> float:
        self._checkZi(zi)

        _checkArray(k, 'k')

        # return _getZiNewton(zi, self._zk, self._rho, k, self._plus)
        rightTail, leftTail = self._getTails()
        xArgs = [self._computeExtremaValues(b) for b in (1, 0, -1)]
        extrema = [Extrema(*args, EXTREMA_MINMAX) for args in xArgs if args != self.NO_EXTREMA_FOUND]
        sortedExtrema = sorted([rightTail, leftTail] + extrema, key=lambda o: o.zi)

        bounds = []
        for i, x in enumerate(sortedExtrema[:-1]):
            if x.value * sortedExtrema[i + 1].value < 0:
                bound = Bound(x, sortedExtrema[i + 1])
                # bound = (x.zi, x.type), (sortedExtrema[i+1].zi, sortedExtrema[i+1].type)
                bounds.append(bound)

        containingBounds = [bound for bound in bounds if zi in bound]
        assert len(containingBounds) == 1, f'zi is contained in more than 1 bound'
        bound = containingBounds[0]

        return _getZiSafely(bound, self._zk, self._rho, k, self._plus)

    def plot(self):
        if __debug__:
            import matplotlib.pyplot as plt

            limit = self._rho * ALMOST_ONE
            plotCount = 1000
            m = int(limit * plotCount)
            x = [-limit] + [i / plotCount for i in range(-m, m)] + [limit]
            f, fPrime, fPPrime = [], [], []

            for zi in x:
                zj = _getZj(zi, self._rho, self._plus)
                f.append(_getZiValue(zi, zj, self._zk, self._k))
                fPrime.append(_getZiPrime(zi, zj, self._zk, self._k))
                fPPrime.append(_getZiPPrime(zi, zj, self._zk, self._k))

            figure = plt.figure()
            ax = figure.add_subplot(111)
            ax.plot(x, f, '-', c='blue')
            ax.plot(x, fPrime, '-', c='dodgerblue')
            ax.plot(x, fPPrime, '-', c='powderblue')

            maxVal = max(f) * 1.2
            minVal = min(f) * 1.2
            ax.set_ylim((minVal, maxVal))
            ax.grid()

            plt.show()
        else:
            raise NotImplemented('matplotlib is not a production dependency')


DUPLICATE_ZERO_EPSILON = 1e-3


class ZeroIntersection:

    def __init__(self, zi, rho, plus):
        assert plus == 1 or plus == -1, f'plus expected to be +/- 1, was {plus}'
        assert zi < rho, f'zi must be less than rho; zi: {zi}, rho: {rho}'

        self._zi = zi
        self._zj = _getZj(zi, rho, plus)
        self._plus = plus

    @property
    def zi(self):
        return self._zi

    @property
    def zj(self):
        return self._zj

    @property
    def plus(self):
        return self._plus

    def __str__(self):
        return str(self._zi)

    def __repr__(self):
        return str((self._zi, self._zj, self._plus))

    def __eq__(self, other):
        if isinstance(other, float):
            return abs(self._zi - float) < DUPLICATE_ZERO_EPSILON
        elif isinstance(other, ZeroIntersection):
            return abs(self._zi - other._zi) < DUPLICATE_ZERO_EPSILON
        return NotImplemented

    def __ne__(self, other):
        return self.__eq__(other)


def _angleDifference(angle):
    if angle > pi:
        return angle - TWOPI
    elif angle < -pi:
        return angle + TWOPI
    return angle


def _signOf(value):
    if value >= 0:
        return 1
    return -1


SMALL = 1e-7


class OrbitPath:
    __slots__ = '_sat', '_geo', '_jd', '_zk', '_rho'

    def __init__(self, sat: Orbitable, geo: GeoPosition, jd: JulianDate):
        self._sat = sat
        self._geo = geo
        self._jd = jd

        self._zk = geo.getZenithVector(jd)[2]
        self._rho = cos(radians(geo.latitude))

    # def _removeDuplicateZeros(self, zeros):
    #     length = len(zeros)
    #     if length != 0 and length != 2 and length != 4:
    #         for i, z0 in enumerate(zeros):
    #             for z1 in zeros[:i] + zeros[i+1:]:
    #                 if z0 == z1:
    #                     zeros.pop(i)
    #                     # zeros state has changed so start search over recursively
    #                     # note: saving the indices to pop and popping after will pop both duplicates and leave
    #                     #   an original missing from the list. recursion is chosen instead of a more complicated
    #                     #   popping algorithm
    #                     self._removeDuplicateZeros(zeros)

    def _getGeneralZeros(self, jd):
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)

        positiveFunction = ZeroFunction(self._zk, self._rho, k, 1)
        negativeFunction = ZeroFunction(self._zk, self._rho, k, -1)

        # zeros that occur in both positive and negative functions refer to different geo-positions
        #   (they have different zjs (one negative on positive) so they are not duplicates)
        zeros = [ZeroIntersection(zi, self._rho, 1) for zi in positiveFunction.computeGeneralZeros()]\
            + [ZeroIntersection(zi, self._rho, -1) for zi in negativeFunction.computeGeneralZeros()]
        # zeros = [(z, _getZj(z, self._rho, 1), 1) for z in positiveFunction.computeGeneralZeros()]\
        #         + [(z, _getZj(z, self._rho, -1), -1) for z in negativeFunction.computeGeneralZeros()]
        return zeros

    NEXT_OCCURRENCE = 1
    PREVIOUS_OCCURRENCE = 2

    def _getTimeTo(self, zi, zj, jd, direction):
        assert direction == self.NEXT_OCCURRENCE or direction == self.PREVIOUS_OCCURRENCE, \
            f'next expected to be 1 or 2, was {direction}'

        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = _angleDifference(lng - radians(self._geo.longitude))
        if direction == self.NEXT_OCCURRENCE and dl < 0:
            dl += TWOPI
        elif direction == self.PREVIOUS_OCCURRENCE and dl > 0:
            dl -= TWOPI

        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    def _checkTime(self, jd):
        a, b, c = _getEllipseVectors(self._sat, jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        num = dot(zeta, gamma - c)
        aDotZ = dot(zeta, a)
        bDotZ = dot(zeta, b)
        den = sqrt(aDotZ*aDotZ + bDotZ*bDotZ)

        return num / den

    def _refineZero(self, zi, jd, plus, direction):
        assert plus == 1 or plus == -1, f'plus expected to be +/- 1, was {plus}'
        assert zi < self._rho, f'zi must be less than rho; zi: {zi}, rho: {self._rho}'
        assert direction == self.NEXT_OCCURRENCE or direction == self.PREVIOUS_OCCURRENCE, \
            f'next expected to be 1 or 2, was {direction}'

        # todo: can we clean this up? we need k to instantiate func, but if we use a 'fake' array can we
        #   instantiate it without messing our computations up? we update the k value before using it anyway
        zj = _getZj(zi, self._rho, plus)
        ellipseVectors = _getEllipseVectors(self._sat, jd)
        k = _getConstants(*ellipseVectors, self._geo, jd)
        func = ZeroFunction(self._zk, self._rho, k, plus)
        time = jd

        # these are named maximum, but they are the maximum which check is still less than 1
        maximumValidCheck = -1e10
        maximumValidZi = zi
        firstIteration = True
        # check must be close to 1 but also less than to avoid domain errors on inverse trig functions
        while not (ALMOST_ONE < (check := self._checkTime(time)) <= 1):
            # if check is close enough but still greater than 1, we need to 'nudge' it in the proper direction
            if 1 < check < 1 + SMALL:
                # doing this on the first iteration will set zi equal to itself and create an infinite loop
                if not firstIteration:
                    zi = (zi + maximumValidZi) / 2
                else:
                    # want the value of zi to be > zero, so move in the direction of the derivative
                    sign = _signOf(func.getPrime(zi, zj))
                    zi += sign * SMALL
            # otherwise just keep iterating towards the intersecting geo-position
            else:
                if maximumValidCheck < check < 1:
                    maximumValidCheck = check
                    maximumValidZi = zi
                ellipseVectors = _getEllipseVectors(self._sat, time)
                k = _getConstants(*ellipseVectors, self._geo, time)
                # update zi to the computed zero for the ZeroFunction with updated k values, with the
                #   previous zi as our initial guess for the zero finding procedure
                zi = func.computeSpecificZero(zi, k)

            # update the geo-position whose zenith vector corresponds to zi, and update the time to that point
            zj = _getZj(zi, self._rho, plus)
            time = self._getTimeTo(zi, zj, jd, direction)
            if firstIteration:
                firstIteration = False

        return time

    def _orderZeros(self, zeros, dls, length, jd):
        assert length == 2 or length == 4, f'expected 2 or 4 zeros, got {length}'

        # sort zeros and dls based on dls
        zeroSorted = [z for _, z in sorted(zip(dls, zeros))]
        dlSorted = sorted(dls)

        # we want to order the zeros to find the next time to the intersecting geo-positions, with
        #   the exception being we want to find the previous zero time if the orbit path is currently
        #   above the geo-position. we only need to sort the zeros accordingly, the sign of the dl is
        #   not significant outside this method, as the direction array specifies which direction to
        #   move in time while searching for intersecting times
        # todo: investigate if this should be in [0, 1) or (-1, 1)
        up = -1 < self._checkTime(jd) < 1
        # find the first positive dl then determine which indices to append to the positive values
        positiveIndices = [i for i, dl in enumerate(dlSorted) if dl > 0]
        if positiveIndices:
            firstPositive = positiveIndices[0]
        else:
            # no positive dls and all positive dls are treated the same, so set firstPositive to 0 if
            #   positiveIndices is empty (see explanation below)
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
        #   therefore the correct ordering of zs should be zs[fp-1:length] + zs[0:fp-1] when up is True
        #   and zs[fp:length] + zs[0:fp] when up is False. Notice that the first and last dls values in
        #   the above outline result in the same output, hence why firstPositive is 0 if positiveIndices
        #   is empty. The same index pattern is valid when the length of the zero array is 2.
        #   The direction for each zero should always be the next occurrence, except when up is True, then
        #   the first zero's previous occurrence should be found.
        if up:
            updatedZeros = zeroSorted[firstPositive-1:length] + zeroSorted[0:firstPositive-1]
            direction = [self.PREVIOUS_OCCURRENCE] + [self.NEXT_OCCURRENCE] * (length - 1)
        else:
            updatedZeros = zeroSorted[firstPositive:length] + zeroSorted[0:firstPositive]
            direction = [self.NEXT_OCCURRENCE] * length

        return updatedZeros, direction

    def computeIntersectionTimes(self, jd):
        zeros = self._getGeneralZeros(jd)
        lenZeros = len(zeros)
        if lenZeros == 0:
            # the path is either always or never visible, so no intersection times exist
            return []
        elif lenZeros != 2 and lenZeros != 4:
            raise ValueError('expected to find 0, 2 or 4 intersections, found %i' % lenZeros)

        # find the delta in longitude from each intersecting geo-position, to current geo-position
        geoLng = radians(self._geo.longitude)
        lngs = [atan2(z.zj, z.zi) - earthOffsetAngle(jd) for z in zeros]
        dls = [_angleDifference(lng - geoLng) for lng in lngs]

        # sort the order and which occurrence of the zeros
        zeroSorted, direction = self._orderZeros(zeros, dls, lenZeros, jd)
        times = [self._refineZero(z.zi, jd, z.plus, d) for z, d in zip(zeroSorted, direction)]
        return times
