from copy import deepcopy
from enum import IntEnum, Enum
from math import sqrt, radians, cos, sin, pi, atan2, acos, asin, tan, atan

from _pyevspace import Vector, rotateMatrixFrom, getMatrixEuler, Angles, ZXZ, dot, vxcl, vang

from sattrack.exceptions import SatelliteAlwaysAbove, NoPassException
from sattrack.orbit.exceptions import RootCountChange
from sattrack.orbit.satellite import trueToMeanAnomaly, smaToMeanMotion
from sattrack.orbit.exceptions import PassedOrbitPathRange
from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.spacetime.juliandate import now
from sattrack.topocentric import getAltitude
from sattrack.util.constants import TWOPI, SIDEREAL_PER_SOLAR
from sattrack.util.conversions import atan3

# Imports solely for type checking
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.orbit.satellite import Orbitable
    from sattrack.spacetime.juliandate import JulianDate
    from sattrack.coordinates import GeoPosition


def computeZj(zi: float, rho: float, sign: int) -> float:
    return sign * sqrt(rho * rho - zi * zi)


def computeZjPrime(zi: float, zj: float) -> float:
    return -zi / zj


def computeZjPPrime(zi: float, zj: float, zjP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj

    return (zjP * zi - zj) / (zj * zj)


def computeZjPPPrime(zi: float, zj: float, zjP: float = None, zjPP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj
    if zjPP is None:
        zjPP = (zjP * zi - zj) / (zj * zj)

    zj2 = zj * zj
    return (zj2 * (zjP + zjPP * zi - zjP) - 2 * zj * zjP * (zjP * zi - zj)) / (zj2 * zj2)


def computeZiValue(zi: float, zj: float, zk: float, k: list[float]) -> float:
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


def computeZiPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj

    term1 = 2 * (k[0] * zi + k[1] * zj * zjP)
    term2 = k[3] * (zi * zjP + zj)
    term3 = k[4] * zk
    term4 = k[5] * zk * zjP

    return term1 + term2 + term3 + term4 + k[6] + k[7] * zjP


def computeZiPPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None, zjPP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj
    if zjPP is None:
        zjPP = (zjP * zi - zj) / (zj * zj)

    term1 = 2 * (k[0] + (k[1] * (zj * zjPP + zjP * zjP)))
    term2 = k[3] * (zi * zjPP + 2 * zjP)
    term3 = k[5] * zk * zjPP
    term4 = k[7] * zjPP

    return term1 + term2 + term3 + term4


def computeZiPPPrime(zi: float, zj: float, zk: float, k: list[float], zjP: float = None, zjPP: float = None,
                     zjPPP: float = None) -> float:
    if zjP is None:
        zjP = -zi / zj
    if zjPP is None:
        zjPP = (zjP * zi - zj) / (zj * zj)
    if zjPPP is None:
        zjPPP = computeZjPPPrime(zi, zj, zjP, zjPP)

    term1 = (2 * k[1]) * (zj * zjPPP + zjP * zjPP + 2 * zjP * zjPP)
    term2 = k[3] * (zi * zjPPP + 3 * zjPP)
    term3 = k[5] * zk * zjPPP
    term4 = k[7] * zjPPP

    return term1 + term2 + term3 + term4


def computeValue(zi: float, zj: float, zk: float, k: list[float], spec: 'FunctionSpec') -> 'Point':
    return FUNCTION_ARRAY[spec](zi, zj, zk, k)


# todo: move this to the Satellite object ?
def computeEllipseVectors(sat: 'Orbitable', jd: 'JulianDate') -> (Vector, Vector, Vector):
    elements = sat.getElements(jd)
    matrix = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    aMag = elements.sma
    bMag = aMag * sqrt(1 - elements.ecc * elements.ecc)
    cMag = aMag * elements.ecc

    a = rotateMatrixFrom(matrix, Vector(aMag, 0, 0))
    b = rotateMatrixFrom(matrix, Vector(0, bMag, 0))
    c = rotateMatrixFrom(matrix, Vector(-cMag, 0, 0))

    return a, b, c


def computeConstants(a: Vector, b: Vector, c: Vector, geo, jd) -> list[float]:
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


def hasSameSign(a: (float, int), b: (float, int)) -> bool:
    # Will return NotImplemented if either a or b is zero.

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


# todo: this might be useful elsewhere, maybe move to util module
def computeAngleDifference(angle: float) -> float:
    validAngle = angle % TWOPI
    if validAngle > pi:
        return validAngle - TWOPI
    return validAngle


# todo: this might be useful elsewhere, maybe move to util module
def signOf(value: float) -> int:
    # Zero will return True.

    if value >= 0:
        return 1
    return -1


class FunctionSpec(IntEnum):
    ZERO_FUNCTION = 0
    FIRST_DERIVATIVE = 1
    SECOND_DERIVATIVE = 2
    THIRD_DERIVATIVE = 3


class ExtremaSpec(IntEnum):
    EXTREMA = 0
    DOMAIN_BOUNDARY = 1


class OrbitPathDirection(IntEnum):
    UNKNOWN = 2
    ASCENDING = 1
    SPECIAL = 0
    DESCENDING = -1

    @classmethod
    def fromValue(cls, primeValue: float):
        """This can only guarantee an accurate direction if the primeValue does not correspond to an
        extrema value. That must be done by ensuring the point indicated by primeValue is also a boundary,
        in which case the OrbitPathDirection should be manually set to SPECIAL.

        Note, a primeValue equal to zero will be set to DESCENDING as a matter of implementation, however
        this type isn't necessarily useful when the derivative is zero, so it probably will not affect its use."""

        if primeValue > 0:
            return cls.ASCENDING
        else:
            return cls.DESCENDING


class OccurrenceDirection(Enum):
    PREVIOUS_OCCURRENCE = 0
    NEXT_OCCURRENCE = 1


# Used to select the function from a FunctionSpec enumeration value.
# This methodology requires the signature of these methods to be equivalent.
FUNCTION_ARRAY = [computeZiValue, computeZiPrime, computeZiPPrime, computeZiPPPrime]

# todo: should we put these in a config file?
# Constants and epsilon values.
DUPLICATE_ZERO_EPSILON = 1e-3
DOMAIN_ONE = 0.999999999999999
BIFURCATE_EPSILON = 1e-10
SPECIAL_EPSILON = 0.1
SWITCH_DT = (1 / 864000) / 2
CHECK_MINIMUM = 0.9999
CHECK_EPSILON = 1 - CHECK_MINIMUM
TIME_DIFFERENCE = 1e-5


class Point:

    __slots__ = '_zi', '_value', '_funcSpec', '_sign'

    def __init__(self, zi: float, value: float, sign: int, funcSpec: FunctionSpec = FunctionSpec.ZERO_FUNCTION):

        self._zi = zi
        self._value = value
        self._funcSpec = funcSpec
        self._sign = sign

    def __str__(self):
        return f'zi: {self._zi}, value: {self._value}, funcSpec: {self._funcSpec}, sign: {self._sign}'

    @property
    def zi(self):
        return self._zi

    @property
    def value(self):
        return self._value

    @property
    def functionSpec(self):
        return self._funcSpec

    @property
    def sign(self):
        return self._sign


class Extrema(Point):
    __slots__ = '_root', '_extremaSpec'

    def __init__(self, zi: float, value: float, sign: int, funcSpec: FunctionSpec = FunctionSpec.ZERO_FUNCTION,
                 extremaSpec: ExtremaSpec = ExtremaSpec.EXTREMA, root: 'FunctionRoot' = None):
        # If boundary is None, extremaSpec is expected to be DOMAIN_BOUNDARY
        super().__init__(zi, value, sign, funcSpec)
        self._root = root
        self._extremaSpec = extremaSpec

    @classmethod
    def fromZi(cls, zi: float, func: 'ZeroFunction', sign: int, root: 'FunctionRoot',
               funcSpec: FunctionSpec = FunctionSpec.ZERO_FUNCTION, extremaSpec: ExtremaSpec = ExtremaSpec.EXTREMA):

        zj = computeZj(zi, func.rho, sign)
        value = FUNCTION_ARRAY[funcSpec](zi, zj, func.zk, func.k)

        rtn = object.__new__(cls)
        rtn.__init__(zi, value, sign, funcSpec, extremaSpec, root)

        return rtn

    @classmethod
    def fromRoot(cls, root: 'FunctionRoot', func: 'ZeroFunction', sign: int):

        rtn = object.__new__(cls)

        funcSpec = FunctionSpec(root.point.functionSpec - 1)
        zj = computeZj(root.point.zi, func.rho, sign)
        value = computeValue(root.point.zi, zj, func.zk, func.k, funcSpec)

        if root.point == root.boundary.left or root.point == root.boundary.right:
            exSpec = ExtremaSpec.DOMAIN_BOUNDARY
        else:
            exSpec = ExtremaSpec.EXTREMA
        rtn.__init__(root.point.zi, value, sign, funcSpec, exSpec, root)

        return rtn

    def __str__(self):
        return super().__str__() + f'extremaSpec: {self._extremaSpec}'

    @property
    def root(self):
        return self._root

    @root.setter
    def root(self, value: 'FunctionRoot'):
        self._root = value

    @property
    def spec(self):
        return self._extremaSpec


class Boundary:
    __slots__ = '_upper', '_lower', '_left', '_right', '_middle'

    def __init__(self, point1: Point, point2: Point):
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

    def __str__(self):
        return str((self._upper, self._lower))

    def __repr__(self):
        return f'Boundary({self.__str__()})'

    def __contains__(self, item: float):
        return self._left.zi <= item <= self._right.zi

    def update(self, point: Point):

        if point.value > 0:
            self._upper = point
        else:
            self._lower = point

        if self._upper.zi < self._lower.zi:
            self._left = self._upper
            self._right = self._lower
        else:
            self._left = self._lower
            self._right = self._upper

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


class Domain(Boundary):

    @classmethod
    def fromZeroFunction(cls, func: 'ZeroFunction', sign: int, spec: FunctionSpec = FunctionSpec.ZERO_FUNCTION):

        rtn = object.__new__(cls)

        method = FUNCTION_ARRAY[spec]
        leftValue = method(-func.limit, computeZj(-func.limit, func.rho, sign), func.zk, func.k)
        rightValue = method(func.limit, computeZj(func.limit, func.rho, sign), func.zk, func.k)
        leftPoint = Extrema(-func.limit, leftValue, sign, spec, ExtremaSpec.DOMAIN_BOUNDARY, None)
        rightPoint = Extrema(func.limit, rightValue, sign, spec, ExtremaSpec.DOMAIN_BOUNDARY, None)

        rtn.__init__(leftPoint, rightPoint)

        return rtn


class FunctionRoot:
    __slots__ = '_point', '_sign', '_direction', '_boundary'

    def __init__(self, point: Point, direction: OrbitPathDirection, boundary: Boundary):
        self._point = point
        self._sign = point.sign
        self._direction = direction
        self._boundary = boundary

    def __str__(self):
        return str(self._point) + f'direction: {self._direction}'

    def __eq__(self, other: 'FunctionRoot'):
        if isinstance(other, float):
            return abs(self._point.zi - other) < DUPLICATE_ZERO_EPSILON
        elif isinstance(other, FunctionRoot):
            return abs(self._point.zi - other._point.zi) < DUPLICATE_ZERO_EPSILON and self._sign == other.sign
        return NotImplemented

    @property
    def point(self):
        return self._point

    @property
    def sign(self):
        return self._sign

    @property
    def direction(self):
        return self._direction

    # Need to be able to set direction to SPECIAL
    @direction.setter
    def direction(self, value: OrbitPathDirection):
        self._direction = value

    @property
    def boundary(self):
        return self._boundary


class RootList(list):

    def __getitem__(self, value):
        if isinstance(value, slice):
            tmp = super().__getitem__(value)
            return RootList(tmp)
        return super().__getitem__(value)

    def __add__(self, other):
        if isinstance(other, RootList):
            tmp = super().__add__(other)
            return RootList(tmp)
        return super().__add__(other)

    @staticmethod
    def _computeRelativeIndex(zi: float, array: list[FunctionRoot], comparator) -> int:
        relativeIndex = 0
        for root in array:
            if comparator(root.point.zi, zi):
                relativeIndex += 1
            else:
                break
        return relativeIndex

    def _orderRoots(self, zi: float, zj: float) -> 'RootList':
        # fixme: can this be done faster, without writing to new list and sorting then copying etc...
        zPos = sorted([z for z in self if z.sign == 1], key=lambda o: o.point.zi, reverse=True)
        zNeg = sorted([z for z in self if z.sign == -1], key=lambda o: o.point.zi)

        if zj > 0 or (zj == 0 and zi > 0):
            relativeIndex = self._computeRelativeIndex(zi, zPos, lambda a, b: a > b)
            orderedZeros = zPos[relativeIndex:] + zNeg + zPos[:relativeIndex]
        else:
            relativeIndex = self._computeRelativeIndex(zi, zNeg, lambda a, b: a < b)
            orderedZeros = zNeg[relativeIndex:] + zPos + zNeg[:relativeIndex]

        return RootList(orderedZeros)

    @staticmethod
    def _sortTwoRoots(orderedRoots: 'RootList', initRoots: 'RootList') -> 'RootList':
        firstRoot = initRoots[0]
        if firstRoot.direction != OrbitPathDirection.SPECIAL and \
                orderedRoots[0].direction != OrbitPathDirection.SPECIAL:
            if firstRoot.direction == orderedRoots[0].direction:
                return orderedRoots
            # return orderedRoots[1:] + orderedRoots[:1]
            return RootList(reversed(orderedRoots))
        # Otherwise, they both must be SPECIAL
        # Find the special boundary and return the order that represents initRoots.
        if initRoots[0].right == initRoots[1].left:
            if orderedRoots[0].right == orderedRoots[1].left:
                return orderedRoots
            # return orderedRoots[1:] + orderedRoots[:1]
            return RootList(reversed(orderedRoots))
        else:
            if orderedRoots[1].left == orderedRoots[0].right:
                return orderedRoots
            return RootList(reversed(orderedRoots))

    @staticmethod
    def _sortFourTwoSame(orderedRoots: 'RootList', firstRoot: FunctionRoot) -> 'RootList':
        # oldPosCount == 2 and newPosCount == 2
        idx = -1
        for i, root in enumerate(orderedRoots):
            if root.direction == firstRoot.direction and root.sign == firstRoot.sign:
                idx = i
                break
        return orderedRoots[idx:] + orderedRoots[:idx]

    @staticmethod
    def _sortFourTwoDifferent(orderedRoots: 'RootList', initRoots: 'RootList', twoToOne: dict) -> 'RootList':
        sameRoot = [root for root in twoToOne['old'] if root.direction == twoToOne['new'][0].direction][0]

        oldIndex = initRoots.index(sameRoot)
        # Only one value in twoToOne['new']
        newIndex = orderedRoots.index(twoToOne['new'][0])
        offset = newIndex - oldIndex

        return orderedRoots[offset:] + orderedRoots[:offset]

    @staticmethod
    def _sortFourThreeDifferent(orderedRoots: 'RootList', initRoots: 'RootList', twoToOne: dict) -> 'RootList':
        sameRoot = [root for root in twoToOne['new'] if root.direction == twoToOne['old'][0].direction][0]

        # Only one value in twoToOne['old']
        oldIndex = initRoots.index(twoToOne['old'][0])
        newIndex = orderedRoots.index(sameRoot)
        offset = newIndex - oldIndex

        return orderedRoots[offset:] + orderedRoots[:offset]

    @staticmethod
    def _sortFourThreeSame(orderedRoots: 'RootList', initRoots: 'RootList', singleRoots: dict) -> 'RootList':
        oldIndex = initRoots.index(singleRoots['old'])
        newIndex = orderedRoots.index(singleRoots['new'])
        offset = newIndex - oldIndex

        return orderedRoots[offset:] + orderedRoots[:offset]

    def _sortFourRoots(self, initRoots: 'RootList', orderedRoots: 'RootList', oldRoots: dict, newRoots: dict) \
            -> 'RootList':
        oldPosCount = len(oldRoots[1])
        newPosCount = len(newRoots[1])

        if oldPosCount == 2:
            if newPosCount == 2:
                return self._sortFourTwoSame(orderedRoots, initRoots[0])
            # New roots go from 2-2 split to 3-1
            else:
                if newPosCount == 1:
                    twoToOne = {'old': oldRoots[1], 'new': newRoots[1]}
                else:  # newPosCount == 3
                    twoToOne = {'old': oldRoots[-1], 'new': newRoots[-1]}

            return self._sortFourTwoDifferent(orderedRoots, initRoots, twoToOne)
        # Old roots are 3-1
        else:
            # New roots go from 3-1 to 2-2
            if newPosCount == 2:
                if oldPosCount == 1:
                    twoToOne = {'old': oldRoots[1], 'new': newRoots[1]}
                else:
                    twoToOne = {'old': oldRoots[-1], 'new': newRoots[-1]}

                return self._sortFourThreeDifferent(orderedRoots, initRoots, twoToOne)
            # We stay 3-1 to 3-1
            else:
                if newPosCount == 1:
                    singleRoots = {'old': oldRoots[1][0], 'new': newRoots[1][0]}
                else:  # newPosCount == 3
                    singleRoots = {'old': oldRoots[-1][0], 'new': newRoots[-1][0]}

                return self._sortFourThreeSame(orderedRoots, initRoots, singleRoots)

    def inOrder(self, zeta: Vector, initRoots: 'RootList' = None) -> 'RootList':
        """The number of roots in self and initRoots must be equal. This should ensure the derived list
        below won't be empty when calling the min function."""

        orderedRoots = self._orderRoots(zeta[0], zeta[1])

        if initRoots is None:
            return orderedRoots

        if len(orderedRoots) == 2:
            return self._sortTwoRoots(orderedRoots, initRoots)

        oldRoots = {
            1: [root for root in initRoots if root.sign == 1],
            -1: [root for root in initRoots if root.sign == -1]
        }
        newRoots = {
            1: [root for root in orderedRoots if root.sign == 1],
            -1: [root for root in orderedRoots if root.sign == -1]
        }

        return self._sortFourRoots(initRoots, orderedRoots, oldRoots, newRoots)


class ZeroFunction:
    __slots__ = '_geo', '_zk', '_rho', '_limit', '_k'

    def __init__(self, geo: 'GeoPosition'):
        self._geo = geo

        lat = radians(geo.latitude)
        self._rho = cos(lat)
        self._zk = sin(lat)
        self._limit = DOMAIN_ONE * self._rho

        self._k = None

    def setState(self, sat: 'Orbitable', jd: 'JulianDate'):
        # Will set the k array manually. Mostly useful for plotting specific instances in time.
        eVecs = computeEllipseVectors(sat, jd)
        self._k = computeConstants(*eVecs, self._geo, jd)

    def _getRootBifurcate(self, sign: int, spec: FunctionSpec) -> FunctionRoot:
        domain = Domain.fromZeroFunction(self, sign, spec)
        bound = deepcopy(domain)

        while abs(bound.left.zi - bound.right.zi) > BIFURCATE_EPSILON:
            middleZi = (bound.left.zi + bound.right.zi) / 2
            zj = computeZj(middleZi, self._rho, sign)
            middleValue = computeValue(middleZi, zj, self._zk, self._k, spec)
            middlePoint = Point(middleZi, middleValue, sign, spec)
            bound.update(middlePoint)

        # This OrbitPathDirection is unknown at the moment.
        return FunctionRoot(bound.left, OrbitPathDirection.UNKNOWN, domain)

    def _computeRoot(self, bound: Boundary, sign: int, spec: FunctionSpec) -> FunctionRoot:
        # fixme: do this with newton's method for better performance?
        boundCopy = deepcopy(bound)

        while abs(boundCopy.left.zi - boundCopy.right.zi) > BIFURCATE_EPSILON:
            middleZi = (boundCopy.left.zi + boundCopy.right.zi) / 2
            zj = computeZj(middleZi, self._rho, sign)
            middleValue = computeValue(middleZi, zj, self._zk, self._k, spec)
            middlePoint = Point(middleZi, middleValue, sign, spec)
            boundCopy.update(middlePoint)

        specUp = FunctionSpec(spec + 1)
        zi = boundCopy.left.zi
        zj = computeZj(zi, self._rho, sign)
        value = FUNCTION_ARRAY[specUp](zi, zj, self._zk, self._k)

        return FunctionRoot(boundCopy.left, OrbitPathDirection.fromValue(value * sign * -1), bound)

    def _getEvenExtrema(self, sign: int) -> list[Extrema]:
        # Get the root from the 2nd derivative and make an extrema value.
        rootDD = self._getRootBifurcate(sign, FunctionSpec.SECOND_DERIVATIVE)
        firstExtrema = Extrema.fromRoot(rootDD, self, sign)

        # Use the tails and the extrema to make left and right bounds of the first derivative.
        domain = Domain.fromZeroFunction(self, sign, FunctionSpec.FIRST_DERIVATIVE)
        leftBound = Boundary(domain.left, firstExtrema)
        rightBound = Boundary(firstExtrema, domain.right)

        # Use the bounds to find the two zeros, and the extrema they represent.
        leftRoot = self._computeRoot(leftBound, sign, FunctionSpec.FIRST_DERIVATIVE)
        rightRoot = self._computeRoot(rightBound, sign, FunctionSpec.FIRST_DERIVATIVE)
        leftExtrema = Extrema.fromRoot(leftRoot, self, sign)
        rightExtrema = Extrema.fromRoot(rightRoot, self, sign)

        return [leftExtrema, rightExtrema]

    def _getOddExtrema(self, sign: int) -> list[Extrema]:
        # Get root from the 3rd derivative as an extrema for the 2nd derivative.
        thirdRoot = self._getRootBifurcate(sign, FunctionSpec.THIRD_DERIVATIVE)
        secondExtrema = Extrema.fromRoot(thirdRoot, self, sign)

        # Compare the signs and the 2nd derivative tails and the extrema.
        domainDD = Domain.fromZeroFunction(self, sign, FunctionSpec.SECOND_DERIVATIVE)
        domainD = Domain.fromZeroFunction(self, sign, FunctionSpec.FIRST_DERIVATIVE)

        # If they are the same:
        if hasSameSign(domainDD.left.value, secondExtrema.value):
            boundary = Boundary(domainD.left, domainD.right)
            root = self._computeRoot(boundary, sign, FunctionSpec.FIRST_DERIVATIVE)
            return [Extrema.fromRoot(root, self, sign)]
        # If they are difference:
        else:
            # Find the roots to the 2nd derivative as extrema for the 1st derivative.
            leftDDBound = Boundary(domainDD.left, secondExtrema)
            rightDDBound = Boundary(secondExtrema, domainDD.right)
            leftDRoot = self._computeRoot(leftDDBound, sign, FunctionSpec.SECOND_DERIVATIVE)
            rightDRoot = self._computeRoot(rightDDBound, sign, FunctionSpec.SECOND_DERIVATIVE)
            leftDExtrema = Extrema.fromRoot(leftDRoot, self, sign)
            rightDExtrema = Extrema.fromRoot(rightDRoot, self, sign)

            # Use tails and extrema of the 1st derivative to check for roots.
            bounds = [Boundary(domainD.left, leftDExtrema), Boundary(leftDExtrema, rightDExtrema),
                      Boundary(rightDExtrema, domainD.right)]
            existingRoots = [bound for bound in bounds if not hasSameSign(bound.left.value, bound.right.value)]

            # Use the tails and valid extrema of the 1st derivative to get roots as extrema values for
            # the zero function.
            roots = [self._computeRoot(bound, sign, FunctionSpec.FIRST_DERIVATIVE) for bound in existingRoots]
            return [Extrema.fromRoot(root, self, sign) for root in roots]

    def _computeExtrema(self, sign: int) -> list[Extrema]:
        # Determine the number of extrema in the zero function based on the derivative of the tails.
        # If they have the same sign, there are an even number of extrema, odd if they differ in sign.
        domainD = Domain.fromZeroFunction(self, sign, FunctionSpec.FIRST_DERIVATIVE)

        if hasSameSign(domainD.left.value, domainD.right.value):
            return self._getEvenExtrema(sign)
        else:
            return self._getOddExtrema(sign)

    @staticmethod
    def _parseForSpecial(bounds: list[Boundary]) -> list[bool]:
        isSpecial = []
        for bound in bounds:
            extrema = min([bound.left, bound.right], key=lambda o: o.value)
            isSpecial.append(abs(extrema.value) < SPECIAL_EPSILON)

        return isSpecial

    def _findRoots(self, sign: int, ignoreSpecial=True) -> list[FunctionRoot]:
        extrema = self._computeExtrema(sign)
        domain = Domain.fromZeroFunction(self, sign, FunctionSpec.ZERO_FUNCTION)

        if not extrema:
            bounds = [domain]
        else:
            bounds = [Boundary(domain.left, extrema[0])] + \
                     [Boundary(extrema[i], extrema[i + 1]) for i in range(len(extrema) - 1)] + \
                     [Boundary(extrema[-1], domain.right)]

        existingRoots = [bound for bound in bounds if not hasSameSign(bound.left.value, bound.right.value)]
        isSpecial = self._parseForSpecial(existingRoots)

        if ignoreSpecial is True:
            rootBounds = [bound for bound, special in zip(existingRoots, isSpecial) if special is False]
            roots = [self._computeRoot(bound, sign, FunctionSpec.ZERO_FUNCTION) for bound in rootBounds]
        else:
            roots = []
            for bound, special in zip(existingRoots, isSpecial):
                root = self._computeRoot(bound, sign, FunctionSpec.ZERO_FUNCTION)
                if special is True:
                    root.direction = OrbitPathDirection.SPECIAL
                roots.append(root)

        return roots

    def computeRoots(self, sat: 'Orbitable', jd: 'JulianDate', ignoreSpecial: bool = True) -> RootList:
        self.setState(sat, jd)
        return RootList(self._findRoots(1, ignoreSpecial) + self._findRoots(-1, ignoreSpecial))

    @property
    def zk(self) -> float:
        return self._zk
    
    @property
    def rho(self) -> float:
        return self._rho
    
    @property
    def limit(self) -> float:
        return self._limit
    
    @property
    def k(self) -> list[float]:
        return self._k

    def plot(self, count: int = 10000, k: RootList = None):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('error: unable to import matplotlib')
            return

        if k is None:
            k = self._k

            m = int(self._limit * count)
            x = [-self._limit] + [i / count for i in range(-m, m)] + [self._limit]
            f, fPrime, fPPrime, fPPPrime = [], [], [], []
            fN, fPrimeN, fPPrimeN, fPPPrimeN = [], [], [], []

            for zi in x:
                zj = computeZj(zi, self._rho, 1)
                f.append(computeZiValue(zi, zj, self._zk, k))
                fPrime.append(computeZiPrime(zi, zj, self._zk, k))
                fPPrime.append(computeZiPPrime(zi, zj, self._zk, k))
                fPPPrime.append(computeZiPPPrime(zi, zj, self._zk, k))
                zj = computeZj(zi, self._rho, -1)
                fN.append(computeZiValue(zi, zj, self._zk, k))
                fPrimeN.append(computeZiPrime(zi, zj, self._zk, k))
                fPPrimeN.append(computeZiPPrime(zi, zj, self._zk, k))
                fPPPrimeN.append(computeZiPPPrime(zi, zj, self._zk, k))

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

            maxValue = max(max(f), max(fN)) * 1.2
            minValue = min(min(f), min(fN)) * 1.2
            ax.set_ylim((minValue, maxValue))
            ax.grid()

            plt.show()


class OrbitPath:
    __slots__ = '_sat', '_geo', '_jd', '_func', '_initCount', '_isAbove', '_zeta', '_initRoots'

    def __init__(self, sat: 'Orbitable', geo: 'GeoPosition'):
        self._sat = sat
        self._geo = geo
        self._func = ZeroFunction(geo)

        self._jd = None
        self._initCount = None
        self._isAbove = None
        self._zeta = None
        self._initRoots = None

    @property
    def orbitable(self) -> 'Orbitable':
        return self._sat

    @property
    def geo(self) -> 'GeoPosition':
        return self._geo

    def isPassPossible(self, time: 'JulianDate') -> bool:
        roots = self._getGeneralRoots(time, ignoreSpecial=True)
        # todo: we should account for geosynchronous orbis here
        return len(roots) != 0

    def getPathCount(self, time: 'JulianDate') -> int:
        roots = self._getGeneralRoots(time, ignoreSpecial=True)
        return divmod(roots, 2)[0]

    def isPathAbove(self, time: 'JulianDate') -> bool:
        return -1 <= self._checkTime(time) <= 1

    def _getGeneralRoots(self, jd: 'JulianDate', ignoreSpecial: bool = True) -> RootList:
        roots = self._func.computeRoots(self._sat, jd, ignoreSpecial)
        if roots:
            return roots.inOrder(self._zeta, self._initRoots)
        return roots

    def _checkTime(self, jd: 'JulianDate') -> float:
        a, b, c = computeEllipseVectors(self._sat, jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        numerator = dot(zeta, gamma - c)
        aDotZ = dot(zeta, a)
        bDotZ = dot(zeta, b)
        denominator = sqrt(aDotZ * aDotZ + bDotZ * bDotZ)

        return numerator / denominator

    def _computeTimeTo(self, zi: float, zj: float, jd: 'JulianDate', direction: OccurrenceDirection):
        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = computeAngleDifference(lng - radians(self._geo.longitude))
        if direction == OccurrenceDirection.NEXT_OCCURRENCE and dl < 0:
            dl += TWOPI
        elif direction == OccurrenceDirection.PREVIOUS_OCCURRENCE and dl > 0:
            dl -= TWOPI

        # Use sidereal day length as one revolution, not solar day length.
        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    def _computeSpecificRoot(self, jd: 'JulianDate', enum: int) -> FunctionRoot:
        roots = self._func.computeRoots(self._sat, jd)
        if len(roots) != self._initCount:
            raise RootCountChange(f'root count changed from {self._initCount} to {len(roots)}',
                                  (self._sat, self._geo, jd))

        orderedRoots = roots.inOrder(self._zeta, self._initRoots)

        return orderedRoots[enum - 1]

    def _switchSides(self, root: FunctionRoot, jd: 'JulianDate') -> float:
        zj = computeZj(root.point.zi, self._func.rho, root.sign)
        eVecs = computeEllipseVectors(self._sat, jd)
        k = computeConstants(*eVecs, self._geo, jd)
        prime = computeZiPrime(root.point.zi, zj, self._func.zk, k)
        moveDirection = root.sign * signOf(prime) * -1

        time = jd
        while self._checkTime(time) > 1:
            time = time.future(moveDirection * SWITCH_DT)

        return time

    def _refineRoot(self, root: FunctionRoot, enum: int, direction: OccurrenceDirection, start: 'JulianDate') \
            -> FunctionRoot:
        time = start
        root = deepcopy(root)
        prevTime = time.future(-1)

        while abs((check := self._checkTime(time)) - 1) > CHECK_EPSILON or abs(time - prevTime) > TIME_DIFFERENCE:
            root = self._computeSpecificRoot(time, enum)
            zj = computeZj(root.point.zi, self._func.rho, root.sign)
            prevTime = time
            time = self._computeTimeTo(root.point.zi, zj, start, direction)

        if check > 1:
            return self._switchSides(root, time)
        return time

    def _computeTimes(self, roots: RootList, directions: list[OccurrenceDirection]) -> list:
        times = []
        continueAgain = False
        skipLast = False

        for i, (root, direction) in enumerate(zip(roots, directions)):
            if skipLast is True and i == len(roots) - 1:
                continue
            if continueAgain:
                continueAgain = False
                continue

            if i == 0:
                jd = self._jd
            else:
                if direction == OccurrenceDirection.NEXT_OCCURRENCE:
                    jd = times[-1]
                else:
                    jd = times[0]

            try:
                time = self._refineRoot(root, i + 1, direction, jd)
                times.append(time)
            except RootCountChange:
                # Right now, if any zeros disappear, we abort looking for that pair.
                # If the current root is ascending, we skip it and its matching descending root:
                if root.direction == OrbitPathDirection.ASCENDING:
                    if i % 2 == 0:
                        continueAgain = True
                        continue
                    elif i == len(roots) - 1:
                        # This means i == 0 is the descending pair, which seems very unlikely to
                        # happen, but in case it does, this branch is here.
                        times.pop(0)
                    # This means i == 1 when len(roots) == 4
                    continueAgain = True
                    continue
                elif root.direction == OrbitPathDirection.DESCENDING:
                    if i % 2 == 1:
                        times.pop(i - 1)
                        continue
                    elif i == 2:
                        times.pop(1)
                        continue
                    # This means i == 0, and we need to skip the last iteration.
                    skipLast = True

        return times

    def _setState(self, jd: 'JulianDate'):
        self._jd = jd
        self._zeta = self._geo.getZenithVector(jd)

        if self._initRoots is not None:
            # In case this method is called twice, _computeGeneralZeros needs self._initRoots to be None
            # to properly call RootList.inOrder() with the correct parameters.
            self._initRoots = None
        roots = self._getGeneralRoots(jd)
        self._initRoots = roots

        self._initCount = len(roots)
        self._isAbove = self._checkTime(jd) < 1

    # todo: I think we can implement the number, but it needs to be thoroughly vetted. Saving the commented code for now
    # def computeOrbitPassTimes(self, jd, number: int = 1):
    def computeOrbitPassTimes(self, jd: 'JulianDate') -> tuple['JulianDate | None']:

        # self._jd = jd
        # self._zeta = self._geo.getZenithVector(jd)
        # self._func = ZeroFunction(self._geo)
        # if number > 0:
        #     self._setState(jd)
        # elif number < 0:
        #     self._setState(jd.future(-1))
        # else:    # number == 0:
        #     raise ValueError('number must be a non-zero integer, not {number}')

        self._setState(jd)
        if self._initCount == 0:
            return []

        # if self._initRoots is not None:
        #     # In case this method is called twice, _computeGeneralZeros needs self._initRoots to be None
        #     # to properly call RootList.inOrder() with the correct parameters.
        #     self._initRoots = None
        # roots = self._getGeneralRoots(jd)
        # self._initRoots = roots
        #
        # self._initCount = len(roots)
        # if self._initCount == 0:
        #     return []
        #
        # self._isAbove = self._checkTime(jd) < 1
        # todo: if we make the first one PREVIOUS, we might be able to find the previous path times
        directions = [OccurrenceDirection.NEXT_OCCURRENCE] * self._initCount
        if self._isAbove:
            directions[-1] = OccurrenceDirection.PREVIOUS_OCCURRENCE

        # runningTimes = []
        # # times = self._computeTimes(roots, directions)
        # times = self._computeTimes(self._initRoots, directions)
        # while True:
        #     # fixme: This might not be the most efficient way of doing this
        #     if len(times) == 0:
        #         # There are no more orbit passes, check if nth pass exists and return it or nothing.
        #         if number > 0:
        #             return runningTimes[number-1] if number <= len(runningTimes) else ()
        #         else:
        #             # runningTimes should be in reverse order based on the first rise time
        #             return runningTimes[number] if -number <= len(runningTimes) else ()
        #     if number > 0:
        #         if len(times) == 2:
        #             runningTimes.append((times[0], times[1]))
        #         else:
        #             runningTimes += [(times[0], times[1]), (times[2], times[3])]
        #         if len(runningTimes) >= number:
        #             return runningTimes[number - 1]
        #         nextTime = times[-1].future(0.0001)
        #     else:   # number < 0
        #         if len(times) == 2:
        #             runningTimes.insert(0, (times[0], times[1]))
        #         else:
        #             runningTimes = [(times[0], times[1]), (times[2], times[3])] + runningTimes
        #         if len(runningTimes) >= -number:
        #             return runningTimes[number]
        #         # fixme: This might not work great, can be greatly improved
        #         nextTime = times[0].future(-1)
        #
        #     self._setState(nextTime)
        #     directions = [OccurrenceDirection.NEXT_OCCURRENCE] * self._initCount
        #     if self._isAbove:
        #         directions[-1] = OccurrenceDirection.PREVIOUS_OCCURRENCE
        #     times = self._computeTimes(self._initRoots, directions)

        times = self._computeTimes(self._initRoots, directions)
        return tuple(sorted(times, key=lambda o: o.value))

    def _computeAnomalies(self, time: 'JulianDate') -> (float, float):
        a, b, c = computeEllipseVectors(self._sat, time)
        zeta = self._geo.getZenithVector(time)
        gamma = self._geo.getPositionVector(time)

        num = dot(zeta, gamma - c)
        zDotA = dot(zeta, a)
        zDotB = dot(zeta, b)
        den = sqrt(zDotA * zDotA + zDotB * zDotB)
        try:
            lhs = acos(num / den)
        except ValueError as e:
            if e.args[0] == 'math domain error':
                raise PassedOrbitPathRange()
            raise e
        rhs = atan3(zDotB, zDotA)

        t1 = (rhs - lhs) % TWOPI
        t2 = (rhs + lhs) % TWOPI

        return t1, t2

    def computePassTimesExec(self, startTime: 'JulianDate', rise: bool, nextOccurrence: bool) -> 'JulianDate':
        # Forward is true to find the next pass, false to find the previous.
        # Rise is true to find the rise time, false to find the set time.
        time = startTime
        prevTime = time.future(-1)
        firstPass = True

        while abs(time - prevTime) > TIME_DIFFERENCE * 10:

            elements = self._sat.getElements(time)
            t1, t2 = self._computeAnomalies(time)
            if rise is True:
                m = trueToMeanAnomaly(t1, elements.ecc)
            else:
                m = trueToMeanAnomaly(t2, elements.ecc)
            dm = computeAngleDifference(m - elements.meanAnomaly)

            # Here we'll adjust dm1 by n multiples of two pi to find the nth pass.
            if firstPass is True:
                if nextOccurrence is True and dm < 0:
                    dm += TWOPI
                elif nextOccurrence is False and dm > 0:
                    dm -= TWOPI
                firstPass = False

            n = smaToMeanMotion(elements.sma, self._sat.body.mu)
            dt = dm / n
            prevTime = time
            time = time.future(dt / 86400)

        return time

    def refineRiseSetTime(self, time: 'JulianDate') -> 'JulianDate':
        updatedTime = time
        alt = getAltitude(self._sat, self._geo, updatedTime)
        while abs(alt) > 4.848e-6:   # 1 arc-second
            topoPosition, topoVelocity = self._sat.getTopocentricState(self._geo, updatedTime)
            velExclude = vxcl(topoVelocity, topoPosition)

            a = asin(topoPosition[2] / topoPosition.mag())
            C = pi / 2
            B = vang(velExclude, -Vector.e3)

            # using cotangent four part formulae
            tmp = (cos(a) * cos(B) + (1 / tan(C)) * sin(B)) / sin(a)
            c = atan(1 / tmp)

            alpha = velExclude.mag() / topoPosition.mag()
            dt = c / alpha

            updatedTime = updatedTime.future(dt / 86400)
            alt = getAltitude(self._sat, self._geo, updatedTime)

        return updatedTime

    def _computeCurrentPass(self, time: 'JulianDate', topoState: Vector) -> ('JulianDate', 'JulianDate'):
        # This should only be called if topoState[0][2] > 0, i.e. the satellite is above the horizon at time.

        # This is an extremely rare case, but we should handle it nonetheless.
        if topoState[0][2] == 0:
            # The sat is rising at time.
            if topoState[1][2] > 0:
                setTime = self.computePassTimesExec(time, rise=False, nextOccurrence=True)
                setTime = self.refineRiseSetTime(setTime)
                riseTime = self.computePassTimesExec(setTime, rise=True, nextOccurrence=False)
                riseTime = self.refineRiseSetTime(riseTime)
            else:  # topoState[1][2] < 0   (not possible to be equal to zero because physics)
                riseTime = self.computePassTimesExec(time, rise=True, nextOccurrence=False)
                riseTime = self.refineRiseSetTime(riseTime)
                setTime = self.computePassTimesExec(riseTime, rise=False, nextOccurrence=True)
                setTime = self.refineRiseSetTime(setTime)
        else:
            # For fidelity, we should find the rise/set time that is 'farthest' from the current time, using
            # the anomalies to determine that.
            t1, t2 = self._computeAnomalies(time)
            elements = self._sat.getElements(time)
            if abs(t2 - elements.trueAnomaly) < abs(t1 - elements.trueAnomaly):
                # We're closer to set time, find rise with time.
                riseTime = self.computePassTimesExec(time, rise=True, nextOccurrence=False)
                riseTime = self.refineRiseSetTime(riseTime)
                setTime = self.computePassTimesExec(riseTime, rise=False, nextOccurrence=True)
                setTime = self.refineRiseSetTime(setTime)
            else:
                # We're closer to rise time, find set with time.
                setTime = self.computePassTimesExec(time, rise=False, nextOccurrence=True)
                setTime = self.refineRiseSetTime(setTime)
                riseTime = self.computePassTimesExec(setTime, rise=True, nextOccurrence=False)
                riseTime = self.refineRiseSetTime(riseTime)

        return riseTime, setTime

    def _computeOrbitPassRange(self, time: 'JulianDate', nextOccurrence: bool) -> ('JulianDate', 'JulianDate'):
        currentRange = self.computeOrbitPassTimes(time)
        nextRange = None
        if len(currentRange) == 0:
            currentRange = None
        elif len(currentRange) == 4:
            if nextOccurrence is True:
                nextRange = currentRange[2:]
                currentRange = currentRange[:2]
            else:
                nextRange = currentRange[:2]
                currentRange = currentRange[2:]
        return currentRange, nextRange

    def isGeosynchronous(self) -> bool:
        period = 1 / SIDEREAL_PER_SOLAR
        period = period * TWOPI / 86400
        delta = period * 0.01

        elements = self._sat.getElements(now())
        n = smaToMeanMotion(elements.sma, self._sat.body.mu)

        return period - delta <= n <= period + delta

    @staticmethod
    def _initRanges(time: 'JulianDate', nextOccurrence: bool, orbitPassTimes: tuple['JulianDate']) \
            -> ('JulianDate', 'tuple[JulianDate] | None', 'tuple[JulianDate] | None'):
        # orbitPassTimes shouldn't be empty here. That case should be handled in the calling method.
        timesLength = len(orbitPassTimes)
        assert timesLength != 0, f'no orbitPassTimes should be handled in calling method'

        if nextOccurrence is True:
            if timesLength == 2:
                # if time is before the valid range, set starting time to beginning of range
                if time < orbitPassTimes[0]:
                    time = orbitPassTimes[0]
                return time, orbitPassTimes, None
            else:   # timesLength == 4
                currentRange = orbitPassTimes[:2]
                nextRange = orbitPassTimes[2:]

                # if time is not during a valid range, set it to be
                if time < orbitPassTimes[0]:
                    time = orbitPassTimes[0]
                elif time < orbitPassTimes[2]:
                    time = orbitPassTimes[2]
                    currentRange = orbitPassTimes[2:]
                    nextRange = None
                return time, currentRange, nextRange
        else:   # nextOccurrence is False
            if timesLength == 2:
                # if time is after the valid range, set starting time to end of range
                if time > orbitPassTimes[1]:
                    time = orbitPassTimes[1]
                return time, orbitPassTimes, None
            else:   # timesLength == 4
                currentRange = orbitPassTimes[2:]
                nextRange = orbitPassTimes[:2]

                # if time is not during a valid range, set it to be
                if time > orbitPassTimes[3]:
                    time = orbitPassTimes[3]
                elif time > orbitPassTimes[1]:
                    time = orbitPassTimes[1]
                    currentRange = orbitPassTimes[:2]
                    nextRange = None
                return time, currentRange, nextRange

    def _checkPassValidity(self, time: 'JulianDate', orbitPassTimes: tuple['JulianDate'], topoState: [Vector, Vector])\
            -> bool:
        if not orbitPassTimes:
            if self.isPathAbove(time):
                # fixme: move this method to the Orbitable class
                if self.isGeosynchronous():
                    # If satellite is up, assume it's always up
                    if topoState[0][2] > 0:
                        raise SatelliteAlwaysAbove(f'{self._sat.name} is geosynchronous satellite that is '
                                                   f'always above {self._geo}')
                    # Else satellite is always down
                    raise NoPassException(f'{self._sat.name} is geosynchronous satellite that is never above '
                                          f'{self._geo}')
                else:
                    # It's safe to look for passes.
                    safe = True
            else:
                raise NoPassException(f'{self._sat.name} orbit path is not visible over {self._geo}')
        else:
            if self.isGeosynchronous():
                # fixme: this can't be correct 100% of the time
                if topoState[0][2] > 0:
                    # return orbitPassTimes
                    raise Exception('calling method should return orbitPassTimes')
                else:
                    raise NoPassException(f'{self._sat.name} is geosynchronous satellite that does not '
                                          f'pass over {self._geo}')
            else:
                # It's safe to look for passes, but with valid ranges.
                safe = False
        return safe

    def _computeNextPassTimes(self, time: 'JulianDate', currentRange: 'tuple[JulianDate]',
                              nextRange: 'tuple[JulianDate] | None') -> ('JulianDate', 'JulianDate'):
        while True:
            try:
                riseTime = self.computePassTimesExec(time, rise=True, nextOccurrence=True)
                break
            except PassedOrbitPathRange:
                if nextRange is not None:
                    currentRange = nextRange
                    nextRange = None
                else:
                    currentRange, nextRange = self._computeOrbitPassRange(currentRange[1].future(0.0001), True)
                    if not currentRange:
                        raise NoPassException(f'{self._sat.name} orbit path is no longer visible over {self._geo}')
                time = currentRange[0]

        riseTime = self.refineRiseSetTime(riseTime)
        setTime = self.computePassTimesExec(riseTime, rise=False, nextOccurrence=True)
        setTime = self.refineRiseSetTime(setTime)

        return riseTime, setTime

    def _computePrevPassTimes(self, time: 'JulianDate', currentRange: 'tuple[JulianDate]',
                              nextRange: 'tuple[JulianDate] | None') -> ('JulianDate', 'JulianDate'):
        while True:
            try:
                setTime = self.computePassTimesExec(time, rise=False, nextOccurrence=False)
                break
            except PassedOrbitPathRange:
                if nextRange is not None:
                    currentRange = nextRange
                    nextRange = None
                else:
                    currentRange, nextRange = self._computeOrbitPassRange(currentRange[1].future(-1), False)
                    if not currentRange:
                        raise NoPassException(f'{self._sat.name} orbit path is no longer visible over {self._geo}')
                    time = currentRange[1]

        setTime = self.refineRiseSetTime(setTime)
        riseTime = self.computePassTimesExec(setTime, rise=True, nextOccurrence=False)
        riseTime = self.refineRiseSetTime(riseTime)

        return riseTime, setTime

    def computeSatellitePassTimes(self, time: 'JulianDate | SatellitePass', nextOccurrence: bool = True)\
            -> ('JulianDate', 'JulianDate'):

        if nextOccurrence is not True and nextOccurrence is not False:
            raise ValueError(f'nextOccurrence must be True or False, not {nextOccurrence}')
        # We can't check if time is a SatellitePass and import his class (recursive import)
        if hasattr(time, 'maxInfo'):
            if nextOccurrence is True:
                time = time.setInfo.time.future(0.0001)
            else:
                time = time.riseInfo.time.future(-0.0001)

        topoState = self._sat.getTopocentricState(self._geo, time)
        orbitPassTimes = self.computeOrbitPassTimes(time)

        try:
            safe = self._checkPassValidity(time, orbitPassTimes, topoState)
        except Exception as e:
            if e.args[0] == 'calling method should return orbitPassTimes':
                return orbitPassTimes
            raise

        # If satellite is currently above, look for that pass.
        if topoState[0][2] >= 0:
            return self._computeCurrentPass(time, topoState)

        if safe is True:
            currentRange = ()
            nextRange = None
        else:
            time, currentRange, nextRange = self._initRanges(time, nextOccurrence, orbitPassTimes)

        # Compute actual pass here.
        if nextOccurrence is True:
            riseTime, setTime = self._computeNextPassTimes(time, currentRange, nextRange)
        else:
            riseTime, setTime = self._computePrevPassTimes(time, currentRange, nextRange)

        return riseTime, setTime
