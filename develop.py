from math import radians, degrees, sin, cos, tan, asin, acos, atan, atan2

from matplotlib import pyplot as plt

from sattrack.api import *
from pyevspace import *

from sattrack.topocentric import _getPVector, _orbitAltitude, _timeToPlane

geo = GeoPosition(38, -98)
tle = TwoLineElement("""ISS (ZARYA)             
1 25544U 98067A   23130.56268595  .00014570  00000+0  26186-3 0  9997
2 25544  51.6402 155.0135 0006311 332.0920 164.8595 15.50089498395958""")
iss = Satellite(tle)
jd = now()



def getValues(sat, geo, jd):
    elements = sat.getElements(jd)
    mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))

    amag = elements.sma
    e = elements.ecc
    cmag = amag * e
    bmag = amag * sqrt(1-e*e)

    a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
    b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
    c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)
    phiPrime = geo.getGeocentricLatitude()
    rho = cos(phiPrime)

    d = gamma.mag() * zeta.mag() * cos(radians(geo.latitude) - phiPrime)
    k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
    k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
    k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
    k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
    k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
    k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
    k7 = 2 * d * c[0]
    k8 = 2 * d * c[1]
    k9 = 2 * d * c[2]
    k10 = d * d

    return {'elements': elements, 'a': a, 'b': b, 'c': c, 'zeta': zeta, 'gamma': gamma, 'rho': rho,
            'd': d, 'k': [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]}

def getZi(zi, zk, rho, k, plus=1):
    k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = k
    zj = plus * sqrt(rho*rho - zi*zi)
    # print('zj:', zj)

    term1 = k1 * zi * zi
    term2 = k2 * zj * zj
    term3 = k3 * zk * zk
    term4 = k4 * zi * zj
    term5 = k5 * zi * zk
    term6 = k6 * zj * zk
    term7 = k7 * zi
    term8 = k8 * zj
    term9 = k9 * zk

    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 - k10


def getZi2(zi, zk, rho, b, plus=1):
    zj = plus * sqrt(rho*rho - zi*zi)

    return zi * b[0] + zj * b[1] + zk * b[2]


def getConstants(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)
    phiPrime = geo.getGeocentricLatitude()
    # rho = cos(phiPrime)
    # zk = zeta[2]

    d = gamma.mag() * zeta.mag() * cos(radians(geo.latitude) - phiPrime)
    k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
    k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
    k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
    k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
    k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
    k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
    k7 = 2 * d * c[0]
    k8 = 2 * d * c[1]
    k9 = 2 * d * c[2]
    k10 = d * d

    return a, b, c, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]

def makeCosPlot(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)
    rho = cos(geo.getGeocentricLatitude())
    a, b, c, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10] = getConstants(sat, geo, jd)
    zk = zeta[2]

    limit = cos(geo.getGeocentricLatitude())
    m = int(limit * 1000)
    x = [i / 1000 for i in range(-m, m)]
    fPos = []
    fNeg = []
    fPrimePos = []
    fPrimeNeg = []
    fPPPos = []
    fPPNeg = []

    for zi in x:
        zjPos = sqrt(rho * rho - zi * zi)
        zjNeg = -zjPos
        term1 = (k1 * zi * zi)
        term2 = (k2 * zjPos * zjPos)
        term3 = (k3 * zk * zk)
        term4 = (k4 * zi * zjPos)
        term5 = (k5 * zi * zk)
        term6 = (k6 * zjPos * zk)
        term7 = (k7 * zi)
        term8 = (k8 * zjPos)
        term9 = (k9 * zk)
        FPos = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 - k10
        FNeg = (k1 * zi * zi) + (k2 * zjNeg * zjNeg) + (k3 * zk * zk) + (k4 * zi * zjNeg) + (k5 * zi * zk) + (k6 * zjNeg * zk) + (k7 * zi) + (k8 * zjNeg) + (k9 * zk) - k10

        try:
            zjPosP = -zi / zjPos
            FPrimePos = (2 * k1 * zi) + (2 * k2 * zjPos * zjPosP) + k4 * (zi * zjPosP + zjPos) + (k5 * zk) + (
                        k6 * zk * zjPosP) + (k7) + (k8 * zjPosP)
        except ZeroDivisionError:
            zjPos = 1e-7
            zjPosP = -zi / zjPos
            FPrimePos = 0
        try:
            zjNegP = -zi / zjNeg
            FPrimeNeg = (2 * k1 * zi) + (2 * k2 * zjNeg * zjNegP) + k4 * (zi * zjNegP + zjNeg) + (k5 * zk) + (
                        k6 * zk * zjNegP) + (k7) + (k8 * zjNegP)
        except ZeroDivisionError:
            zjNeg = 1e-7
            zjNegP = -zi / zjNeg
            FPrimeNeg = 0

        try:
            zjPosPP = (zjPosP * zi - zjPos) / (zjPos * zjPos)
            FPPPos = 2 * k1 + 2 * k2 * (zjPos * zjPosPP + zjPosP*zjPosP) + k4 * (zi * zjPosPP + 2 * zjPosP) + k6 * zk * zjPosPP + k8 * zjPosPP
        except ZeroDivisionError:
            FPPPos = 0
        try:
            zjNegPP = (zjNegP * zi - zjNeg) / (zjNeg * zjNeg)
            FPPNeg = 2 * k1 + 2 * k2 * (zjNeg * zjNegPP + zjNegP*zjNegP) + k4 * (zi * zjNegPP + 2 * zjNegP) + k6 * zk * zjNegPP + k8 * zjNegPP
        except ZeroDivisionError:
            FPPNeg = 0

        fPos.append(FPos)
        fNeg.append(FNeg)
        fPrimePos.append(FPrimePos)
        fPrimeNeg.append(FPrimeNeg)
        fPPPos.append(FPPPos)
        fPPNeg.append(FPPNeg)

    # truePosPrimeP = [0]
    # dt = 0.001
    # for i in range(len(fPos) - 1):
    #     # zi = x[i]
    #     df = fPrimeNeg[i+1] - fPrimeNeg[i]
    #     # df = getZi(zi+dt, zk, rho, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10], 1) - zi
    #     truePosPrimeP.append(df/dt)

    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.plot(x, fPos, '-', c='blue')
    ax.plot(x, fNeg, '-', c='red')
    ax.plot(x, fPrimePos, '-', c='dodgerblue')
    ax.plot(x, fPrimeNeg, '-', c='coral')

    maxVal = max(max(fPos), max(fNeg)) * 1.2
    minVal = min(min(fPos), min(fNeg)) * 1.2
    ax.set_ylim((minVal, maxVal))
    ax.grid()

    figure2 = plt.figure()
    ax2 = figure2.add_subplot(111)
    ax2.plot(x, fPrimePos, '-', c='dodgerblue')
    ax2.plot(x, fPrimeNeg, '-', c='coral')
    ax2.plot(x, fPPPos, '-', c='powderblue')
    ax2.plot(x, fPPNeg, '-', c='pink')
    # ax2.plot(x, truePosPrimeP, '-', c='lime')
    ax2.set_ylim((minVal, maxVal))
    ax2.grid()

    plt.show()


def getZiP(zi, zj, zk, k):
    zjP = -zi / zj
    return (2 * k[0] * zi) + (2 * k[1] * zj * zjP) + (k[3] * (zi * zjP + zj)) + (k[4] * zk) + (k[5] * zk * zjP) + (k[6]) + k[7] * zjP


def getZiPP(zi, zj, zk, k):
    zjP = -zi / zj
    zjPP = (zjP * zi - zj) / (zj*zj)
    return (2 * k[0]) + (2 * k[1] * (zj * zjPP + zjP * zjP)) + (k[3] * (zi * zjPP + 2 * zjP)) + (k[5] * zk * zjPP) + (k[7] * zjPP)


tle = TwoLineElement('''NOAA 1 [-]
1 04793U 70106A   23134.16528281 -.00000029  00000+0  91342-4 0  9997
2 04793 101.5043 194.0017 0031725  39.1298  23.1377 12.54015040399041''')
sso = Satellite(tle)

tle = TwoLineElement('''MOLNIYA 2-9             
1 07276U 74026A   23134.29008326  .00000179  00000+0  00000+0 0  9999
2 07276  64.2798 182.1851 6382064 276.7617  19.3773  2.45098251257494''')
molniya = Satellite(tle)


def computeExtremaValues(sat, geo, jd, plus=1, start=True):
    assert plus == 1 or plus == -1
    assert start is True or start is False or start == 0

    a, b, c, k = getConstants(sat, geo, jd)
    rho = cos(geo.getGeocentricLatitude())
    zeta = geo.getZenithVector(jd)
    zk = zeta[2]

    # frac = 0.99
    # dt = 0.01
    if start is True:
        zi1 = rho * 0.99999
        zj = plus * sqrt(rho*rho - zi1*zi1)
        ziP = getZiP(zi1, zj, zk, k)
        ziPP = getZiPP(zi1, zj, zk, k)
        # reduce initial guess until it is either positive and increasing or negative and decreasing
        while ziP > 0 and ziPP < 0 or ziP < 0 and ziPP > 0:
            zi1 *= 0.99
            ziP = getZiP(zi1, zj, zk, k)
            ziPP = getZiPP(zi1, zj, zk, k)
    elif start is False:
        zi1 = -rho * 0.99999
        zj = plus * sqrt(rho*rho - zi1*zi1)
        ziP = getZiP(zi1, zj, zk, k)
        ziPP = getZiPP(zi1, zj, zk, k)
        # reduce initial guess until it is either positive and decreasing or negative and increasing
        while ziP > 0 and ziPP > 0 or ziP < 0 and ziPP < 0:
            zi1 *= 0.99
            ziP = getZiP(zi1, zj, zk, k)
            ziPP = getZiPP(zi1, zj, zk, k)
    else:
        zi1 = 0
    zi0 = -2

    while abs(zi0 - zi1) > 1e-7:
        zi0 = zi1
        try:
            zj = plus * sqrt(rho*rho - zi1*zi1)
        except ValueError:
            return None, None, None
        zi1 = zi0 - getZiP(zi0, zj, zk, k) / getZiPP(zi0, zj, zk, k)

    return zi1, getZi(zi1, zk, rho, k, plus), getZiPP(zi1, zj, zk, k)


class Extrema:
    __slots__ = '_zi', '_value', '_ziPP'

    def __init__(self, zi, value, ziPP):
        self._zi = zi
        self._value = value
        self._ziPP = ziPP

    def __str__(self):
        return str((self._zi, self._value, self._ziPP))

    def __repr__(self):
        return str((self._zi, self._value, self._ziPP))

    @property
    def minimum(self):
        return self._ziPP > 0

    @property
    def maximum(self):
        return self._ziPP < 0
    
    @property
    def zi(self):
        return self._zi
    
    @property
    def value(self):
        return self._value

    def __eq__(self, other):
        return abs(self._zi - other.zi) <= 1e-5

    def __ne__(self, other):
        return abs(self._zi - other.zi) > 1e-5

    def __hash__(self):
        return hash((self._zi, self._value, self._ziPP))


def getExtrema(sat, geo, jd):
    posRight = Extrema(*computeExtremaValues(sat, geo, jd, 1, True))
    posLeft = Extrema(*computeExtremaValues(sat, geo, jd, 1, False))
    posZero = Extrema(*computeExtremaValues(sat, geo, jd, 1, 0))
    negRight = Extrema(*computeExtremaValues(sat, geo, jd, -1, True))
    negLeft = Extrema(*computeExtremaValues(sat, geo, jd, -1, False))
    negZero = Extrema(*computeExtremaValues(sat, geo, jd, -1, 0))

    posList = []
    for i in (posRight, posZero, posLeft):
        if i.zi is not None:
            if posList:
                if all(abs(i.zi - j.zi) > 1e-5 for j in posList):
                    posList.append(i)
            else:
                posList.append(i)
    negList = []
    for i in (negRight, negZero, negLeft):
        if i.zi is not None:
            if negList:
                if all(abs(i.zi - j.zi) > 1e-5 for j in negList):
                    negList.append(i)
            else:
                negList.append(i)

    return {'pos': posList, 'neg': negList}


def getGlobalExtrema(sat, geo, jd):
    extrema = getExtrema(sat, geo, jd)
    pos = extrema['pos']
    neg = extrema['neg']

    posMaxList = [i for i in pos if i.maximum]
    posMinList = [i for i in pos if i.minimum]
    negMaxList = [i for i in neg if i.maximum]
    negMinList = [i for i in neg if i.minimum]

    if posMaxList:
        posMax = max(posMaxList, key=lambda o: o.value)
    else:
        posMax = None
    if posMinList:
        posMin = min(posMinList, key=lambda o: o.value)
    else:
        posMin = None
    if negMaxList:
        negMax = max(negMaxList, key=lambda o: o.value)
    else:
        negMax = None
    if negMinList:
        negMin = min(negMinList, key=lambda o: o.value)
    else:
        negMin = None

    return {'pos': (posMin, posMax), 'neg': (negMin, negMax)}


# make plot object, properties: intersections(zeros),
    # plot has function objects

# function object, properties: odd, even, zeros, extrema, minimum, maximum


# class will take satellite vectors, geo values and plus/minus value and will spit out general zi zeros
# then (the main purpose) is to update the class with time specific satellite vectors and refine towards a specific
#   zi zero. Two or four instances of this class will be used in the OrbitalPath class to find the times.
class ZeroFunction:

    def __init__(self, zk, rho, k, plus):
        assert plus == 1 or plus == -1

        self._zk = zk
        self._rho = rho
        self._k = k
        self._plus = plus

    def _getZj(self, zi):
        return self._plus * sqrt(self._rho*self._rho - zi*zi)

    def _computeExtremaValues(self, start):
        assert start is True or start is False or start == 0

        if start is True:
            zi1 = self._rho * 0.99999
            zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            ziP = getZiP(zi1, zj, self._zk, self._k)
            ziPP = getZiPP(zi1, zj, self._zk, self._k)
            # reduce initial guess until it is either positive and increasing or negative and decreasing
            while ziP > 0 and ziPP < 0 or ziP < 0 and ziPP > 0:
                zi1 *= 0.99
                ziP = getZiP(zi1, zj, self._zk, self._k)
                ziPP = getZiPP(zi1, zj, self._zk, self._k)
        elif start is False:
            zi1 = -self._rho * 0.99999
            zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            ziP = getZiP(zi1, zj, self._zk, self._k)
            ziPP = getZiPP(zi1, zj, self._zk, self._k)
            # reduce initial guess until it is either positive and decreasing or negative and increasing
            while ziP > 0 and ziPP > 0 or ziP < 0 and ziPP < 0:
                zi1 *= 0.99
                ziP = getZiP(zi1, zj, self._zk, self._k)
                ziPP = getZiPP(zi1, zj, self._zk, self._k)
        else:
            zi1 = 0
        zi0 = -2

        while abs(zi0 - zi1) > 1e-7:
            zi0 = zi1
            try:
                zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            except ValueError:
                return None, None, None
            zi1 = zi0 - getZiP(zi0, zj, self._zk, self._k) / getZiPP(zi0, zj, self._zk, self._k)

        return zi1, getZi(zi1, self._zk, self._rho, self._k, self._plus), getZiPP(zi1, zj, self._zk, self._k)

    def _getTails(self):
        ziRight = self._rho * 0.99999
        ziLeft = -self._rho * 0.99999
        rightTail = Extrema(ziRight, self.getValue(ziRight), self.getPPrime(ziRight))
        leftTail = Extrema(ziLeft, self.getValue(ziLeft), self.getPPrime(ziLeft))

        return rightTail, leftTail

    def _getValue(self, zi, zj, k):
        # if zj is None:
        #     zj = self._getZj(zi)

        term1 = k[0] * zi * zi
        term2 = k[1] * zj * zj
        term3 = k[2] * self._zk * self._zk
        term4 = k[3] * zi * zj
        term5 = k[4] * zi * self._zk
        term6 = k[5] * zj * self._zk
        term7 = k[6] * zi
        term8 = k[7] * zj
        term9 = k[8] * self._zk

        return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 - k[9]

    def _getPrime(self, zi, zj, k):
        zjP = -zi / zj
        return (2 * k[0] * zi) + (2 * k[1] * zj * zjP) + (k[3] * (zi * zjP + zj)) + (k[4]
                                        * self._zk) + (k[5] * self._zk * zjP) + (k[6]) + k[7] * zjP

    def _getPPrime(self, zi, zj, k):
        zjP = -zi / zj
        zjPP = (zjP * zi - zj) / (zj * zj)
        return (2 * k[0]) + (2 * k[1] * (zj * zjPP + zjP * zjP)) + (k[3] * (zi * zjPP + 2 * zjP)) + (
                k[5] * self._zk * zjPP) + (k[7] * zjPP)

    def getValue(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)

        return self._getValue(zi, zj, self._k)

    def getPrime(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)

        return self._getPrime(zi, zj, self._k)

    def getPPrime(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)

        return self._getPPrime(zi, zj, self._k)

    def plot(self):
        m = int(self._rho * .99999 * 1000)
        x = [i / 1000 for i in range(-m, m)]
        f = []
        fPrime = []
        fPPrime = []

        for zi in x:
            zj = self._getZj(zi)

            f.append(self.getValue(zi, zj))
            fPrime.append(self.getPrime(zi, zj))
            fPPrime.append(self.getPPrime(zi, zj))

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

    def computeGeneralZeros(self):
        rightTail, leftTail = self._getTails()

        xArgs = [self._computeExtremaValues(s) for s in (True, 0, False)]
        xt = [Extrema(*args) for args in xArgs if args != (None, None, None)]
        extrema = set(xt)

        sortedExtrema = sorted([rightTail, leftTail] + list(extrema), key=lambda o: o.zi)
        zeros = []
        for i, x in enumerate(sortedExtrema[:-1]):
            if x.value * sortedExtrema[i+1].value < 0:
                z = x.zi, sortedExtrema[i+1].zi
                zeros.append(z)

        intersections = []
        for z in zeros:
            zi0 = -2
            if self.getPrime(z[0]) > 1e3:
                zi1 = z[0]
            elif self.getPrime(z[1]) > 1e3:
                zi1 = z[1]
            else:
                zi1 = (z[0] + z[1]) / 2
            while abs(zi0 - zi1) > 1e-7:
                zi0 = zi1
                zi1 = zi0 - self.getValue(zi0) / self.getPrime(zi0)
            intersections.append(zi1)

        return intersections

    def computeSpecificZero(self, zi, k):
        zi0 = -2
        zi1 = zi
        while abs(zi0 - zi1) > 1e-7:
            zi0 = zi1
            zj = self._plus * sqrt(self._rho*self._rho - zi0*zi0)
            zi1 = zi0 - self._getValue(zi0, zj, k) / self._getPrime(zi0, zj, k)

        return zi1


def getTimeTo(zi, geo, jd, plus=1):
    rho = cos(radians(geo.latitude))
    zj = plus * sqrt(rho*rho - zi*zi)

    lng = atan2(zj, zi) - earthOffsetAngle(jd)
    dl = lng - radians(geo.longitude)
    dt = SIDEREAL_PER_SOLAR * dl / TWOPI

    return jd.future(dt)


def checkZi(zi, sat, geo, jd, plus=1):
    rho = cos(radians(geo.latitude))
    zj = plus * sqrt(rho * rho - zi * zi)

    a, b, c, k = getConstants(sat, geo, jd)
    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)

    num = dot(zeta, gamma - c)
    den = sqrt(dot(zeta, a)**2 + dot(zeta, b)**2)

    return num, den, num / den


class Function:
    __slots__ = '_zk', '_rho', '_a', '_b', '_c', '_k', '_plus', '_rightTail', '_leftTail', '_extrema', \
                '_zeros', '_zeroCount', '_intersections'

    def __init__(self, zk, rho, a, b, c, k, plus):
        self._zk = zk
        self._rho = rho
        self._a = a
        self._b = b
        self._c = c
        self._k = k
        self._plus = plus

        ziRight = rho * 0.99999
        ziLeft = -rho * 0.99999
        self._rightTail = Extrema(ziRight, self.getValue(ziRight), self.getPPrime(ziRight))
        self._leftTail = Extrema(ziLeft, self.getValue(ziLeft), self.getPPrime(ziLeft))

        self._updateExtrema()
        self._updateZeros()
        self._updateIntersections()

    def _getZj(self, zi):
        return self._plus * sqrt(self._rho*self._rho - zi*zi)

    def getValue(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)

        term1 = self._k[0] * zi * zi
        term2 = self._k[1] * zj * zj
        term3 = self._k[2] * self._zk * self._zk
        term4 = self._k[3] * zi * zj
        term5 = self._k[4] * zi * self._zk
        term6 = self._k[5] * zj * self._zk
        term7 = self._k[6] * zi
        term8 = self._k[7] * zj
        term9 = self._k[8] * self._zk

        return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 - self._k[9]

    def getPrime(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)
        zjP = -zi / zj
        return (2 * self._k[0] * zi) + (2 * self._k[1] * zj * zjP) + (self._k[3] * (zi * zjP + zj)) + (self._k[4]
                                        * self._zk) + (self._k[5] * self._zk * zjP) + (self._k[6]) + self._k[7] * zjP

    def getPPrime(self, zi, zj=None):
        if zj is None:
            zj = self._getZj(zi)
        zjP = -zi / zj
        zjPP = (zjP * zi - zj) / (zj * zj)
        return (2 * self._k[0]) + (2 * self._k[1] * (zj * zjPP + zjP * zjP)) + (self._k[3] * (zi * zjPP + 2 * zjP)) + (
                    self._k[5] * self._zk * zjPP) + (self._k[7] * zjPP)

    def plot(self):
        m = int(self._rho * .99999 * 1000)
        x = [i / 1000 for i in range(-m, m)]
        f = []
        fPrime = []
        fPPrime = []

        for zi in x:
            zj = self._getZj(zi)

            f.append(self.getValue(zi, zj))
            fPrime.append(self.getPrime(zi, zj))
            fPPrime.append(self.getPPrime(zi, zj))

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

    def _computeExtremaValues(self, start):
        assert start is True or start is False or start == 0

        if start is True:
            zi1 = self._rho * 0.99999
            zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            ziP = getZiP(zi1, zj, self._zk, self._k)
            ziPP = getZiPP(zi1, zj, self._zk, self._k)
            # reduce initial guess until it is either positive and increasing or negative and decreasing
            while ziP > 0 and ziPP < 0 or ziP < 0 and ziPP > 0:
                zi1 *= 0.99
                ziP = getZiP(zi1, zj, self._zk, self._k)
                ziPP = getZiPP(zi1, zj, self._zk, self._k)
        elif start is False:
            zi1 = -self._rho * 0.99999
            zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            ziP = getZiP(zi1, zj, self._zk, self._k)
            ziPP = getZiPP(zi1, zj, self._zk, self._k)
            # reduce initial guess until it is either positive and decreasing or negative and increasing
            while ziP > 0 and ziPP > 0 or ziP < 0 and ziPP < 0:
                zi1 *= 0.99
                ziP = getZiP(zi1, zj, self._zk, self._k)
                ziPP = getZiPP(zi1, zj, self._zk, self._k)
        else:
            zi1 = 0
        zi0 = -2

        while abs(zi0 - zi1) > 1e-7:
            zi0 = zi1
            try:
                zj = self._plus * sqrt(self._rho * self._rho - zi1 * zi1)
            except ValueError:
                return None, None, None
            zi1 = zi0 - getZiP(zi0, zj, self._zk, self._k) / getZiPP(zi0, zj, self._zk, self._k)

        return zi1, getZi(zi1, self._zk, self._rho, self._k, self._plus), getZiPP(zi1, zj, self._zk, self._k)

    def _updateTails(self):
        ziRight = self._rho * 0.99999
        ziLeft = -self._rho * 0.99999
        self._rightTail = Extrema(ziRight, self.getValue(ziRight), self.getPPrime(ziRight))
        self._leftTail = Extrema(ziLeft, self.getValue(ziLeft), self.getPPrime(ziLeft))

    def _updateExtrema(self):
        self._updateTails()
        rightArgs = self._computeExtremaValues(True)
        right = Extrema(*rightArgs) if rightArgs[0] is not None else None
        zeroArgs = self._computeExtremaValues(0)
        zero = Extrema(*zeroArgs) if zeroArgs[0] is not None else None
        leftArgs = self._computeExtremaValues(False)
        left = Extrema(*leftArgs) if leftArgs[0] is not None else None

        xt = []
        for i in (right, zero, left):
            if i is not None:
                if all([i != x for x in xt]):
                    xt.append(i)

        self._extrema = xt

    def _updateZeros(self):
        sortedExtrema = sorted([self._rightTail, self._leftTail] + self._extrema, key=lambda o: o.zi)
        zeros = []
        for i, x in enumerate(sortedExtrema[:-1]):
            if x.value * sortedExtrema[i+1].value < 0:
                z = x.zi, sortedExtrema[i+1].zi
                zeros.append(z)
        self._zeros = zeros
        self._zeroCount = len(zeros)

    def _updateIntersections(self):
        intersections = []
        for z in self._zeros:
            zi0 = -2
            if self.getPrime(z[0]) > 1e-3:
                zi1 = z[0]
            elif self.getPrime(z[1]) > 1e-3:
                zi1 = z[1]
            else:
                zi1 = (z[0] + z[1]) / 2
            while abs(zi0 - zi1) > 1e-7:
                zi0 = zi1
                zi1 = zi0 - self.getValue(zi0) / self.getPrime(zi0)
            intersections.append(zi1)

        self._intersections = intersections

    @property
    def extrema(self):
        return self._extrema

    def updateSatValues(self, a, b, c, k):
        self._a = a
        self._b = b
        self._c = c
        self._k = k

    def updateGeoValues(self, zk, rho):
        self._zk = zk
        self._rho = rho

    def update(self):
        self._updateExtrema()
        self._updateZeros()
        self._updateIntersections()

    @property
    def zeroCount(self):
        return self._zeroCount

    @property
    def zeros(self):
        return self._zeros

    @property
    def intersections(self):
        return self._intersections

    def updateSpecific(self, a, b, c, k, zi):
        self._a = a
        self._b = b
        self._c = c
        self._k = k


class OrbitPath2:

    def __init__(self, sat, geo, jd):
        self._sat = sat
        self._geo = geo
        self._jd = jd

        self._zk = geo.getZenithVector(jd)[2]
        self._rho = cos(radians(geo.latitude))
        self._mask = -1     # start in error state

    def _getEllipseVectors(self, jd):
        elements = self._sat.getElements(jd)
        mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
        amag = elements.sma
        e = elements.ecc
        cmag = amag * e
        bmag = amag * sqrt(1 - e * e)
        a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
        b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
        c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

        return a, b, c

    def _getConstants(self, a, b, c, jd):
        # a, b, c = getEllipseVectors(self._sat, jd)
        gamma = geo.getPositionVector(jd)
        phiPrime = geo.getGeocentricLatitude()

        d = gamma.mag() * cos(radians(geo.latitude) - phiPrime)
        k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
        k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
        k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
        k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
        k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
        k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
        k7 = 2 * d * c[0]
        k8 = 2 * d * c[1]
        k9 = 2 * d * c[2]
        k10 = d * d

        # return a, b, c, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]
        return [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]

    def _getZj(self, zi, plus):
        assert plus == 1 or plus == -1
        return plus * sqrt(self._rho*self._rho - zi*zi)

    def _getGeneralZeros(self, jd):
        # zk = self._geo.getZenithVector(jd)
        # rho = cos(radians(geo.latitude))
        eVecs = self._getEllipseVectors(jd)
        k = self._getConstants(*eVecs, jd)

        posFunc = ZeroFunction(self._zk, self._rho, k, 1)
        negFunc = ZeroFunction(self._zk, self._rho, k, -1)

        # todo: do we need to worry about duplicates in both the positive and negative functions? gut says no
        zeros = [(z, self._getZj(z, 1), 1) for z in posFunc.computeGeneralZeros()] \
                + [(z, self._getZj(z, -1), -1) for z in negFunc.computeGeneralZeros()]
        return zeros

    def _angleDifference(self, angle):
        if angle > pi:
            return angle - TWOPI
        elif angle < -pi:
            return angle + TWOPI
        return angle

    def _getTimeTo(self, zi, zj, jd, next):
        lng = atan2(zj, zi) - earthOffsetAngle(jd)
        dl = self._angleDifference(lng - radians(geo.longitude))
        if next and dl < 0:
            dl += TWOPI
        elif not next and dl > 0:
            dl -= TWOPI

        dt = SIDEREAL_PER_SOLAR * dl / TWOPI
        return jd.future(dt)

    def _checkTime(self, jd):
        a, b, c = self._getEllipseVectors(jd)
        zeta = self._geo.getZenithVector(jd)
        gamma = self._geo.getPositionVector(jd)

        num = dot(zeta, gamma - c)
        den = sqrt(dot(zeta, a)**2 + dot(zeta, b)**2)

        return num / den

    def _refineZero(self, zi, jd, plus, next):
        assert plus == 1 or plus == -1
        # zj = plus * sqrt(self._rho*self._rho - zi*zi)
        eVecs = self._getEllipseVectors(jd)
        k = self._getConstants(*eVecs, jd)
        func = ZeroFunction(self._zk, self._rho, k, plus)

        check = 2
        time = jd
        while not (0.9999999 < check <= 1):
            eVecs = self._getEllipseVectors(time)
            k = self._getConstants(*eVecs, time)
            zi = func.computeSpecificZero(zi, k)
            zj = plus * sqrt(self._rho*self._rho - zi*zi)
            time = self._getTimeTo(zi, zj, jd, next)
            check = self._checkTime(time)

        return time

    def computeIntersectionTimes(self, jd):
        zeros = self._getGeneralZeros(jd)
        lenZeros = len(zeros)
        if lenZeros == 0:
            return []
        elif lenZeros != 2 and lenZeros != 4:
            raise ValueError('expected to find 0, 2, or 4 intersections, found %i' % lenZeros)

        lngs = [atan2(zj, zi) - earthOffsetAngle(jd) for zi, zj, _ in zeros]
        geoLng = radians(self._geo.longitude)
        dls = [self._angleDifference(lng - geoLng) for lng in lngs]

        zeroSorted = [z for _, z in sorted(zip(dls, zeros))]
        # dlSorted = sorted(dls)

        altTerm = self._checkTime(jd)
        if -1 < altTerm < 1:
            if lenZeros == 2:
                direction = [False, True]
            else:
                # dlSorted = dlSorted[1:] + [dlSorted[0] + TWOPI]
                zeroSorted = zeroSorted[1:] + [zeroSorted[0]]
                direction = [False, True, True, True]
        else:
            halfLen = lenZeros >> 1
            # dlSorted = dlSorted[halfLen:] + [i + TWOPI for i in dlSorted[:halfLen]]
            zeroSorted = zeroSorted[halfLen:] + zeroSorted[:halfLen]
            direction = [True] * lenZeros

        times = [self._refineZero(zi, jd, p, d) for (zi, _, p), d in zip(zeroSorted, direction)]

        return times


INTERSECTION_EXISTS = 0b001  # 1 yes, 0 no
INTERSECTION_COUNT = 0b010  # 1 - 4, 0 - 2
UNDER_ORBIT = 0b100  # 1 - yes, 0 no
TWO_INTERSECTIONS = INTERSECTION_EXISTS & ~INTERSECTION_COUNT
FOUR_INTERSECTIONS = INTERSECTION_EXISTS | INTERSECTION_COUNT


class OrbitPath:

    def __init__(self, sat, geo, jd):
        self._sat = sat
        self._geo = geo
        self._jd = jd

        a, b, c, k = self._getConstants(jd)
        zeta = geo.getZenithVector(jd)
        self._rho = cos(geo.getGeocentricLatitude())

        self._posFunction = Function(zeta[2], self._rho, a, b, c, k, 1)
        self._negFunction = Function(zeta[2], self._rho, a, b, c, k, -1)

        # posIntersections = [Intersection(zi, zeta[2], 1) for zi in self._posFunction.intersections]
        # negIntersections = [Intersection(zi, zeta[2], -1) for zi in self._negFunction.intersections]
        # self._intersections = posIntersections + negIntersections

        self.computeTimes(jd)

    # @staticmethod
    # def _orderIntersections(intersections, lng):
    #     pass

    def _getConstants(self, jd):
        a, b, c = getEllipseVectors(self._sat, jd)
        gamma = geo.getPositionVector(jd)
        phiPrime = geo.getGeocentricLatitude()

        d = gamma.mag() * cos(radians(geo.latitude) - phiPrime)
        k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
        k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
        k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
        k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
        k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
        k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
        k7 = 2 * d * c[0]
        k8 = 2 * d * c[1]
        k9 = 2 * d * c[2]
        k10 = d * d

        return a, b, c, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]

    def _getEllipseVectors(self):
        elements = self._sat.getElements(self._jd)
        mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
        amag = elements.sma
        e = elements.ecc
        cmag = amag * e
        bmag = amag * sqrt(1 - e * e)
        a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
        b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
        c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

        return a, b, c

    @property
    def intersections(self):
        return self._posFunction.intersections + self._negFunction.intersections

    @property
    def positive(self):
        return self._posFunction

    @property
    def negative(self):
        return self._negFunction

    @staticmethod
    def _orderIntersections(lngs, ints):
        intsSorted = [x for _, x in zip(lngs, ints)]
        lngsSorted = sorted(lngs)
        return lngsSorted, intsSorted

    @staticmethod
    def _angleDifference(diff):
        if diff > pi:
            return diff - TWOPI
        elif diff < -pi:
            return diff + TWOPI
        return diff

    def _timeTo(self, dLng, jd0):
        dt = dLng / TWOPI
        # convert solar day to sidereal day
        return jd0.future(dt * SIDEREAL_PER_SOLAR)

    def computeTimes(self, jd):
        mask = 0
        rawZis = self._posFunction.intersections + self._negFunction.intersections
        if len(rawZis) == 0:
            return []
        elif len(rawZis) == 2:
            mask = TWO_INTERSECTIONS
        elif len(rawZis) == 4:
            mask = FOUR_INTERSECTIONS
        else:
            raise ValueError('expected to find 0, 2 or 4 intersections, found %i' % len(rawZis))
        zjSigns = [1 for i in range(len(self._posFunction.intersections))] + [-1 for i in range(len(self._negFunction.intersections))]
        zjs = [zjSign * sqrt(self._rho*self._rho - zi*zi) for zi, zjSign in zip(rawZis, zjSigns)]

        # zjs = []
        # for i, zi in enumerate(rawZis):
        #     zj = zjSigns[i] * sqrt(self._rho*self._rho - zi*zi)
        #     zjs.append(zj)

        geoLng = radians(self._geo.longitude)
        lngs = [atan2(zj, zi) for zi, zj in zip(rawZis, zjs)]
        dls = [self._angleDifference(lng - geoLng) for lng in lngs]

        dlSorted, ziSorted = self._orderIntersections(dls, rawZis)

        # correctly 'order' the sorted dl's, based on wanting future or past time(s).
        # if orbit path is above, find most recent negative dl, else find next 2 positive

        a, b, c, k = self._getConstants(jd)
        zeta = geo.getZenithVector(jd)
        gamma = geo.getPositionVector(jd)
        altTerm = dot(zeta, gamma - c) / sqrt(dot(zeta, a)**2 + dot(zeta, b)**2)
        if -1 < altTerm < 1:
            mask |= UNDER_ORBIT
            # find previous most recent negative dl, next for the rest
            if mask & FOUR_INTERSECTIONS:
            # if len(dlSorted) == 4:
                dlSorted = dlSorted[1:] + [dlSorted[0] + TWOPI]
                ziSorted = ziSorted[1:] + [ziSorted[0]]
            # if two intersections, should be ordered how we want it
        else:
            halfLen = int(len(dlSorted) / 2)
            dlSorted = dlSorted[halfLen:] + [i + TWOPI for i in dlSorted[:halfLen]]
            ziSorted = ziSorted[halfLen:] + ziSorted[:halfLen]

        times = [self._timeTo(dl, jd) for dl in dlSorted]

        stopRecurse = []
        for i, dl, zi, jd in enumerate(zip(dlSorted, ziSorted, times)):
            if i in stopRecurse:
                continue

            a, b, c, k = self._getConstants(jd)


class Plot:

    def __init__(self, sat, geo, jd):
        self._sat = sat
        self._geo = geo
        self._jd = jd

        a, b, c, k = self._getConstants()
        zeta = geo.getZenithVector(jd)
        rho = cos(geo.getGeocentricLatitude())

        self._posFunction = Function(zeta[2], rho, a, b, c, k, 1)
        self._negFunction = Function(zeta[2], rho, a, b, c, k, -1)

    def _getConstants(self):
        a, b, c = getEllipseVectors(self._sat, self._jd)
        gamma = geo.getPositionVector(jd)
        phiPrime = geo.getGeocentricLatitude()

        d = gamma.mag() * cos(radians(geo.latitude) - phiPrime)
        k1 = a[0] * a[0] + b[0] * b[0] - c[0] * c[0]
        k2 = a[1] * a[1] + b[1] * b[1] - c[1] * c[1]
        k3 = a[2] * a[2] + b[2] * b[2] - c[2] * c[2]
        k4 = 2 * (a[0] * a[1] + b[0] * b[1] - c[0] * c[1])
        k5 = 2 * (a[0] * a[2] + b[0] * b[2] - c[0] * c[2])
        k6 = 2 * (a[1] * a[2] + b[1] * b[2] - c[1] * c[2])
        k7 = 2 * d * c[0]
        k8 = 2 * d * c[1]
        k9 = 2 * d * c[2]
        k10 = d * d

        return a, b, c, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10]

    def _getEllipseVectors(self):
        elements = self._sat.getElements(self._jd)
        mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
        amag = elements.sma
        e = elements.ecc
        cmag = amag * e
        bmag = amag * sqrt(1 - e * e)
        a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
        b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
        c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

        return a, b, c

    @property
    def intersections(self):
        return self._posFunction.intersections + self._negFunction.intersections

    @property
    def positive(self):
        return self._posFunction

    @property
    def negative(self):
        return self._negFunction


SIDEREAL_PER_SOLAR = EARTH_SIDEREAL_PERIOD / 86400


class Intersection:

    def __init__(self, zi, plus):
        self._zi = zi
        # self._rho = sqrt(1 - zk*zk)
        assert plus == 1 or plus == -1
        self._plus = plus

    @property
    def zi(self):
        return self._zi

    @property
    def plus(self):
        return self._plus

    # def update(self, zero):
    #     self._zi = zero

    # def computeTimeTo(self, zi, geoLng, jd, next=True):
    #     '''next is True to find time to next lng pass, False is previous'''
    #     zj = self._plus * sqrt(self._rho*self._rho - zi*zi)
    #
    #     lng = (atan2(zj, zi) - earthOffsetAngle(jd)) % TWOPI
    #     if lng > pi:
    #         lng -= TWOPI
    #
    #     dLng = lng - geoLng
    #     if dLng > pi:
    #         dLng -= TWOPI
    #     elif dLng < -pi:
    #         dLng += TWOPI
    #
    #     if next is True and dLng < 0:
    #         dLng += TWOPI
    #     dt = SIDEREAL_PER_SOLAR * dLng / 360
    #
    #     return jd.future(dt)


def makeZetaPlot(sat, geo, jd, plus=True, ylim=None, filename=None):
    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)
    elements = sat.getElements(jd)
    mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    amag = elements.sma
    e = elements.ecc
    cmag = amag * e
    bmag = amag * sqrt(1 - e * e)
    a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
    b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
    c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

    rho = cos(geo.getGeocentricLatitude())
    d = dot(zeta, gamma)

    k1 = a[0]*a[0] + b[0]*b[0] - c[0]*c[0]
    k2 = a[1]*a[1] + b[1]*b[1] - c[1]*c[1]
    k3 = a[2]*a[2] + b[2]*b[2] - c[2]*c[2]
    k4 = 2 * (a[0]*a[1] + b[0]*b[1] - c[0]*c[1])
    k5 = 2 * (a[0]*a[2] + b[0]*b[2] - c[0]*c[2])
    k6 = 2 * (a[1]*a[2] + b[1]*b[2] - c[1]*c[2])
    k7 = 2 * d * c[0]
    k8 = 2 * d * c[1]
    k9 = 2 * d * c[2]
    k10 = d * d

    limit = cos(geo.getGeocentricLatitude())
    m = int(limit * 1000)
    x = [i/1000 for i in range(-m, m)]

    positiveArr = []
    negativeArr = []
    positiveArr2 = []
    negativeArr2 = []
    for i in x:
        fpos = getZi(i, zeta[2], rho, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10], 1)
        fneg = getZi(i, zeta[2], rho, [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10], -1)
        positiveArr.append(fpos)
        negativeArr.append(fneg)
        # print('zk:', zeta[2])
        # print('a:', a)
        # print('b:', b)
        # print('c:', c)
        # print('zi:', i)
        # print('fpos:', fpos)
        # print('fneg:', fneg)
        # return
        positiveArr2.append(getZi2(i, zeta[2], rho, b, 1))
        negativeArr2.append(getZi2(i, zeta[2], rho, b, -1))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, positiveArr, s=3, c='lime', label='positive zero function')
    ax.scatter(x, negativeArr, s=3, c='darkorange', label='negative zero function')
    # ax.scatter(x, positiveArr2, s=3)
    # ax.scatter(x, negativeArr2, s=3)

    # fig, ax = plt.subplots()
    # ax.plot(x, ans)
    ax.set(xlabel='$^\zeta$i', ylabel='f($^\zeta$i)',
           title=f'Zero function for {sat.name.replace(" ", "-")} at {geo.latitude}$^\circ$ latitude.')
    ax.grid()
    if ylim:
        plt.ylim(ylim)
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)


def makeZetaPlotTrue(sat, geo, jd, plus=True, ylim=None, filename=None):
    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)
    elements = sat.getElements(jd)
    mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    amag = elements.sma
    e = elements.ecc
    cmag = amag * e
    bmag = amag * sqrt(1 - e * e)
    a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
    b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
    c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

    rho = cos(geo.getGeocentricLatitude())
    d = dot(zeta, gamma)

    limit = cos(geo.getGeocentricLatitude())
    m = int(limit * 1000)
    x = [i/1000 for i in range(-m, m)]

    def f(zi, sign):
        zj = sign * sqrt(rho*rho - zi*zi)
        z = Vector(zi, zj, zeta[2])
        theta = atan2(zj, zi)
        rho2 = gamma.mag() * rho
        gi = rho2 * cos(theta)
        gj = rho2 * sin(theta)
        g = Vector(gi, gj, gamma[2])

        num = dot(z, g-c)
        den = sqrt(dot(z, a)**2 + dot(z, b)**2)
        return num / den - 1

    positiveArr = []
    negativeArr = []
    positiveArr2 = []
    negativeArr2 = []
    for i in x:

        positiveArr.append(f(i, 1))
        negativeArr.append(f(i, -1))
        positiveArr2.append(getZi2(i, zeta[2], rho, b, 1))
        negativeArr2.append(getZi2(i, zeta[2], rho, b, -1))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, positiveArr, s=3, c='lime', label='positive zero function')
    ax.scatter(x, negativeArr, s=3, c='darkorange', label='negative zero function')
    ax.scatter(x, positiveArr2, s=3, c='red', label='positive zero function2')
    ax.scatter(x, negativeArr2, s=3, c='blue', label='negative zero function2')

    # fig, ax = plt.subplots()
    # ax.plot(x, ans)
    ax.set(xlabel='$^\zeta$i', ylabel='f($^\zeta$i)',
           title=f'Zero function for {sat.name.replace(" ", "-")} at {geo.latitude}$^\circ$ latitude.')
    ax.grid()
    if ylim:
        plt.ylim(ylim)
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)


def complexPlot(sat, geo, jd):
    cosArray = []
    tanArray = []
    angleArrayPos = []
    angleArrayNeg = []
    x = [i / 1440 for i in range(1440)]

    for i in x:
        time = jd.future(i)
        a, b, c = getEllipseVectors(sat, time)
        zeta = geo.getZenithVector(time)
        gamma = geo.getPositionVector(time)

        cosNum = dot(zeta, gamma - c)
        cosDen = sqrt(dot(zeta, a) ** 2 + dot(zeta, b) ** 2)
        # cosTerm = cosNum / cosDen

        tanNum = dot(zeta, b)
        tanDen = dot(zeta, a)
        tanTerm = atan2(tanNum, tanDen)

        # cosArray.append(cosTerm)
        if cosNum / cosDen <= 1:
        # if cosTerm <= 1:
            cosTerm = acos(cosNum / cosDen)
            pos = tanTerm + acos(cosTerm)
            if pos < 0:
                pos += TWOPI
            neg = tanTerm - acos(cosTerm)
            if neg < 0:
                neg += TWOPI
            angleArrayPos.append(pos)
            angleArrayNeg.append(neg)
        else:
            angleArrayPos.append(-10.0)
            angleArrayNeg.append(-10.0)
            cosTerm = -10
            tanTerm = -10
        cosArray.append(cosTerm)
        tanArray.append(tanTerm)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, tanArray, s=2, c='darkorange', label='atan term')
    ax.scatter(x, angleArrayPos, s=2, c='red', label='acos + atan')
    ax.scatter(x, angleArrayNeg, s=2, c='blue', label='-acos + atan')

    axCos = ax.twinx()
    axCos.scatter(x, cosArray, s=2, c='lime', label='acos term')
    axCos.set(ylim=(-1, 1))

    ax.set(xlabel='dt (solar days)', ylabel='angle (radians)',
           title=f'Angle terms for {sat.name.replace(" ", "-")} at {geo.latitude}$^\circ$ latitude.', ylim=(-7, 7))
    ax.grid()
    plt.legend()
    plt.show()


def complexPlot2(sat, geo, jd):
    x = [i / 1440 for i in range(1440)]
    cosNumArr = []
    cosDenArr = []
    altArr = []

    for i in x:
        time = jd.future(i)
        a, b, c = getEllipseVectors(sat, time)
        zeta = geo.getZenithVector(time)
        gamma = geo.getPositionVector(time)

        cosNum = dot(zeta, gamma - c)
        cosDen = sqrt(dot(zeta, a) ** 2 + dot(zeta, b) ** 2)

        # alt = _orbitAltitude(sat, geo, time)
        # altArr.append(alt)
        alt = getAltitudeBisect(sat, geo, time)
        altArr.append(alt)

        cosNumArr.append(cosNum)
        cosDenArr.append(cosDen)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, cosNumArr, s=2, c='lime', label='acos numerator term')
    ax.scatter(x, cosDenArr, s=2, c='darkorange', label='acos denominator term')

    ax2 = ax.twinx()
    ax2.scatter(x, altArr, s=2, c='blue', label='orbital plane altitude')

    ax.set(xlabel='dt (solar days)', ylabel='angle (radians)',
           title=f'Angle terms for {sat.name.replace(" ", "-")} at {geo.latitude}$^\circ$ latitude.')
    ax.grid()
    # ax2.set(ylabel='altitude (radians)')
    # ax2.grid()
    plt.legend()
    plt.show()


def getEllipseVectors(sat, jd):
    elements = sat.getElements(jd)
    mat = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
    amag = elements.sma
    e = elements.ecc
    cmag = amag * e
    bmag = amag * sqrt(1 - e * e)
    a = rotateMatrixFrom(mat, Vector(amag, 0, 0))
    b = rotateMatrixFrom(mat, Vector(0, bmag, 0))
    c = rotateMatrixFrom(mat, Vector(-cmag, 0, 0))

    return a, b, c


def getAltitudeFuncValues(sat, geo, jd, t):
    a, b, c = getEllipseVectors(sat, jd)

    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)

    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma

    rhoMag = rho.mag()
    rhoP = b * cos(t) - a * sin(t)
    alpha = dot(rho, zeta) / rhoMag
    sqt = sqrt(1 - alpha * alpha)

    term1 = rhoP / (rhoMag * sqt)
    term2 = rho * dot(rho, rhoP) / (rhoMag * rhoMag * rhoMag * sqt)
    f = dot(zeta, term1 - term2)

    rhoMag2 = rho.mag2()
    rhoPP = -a * cos(t) - b * sin(t)
    rhoMagP = dot(rho, rhoP) / rhoMag
    alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP)

    mP = -gamma * ((rhoMag * -alphaP) / (2 * sqt) + rhoMagP * sqt) / (rhoMag * rhoMag * (1 - alpha * alpha))
    term1 = rhoP * dot(rho, rhoP) + rho * (dot(rho, rhoPP) + dot(rhoP, rhoP))
    term2 = rho * rhoMagP * dot(rho, rhoP)
    term3 = rho * dot(rho, rhoP) / rhoMag
    nP = (term1 / rhoMag) - (term2 / rhoMag2)
    fp = dot(zeta / (rhoMag * sqt), -nP) + dot(mP, rhoP - term3)

    return f, fp


def getAltitudeFuncValues2(a, b, c, zeta, gamma, t):
    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma

    rhoMag = rho.mag()
    rhoP = b * cos(t) - a * sin(t)
    alpha = dot(rho, zeta) / rhoMag
    sqt = sqrt(1 - alpha * alpha)

    term1 = rhoP / (rhoMag * sqt)
    term2 = rho * dot(rho, rhoP) / (rhoMag * rhoMag * rhoMag * sqt)
    f = dot(zeta, term1 - term2)

    rhoMag2 = rho.mag2()
    rhoPP = -a * cos(t) - b * sin(t)
    rhoMagP = dot(rho, rhoP) / rhoMag
    alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP)

    mP = -gamma * ((rhoMag * -alphaP) / (2 * sqt) + rhoMagP * sqt) / (rhoMag * rhoMag * (1 - alpha * alpha))
    term1 = rhoP * dot(rho, rhoP) + rho * (dot(rho, rhoPP) + dot(rhoP, rhoP))
    term2 = rho * rhoMagP * dot(rho, rhoP)
    term3 = rho * dot(rho, rhoP) / rhoMag
    nP = (term1 / rhoMag) - (term2 / rhoMag2)
    fp = dot(zeta / (rhoMag * sqt), -nP) + dot(mP, rhoP - term3)

    return f, fp


def getAltitudeFuncValues3(a, b, c, zeta, gamma, t):
    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma

    rhoMag = rho.mag()
    rhoMag2 = rho.mag2()
    alpha = dot(rho, zeta) / rhoMag
    sqt = sqrt(1 - alpha * alpha)
    rhoP = b * cos(t) - a * sin(t)
    rhoMagP = dot(rho, rhoP) / rhoMag

    alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta)) / rhoMag2
    f = alphaP / sqt

    rhoPP = -a * cos(t) - b * sin(t)
    rhoMagPP = (rhoMag * (dot(rho, rhoPP) + dot(rhoP, rhoP)) - rhoMagP * dot(rho, rhoP)) / rhoMag2
    aPPterm1 = rhoMag * dot(rhoPP, zeta) + rhoMagPP * dot(rho, zeta)
    aPPterm2 = 2 * rhoMag * rhoMagP * (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta))
    alphaPP = (rhoMag2 * aPPterm1 - aPPterm2) / (rhoMag2 * rhoMag2)

    term1 = alphaPP * sqt
    term2 = alpha * alphaP * alphaP / sqt
    fP = (term1 + term2) / (1 - alpha * alpha)

    return f, fP


def getAltitudeOld(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)
    Fs = [abs(getAltitudeFuncValues2(a, b, c, zeta, gamma, t/100)[0]) for t in range(628)]
    t = Fs.index(min(Fs)) / 100
    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma
    if dot(rho, zeta) < 0:
        t += pi
        r = c + a * cos(t) + b * sin(t)
        rho = r - gamma
    return asin(dot(rho, zeta) / rho.mag())
    # t0 = -1
    # t1 = 0
    #
    # while abs(t1 - t0) > 1e-5:
    #     t0 = t1
    #     f, fp = getAltitudeFuncValues2(a, b, c, zeta, gamma, t1)
    #     t1 = t0 - f / fp
    #
    # t1 %= TWOPI
    # r = c + a * cos(t1) + b * sin(t1)
    # rho = r - gamma
    # if dot(rho, zeta) < 0:
    #     t0 = -1
    #     t1 += pi
    #     while abs(t1 - t0) > 1e-5:
    #         t0 = t1
    #         f, fp = getAltitudeFuncValues2(a, b, c, zeta, gamma, t1)
    #         t1 = t0 - f / fp
    #
    # r = c + a * cos(t1) + b * sin(t1)
    # rho = r - gamma
    # return asin(dot(rho, zeta) / rho.mag())


def makeAltitudeFPlot(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)

    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)

    x = [i / 100 for i in range(628)]
    F = []
    FP = []

    for t in x:
        r = c + a * cos(t) + b * sin(t)
        rho = r - gamma

        rhoP = b * cos(t) - a * sin(t)
        rhoMag = rho.mag()
        alpha = dot(rho, zeta) / rhoMag
        rhoMag2 = rho.mag2()
        sqt = sqrt(1 - alpha * alpha)
        rhoPP = -a * cos(t) - b * sin(t)
        rhoMagP = dot(rho, rhoP) / rhoMag
        rhoMagPP = (rhoMag * (dot(rho, rhoPP) + dot(rhoP, rhoP)) - rhoMagP * dot(rho, rhoP)) / rhoMag2

        m = zeta / rhoMag2 * sqt
        n = rhoP * rhoMag - rho * rhoMagP
        f = dot(m, n)

        alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta)) / rhoMag2
        term1 = (-rhoMag2 * alpha * alphaP) / sqt
        # term1 = -rhoMag2 * rhoMag2 * alpha * f
        term2 = 2 * sqt * rhoMag * rhoMagP
        term3 = rhoMag2 * rhoMag2 * (1 - alpha * alpha)
        mP = -zeta * (term1 + term2) / term3
        nP = rhoPP * rhoMag - rho * rhoMagPP
        fP = dot(m, nP) + dot(mP, n)

        # rhoMag = rho.mag()
        # rhoP = b * cos(t) - a * sin(t)
        # alpha = dot(rho, zeta) / rhoMag
        # sqt = sqrt(1 - alpha * alpha)
        #
        # term1 = rhoP / (rhoMag * sqt)
        # term2 = rho * dot(rho, rhoP) / (rhoMag * rhoMag * rhoMag * sqt)
        # f = dot(zeta, term1 - term2)
        #
        # rhoMag2 = rho.mag2()
        # rhoPP = -a * cos(t) - b * sin(t)
        # rhoMagP = dot(rho, rhoP) / rhoMag
        # alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP)
        #
        # mP = -gamma * ((rhoMag * -alphaP) / (2 * sqt) + rhoMagP * sqt) / (rhoMag * rhoMag * (1 - alpha * alpha))
        # term1 = rhoP * dot(rho, rhoP) + rho * (dot(rho, rhoPP) + dot(rhoP, rhoP))
        # term2 = rho * rhoMagP * dot(rho, rhoP)
        # term3 = rho * dot(rho, rhoP) / rhoMag
        # nP = (term1 / rhoMag) - (term2 / rhoMag2)
        # fp = dot(zeta / (rhoMag * sqt), -nP) + dot(mP, rhoP - term3)

        F.append(f)
        FP.append(fP)

    trueFP = [0]
    dt = 1 / 100
    for i in range(len(F) - 1):
        dF = F[i+1] - F[i]
        trueFP.append(dF / dt)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, F, '-', c='blue')
    ax.plot(x, trueFP, '-', c='darkorange')
    ax2 = ax.twinx()
    ax2.plot(x, FP, '-', c='red')
    ax.grid()
    plt.show()


def makePArr(arr, dt = 1/100):
    trueVal = [0]
    for i in (range(627)):
        dVal = arr[i+1] - arr[i]
        trueVal.append(dVal / dt)
    return trueVal


def plotArrs(x, arr, trueP, arrP):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, arr, '-', c='blue')
    ax.plot(x, trueP, '-', c='darkorange')
    ax.plot(x, arrP, '-', c='red')

    ax.grid()
    plt.show()


def makeAltitudeFPlot2(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)

    zeta = geo.getZenithVector(jd)
    gamma = geo.getPositionVector(jd)

    x = [i / 100 for i in range(628)]
    F = []
    FP = []

    for t in x:
        r = c + a * cos(t) + b * sin(t)
        rho = r - gamma

        rhoMag = rho.mag()
        rhoMag2 = rho.mag2()
        alpha = dot(rho, zeta) / rhoMag
        sqt = sqrt(1 - alpha * alpha)
        rhoP = b * cos(t) - a * sin(t)
        rhoMagP = dot(rho, rhoP) / rhoMag

        alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta)) / rhoMag2
        f = alphaP / sqt

        rhoPP = -a * cos(t) - b * sin(t)
        rhoMagPP = (rhoMag * (dot(rho, rhoPP) + dot(rhoP, rhoP)) - rhoMagP * dot(rho, rhoP)) / rhoMag2
        aPPterm1 = rhoMag * dot(rhoPP, zeta) - rhoMagPP * dot(rho, zeta)
        aPPterm2 = 2 * rhoMag * rhoMagP * (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta))
        alphaPP = (rhoMag2 * aPPterm1 - aPPterm2) / (rhoMag2 * rhoMag2)

        term1 = alphaPP * sqt
        term2 = alpha * alphaP * alphaP / sqt
        fP = (term1 + term2) / (1 - alpha * alpha)

        F.append(f)
        FP.append(fP)

    trueFP = []
    dt = 1e-6
    for i in range(628):
        dF = getAltitudeFuncValues3(a, b, c, zeta, gamma, i/100 + dt)[0] - getAltitudeFuncValues3(a, b, c, zeta, gamma, i/100)[0]
        trueFP.append(dF / dt)

    # trueFP = [0]
    # dt = 1 / 100
    # for i in range(len(F) - 1):
    #     dF = F[i+1] - F[i]
    #     trueFP.append(dF / dt)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(x, F, '-', c='blue')
    ax.plot(x, trueFP, '-', c='darkorange')
    ax.plot(x, FP, '-', c='red')
    # ax2 = ax.twinx()
    # ax2.plot(x, FP, '-', c='red')
    ax.grid()
    plt.show()


def getAltitudeFuncValues4(a, b, c, zeta, gamma, t):
    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma

    rhoMag = rho.mag()
    rhoMag2 = rho.mag2()
    alpha = dot(rho, zeta) / rhoMag
    sqt = sqrt(1 - alpha * alpha)
    rhoP = b * cos(t) - a * sin(t)
    rhoMagP = dot(rho, rhoP) / rhoMag

    alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta)) / rhoMag2
    f = alphaP / sqt

    rhoPP = -a * cos(t) - b * sin(t)
    rhoMagPP = (rhoMag * (dot(rho, rhoPP) + dot(rhoP, rhoP)) - rhoMagP * dot(rho, rhoP)) / rhoMag2
    aPPterm1 = rhoMag * dot(rhoPP, zeta) - rhoMagPP * dot(rho, zeta)
    aPPterm2 = 2 * rhoMag * rhoMagP * (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta))
    alphaPP = (rhoMag2 * aPPterm1 - aPPterm2) / (rhoMag2 * rhoMag2)

    term1 = alphaPP * sqt
    term2 = alpha * alphaP * alphaP / sqt
    fP = (term1 + term2) / (1 - alpha * alpha)

    return f, fP


def getAltitude(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)

    t0 = -1
    t1 = 0
    rho = 0

    # this needs to be a while loop because steep functions can cause drastic changes while checking the t1 + pi
    # rhoDotZeta = -1
    # while rhoDotZeta < 0:
    #     t1 += pi
    while abs(t1 - t0) > 1e-5:
        t0 = t1
        f, fp = getAltitudeFuncValues4(a, b, c, zeta, gamma, t0)
        t1 = t0 - f / fp

    t1 %= TWOPI
    t2 = t1 + pi

    r1 = c + a * cos(t1) + b * sin(t1)
    rho1 = r1 - gamma
    rhoDotZeta1 = dot(rho1, zeta)

    r2 = c + a * cos(t2) + b * sin(t2)
    rho2 = r2 - gamma
    rhoDotZeta2 = dot(rho2, zeta)

    if rhoDotZeta1 > rhoDotZeta2:
        return asin(rhoDotZeta1 / rho1.mag())
    else:
        return asin(rhoDotZeta2 / rho2.mag())


def getAltitudeFunc(a, b, c, zeta, gamma, t):
    r = c + a * cos(t) + b * sin(t)
    rho = r - gamma

    rhoMag = rho.mag()
    rhoMag2 = rho.mag2()
    alpha = dot(rho, zeta) / rhoMag
    sqt = sqrt(1 - alpha * alpha)
    rhoP = b * cos(t) - a * sin(t)
    rhoMagP = dot(rho, rhoP) / rhoMag

    alphaP = (rhoMag * dot(rhoP, zeta) - rhoMagP * dot(rho, zeta)) / rhoMag2
    f = alphaP / sqt
    return f


def getAltitudeSecant(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)

    t0 = 0
    t1 = pi
    f0 = getAltitudeFunc(a, b, c, zeta, gamma, t0)
    f1 = getAltitudeFunc(a, b, c, zeta, gamma, t1)
    t2 = 0

    while abs(t1 - t0) > 1e-5:
        t2 = t1 - f1 * (t1 - t0) / (f1 - f0)
        t0 = t1
        t1 = t2

    r = c + a * cos(t2) + b * sin(t2)
    rho = r - gamma
    return asin(dot(rho, zeta) / rho.mag())


def compareSign(a, b):
    return (a * b) > 0


def getSign(a):
    if a < 0:
        return -1
    else:
        return 1


def getAltitudeBisect(sat, geo, jd):
    a, b, c = getEllipseVectors(sat, jd)
    gamma = geo.getPositionVector(jd)
    zeta = geo.getZenithVector(jd)

    t0 = 0
    t1 = pi
    prevF0Sign = getSign(getAltitudeFunc(a, b, c, zeta, gamma, t0))
    prevF1Sign = getSign(getAltitudeFunc(a, b, c, zeta, gamma, t1))
    t = pi / 2
    ft = getAltitudeFunc(a, b, c, zeta, gamma, pi / 2)
    # tSign = getSign(ft)
    while abs(ft) > 0.01:
        t = (t0 + t1) / 2
        ft = getAltitudeFunc(a, b, c, zeta, gamma, t)
        tSign = getSign(ft)
        f0Sign = getSign(getAltitudeFunc(a, b, c, zeta, gamma, t0))
        if f0Sign == tSign:
            t0 = t
        else:
            t1 = t

        # f0Sign = getSign(getAltitudeFunc(a, b, c, zeta, gamma, t0))
        # f1Sign = getSign(getAltitudeFunc(a, b, c, zeta, gamma, t1))
        # if f0Sign == tSign:
        #     t0 = t
        # else:
        #     t1 = t
        # t = (t0 + t1) / 2
        # ft = getAltitudeFunc(a, b, c, zeta, gamma, t)
        # tSign = getSign(ft)

    r1 = c + a * cos(t)  + b * sin(t)
    rho1 = r1 - gamma
    r2 = c + a * cos(t + pi) + b * sin(t + pi)
    rho2 = r2 - gamma
    rdz1 = dot(rho1, zeta)
    rdz2 = dot(rho2, zeta)

    if rdz1 > rdz2:
        return asin(rdz1 / rho1.mag())
    else:
        return asin(rdz2 / rho2.mag())

    # return asin(dot(rho, zeta) / rho.mag())


def makeAltitudePlot(sat, geo, jd):
    x = [i / 1000 for i in range(1000)]
    # alts = [getAltitude(sat, geo, jd.future(i)) for i in x]
    alts = []
    for i in x:
        alts.append(getAltitudeBisect(sat, geo, jd.future(i)))
        print(i)
    # alts2 = [_orbitAltitude(sat, geo, jd.future(i)) for i in x]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, alts, '-', c='blue')
    # ax.plot(x, alts2, '-', c='red')
    ax.grid()
    plt.show()
