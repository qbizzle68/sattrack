from math import sin, cos, radians

from sattrack.util.constants import DELTAT, TWOPI
from sattrack.bodies._tables import SIN_TABLE, NUTATION_TABLE

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate


class JulianTimes:
    __slots__ = '_jd', '_jde', '_jc', '_jce', '_jme'

    def __init__(self, time: 'JulianDate', deltat=None):
        if deltat is None:
            deltat = DELTAT

        self._jd = time.value
        self._jde = time.value + deltat / 86400
        self._jc = (time.value - 2451545) / 36525
        self._jce = (self._jde - 2451545) / 36525
        self._jme = self._jce / 10

    @property
    def JD(self):
        return self._jd

    @property
    def JDE(self):
        return self._jde

    @property
    def JC(self):
        return self._jc

    @property
    def JCE(self):
        return self._jce

    @property
    def JME(self):
        return self._jme


def computeMeanElongationMoon(time: 'JulianDate') -> float:
    # returns in radians

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    return 5.1984694602504185 + JCE * (7771.377146170642 + JCE * (-3.34090925416754e-05 + (JCE * 9.21144458867353e-08)))


def computeMeanAnomalySun(time: 'JulianDate') -> float:
    # returns in radians

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    return 6.240036788719593 + JCE * (628.3019560241842 + JCE * (-2.79776279094691e-06 + (JCE * -5.81776417331443e-08)))


def computeMeanAnomalyMoon(time: 'JulianDate') -> float:
    # returns in radians

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    return 2.355548369303256 + JCE * (8328.691422882925 + JCE * (0.00015179477570445083 + (JCE * 3.10280755910103e-07)))


def computeLatitudeArgumentMoon(time: 'JulianDate') -> float:
    # returns in radians

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    return 1.6279019291238244 + JCE * (8433.466158317484 + JCE * (-6.42717497046911e-05 + (JCE * 5.33299493382934e-08)))


def computeRaanMoon(time: 'JulianDate') -> float:
    # returns in radians

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    return 2.1824385855759 + JCE * (-33.757045936662394 + JCE * (3.614227815029857e-05 + (JCE * 3.878509448876287e-08)))


def computeNutationDeltas(time: 'JulianDate') -> (float, float):
    # Returns nutation in longitude and nutation in obliquity in radians.

    if isinstance(time, JulianTimes):
        JCE = time.JCE
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525

    deltaPsi = 0
    deltaEpsilon = 0

    x0 = computeMeanElongationMoon(time)
    x1 = computeMeanAnomalySun(time)
    x2 = computeMeanAnomalyMoon(time)
    x3 = computeLatitudeArgumentMoon(time)
    x4 = computeRaanMoon(time)
    xTerms = (x0, x1, x2, x3, x4)

    for (a, b, c, d), yRow in zip(NUTATION_TABLE, SIN_TABLE):
        sinArg = 0
        cosArg = 0
        for xi, yi in zip(xTerms, yRow):
            sinArg += xi * yi
            cosArg += xi * yi

        dPsi_i = (a + b * JCE) * sin(sinArg)
        dEps_i = (c + d * JCE) * cos(cosArg)
        deltaPsi += dPsi_i
        deltaEpsilon += dEps_i

    return radians(deltaPsi / 36000000), radians(deltaEpsilon / 36000000)


def computeMeanObliquity(time: 'JulianDate') -> float:
    # Returns the mean obliquity of the ecliptic in radians.
    
    if isinstance(time, JulianTimes):
        JME = time.JME
    else:
        JDE = time.value + (DELTAT / 86400)
        JCE = (JDE - 2451545) / 36525
        JME = JCE / 10
    U = JME / 10

    epsilon0 = 23.439291111 + U * \
        (-1.300258333 + U *
            (-0.000430556 + U *
             (0.555347222 + U *
              (-0.01427222 + U *
               (-0.069352778 + U *
                (-0.010847222 + U *
                 (1.977778e3 + U *
                  (7.741667e3 + U *
                   (1.608333e3 + U * 6.805556e4)))))))))

    return radians(epsilon0)


def computeTrueObliquity(time: 'JulianDate') -> float:
    # Returns the true obliquity of the ecliptic in radians.

    meanObliquity = computeMeanObliquity(time)
    _, nutationObliquity = computeNutationDeltas(time)

    return meanObliquity + nutationObliquity


def computeMeanSiderealTime(time: 'JulianDate | JulianTimes') -> float:
    if isinstance(time, JulianTimes):
        JD = time.JD
        JC = time.JC
    else:
        JD = time.value
        JC = (JD - 2451545) / 36525

    rtn = 4.894961212735793 + 6.300388098984957 * (JD - 2451545) \
        + JC * JC * (6.770708127139162e-06 - JC * 4.508729661571505e-10)

    return rtn % TWOPI


def computeApparentSiderealTime(time: 'JulianDate | JulianTimes') -> float:
    deltaPsi, deltaEpsilon = computeNutationDeltas(time)
    meanObliquity = computeMeanObliquity(time)
    trueObliquity = meanObliquity + deltaEpsilon
    meanSiderealTime = computeMeanSiderealTime(time)

    return (meanSiderealTime + deltaPsi * cos(trueObliquity)) % TWOPI


def getEarthOffsetAngle(time: 'JulianDate') -> float:
    return computeApparentSiderealTime(time)
