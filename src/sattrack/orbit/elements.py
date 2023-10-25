from math import cos, sin, radians, degrees, sqrt, atan2, floor, pi
from typing import TYPE_CHECKING

from pyevspace import Vector

from sattrack.util.constants import TWOPI, EARTH_MU, SECONDS_PER_DAY
from sattrack.util.helpers import atan3
from sattrack.orbit.sgp4 import elementsFromState

if TYPE_CHECKING:
    from sattrack.orbit.sgp4 import TwoLineElement
    from sattrack.core.juliandate import JulianDate


# Using this method produces Earth shadow times almost identical to those produced by Stellarium, the sgp4 method,
# creates values several minutes off. Therefore, based on empirical evidence, this seems to be the correct methodology.
def elementsFromTle(tle: 'TwoLineElement', time: 'JulianDate') -> (float, float, float, float, float, float):
    """Computes the classic orbital elements directly from values of a two-line element set. These values are less
    accurate instantaneously, however they are much more stable than the oscillating values computed from the SGP4
    module.

    Return items are raan, inc, aop, ecc, sma, meanAnomaly."""

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


class Elements:
    """An object containing the classical orbital elements of a satellite at a given time. All values of the constructor
    are in radians, for degrees use the fromDegrees() class method. A true anomaly can be set if known or will be
    computed otherwise. The epoch parameter is used to compute future or past anomalies."""

    __slots__ = '_raan', '_inc', '_aop', '_ecc', '_sma', '_meanAnomaly', '_trueAnomaly', '_epoch'

    def __init__(self, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                 semiMajorAxis: float, meanAnomaly: float, epoch: 'JulianDate', trueAnomaly: float = None):

        if trueAnomaly is not None:
            self._trueAnomaly = trueAnomaly % TWOPI
        else:
            self._trueAnomaly = meanToTrueAnomaly(meanAnomaly, eccentricity)

        self._raan = raan % TWOPI
        self._inc = inclination
        self._aop = argumentOfPeriapsis % TWOPI
        self._ecc = eccentricity
        self._sma = semiMajorAxis
        self._meanAnomaly = meanAnomaly % TWOPI
        self._epoch = epoch

    @classmethod
    def fromDegrees(cls, raan: float, inclination: float, argumentOfPeriapsis: float, eccentricity: float,
                    semiMajorAxis: float, meanAnomaly: float, epoch: 'JulianDate', trueAnomaly: float = None):
        """Create an Elements object with parameter units of degrees."""

        return cls(radians(raan), radians(inclination), radians(argumentOfPeriapsis), eccentricity,
                   semiMajorAxis, radians(meanAnomaly), epoch, radians(trueAnomaly))

    @classmethod
    def fromTle(cls, tle: 'TwoLineElement', epoch: 'JulianDate'):
        """Create an Elements object from the components of a TLE. The epoch parameter is used to adjust the mean
        elements from the TLE for more accurate values. These elements are less accurate instantaneously, but they are
        much more stable than the oscillating elements computed from the SGP4 algorithm."""

        elements = elementsFromTle(tle, epoch)
        trueAnomaly = meanToTrueAnomaly(elements[-1], elements[3])

        return cls(*elements, epoch, trueAnomaly)

    @classmethod
    def fromState(cls, position: Vector, velocity: Vector, epoch: 'JulianDate', MU: float = EARTH_MU):
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

    @property
    def inc(self):
        """Returns inclination in radians."""
        return self._inc

    @property
    def aop(self):
        """Returns argument of periapsis in radians."""
        return self._aop

    @property
    def ecc(self):
        """Returns eccentricity of the orbit."""
        return self._ecc

    @property
    def sma(self):
        """Returns the semi-major axis in kilometers."""
        return self._sma

    @property
    def meanAnomaly(self):
        """Returns the mean anomaly in radians."""
        return self._meanAnomaly

    @property
    def epoch(self):
        """Returns the epoch as a JulianDate."""
        return self._epoch

    @property
    def trueAnomaly(self):
        """Returns the true anomaly in radians."""
        return self._trueAnomaly

    def updateAnomaly(self, anomaly: str, angle: float, epoch: 'JulianDate'):
        anomalyArg = anomaly.lower()

        if anomalyArg == 'true':
            self._trueAnomaly = angle
            self._meanAnomaly = trueToMeanAnomaly(angle, self._ecc)
        elif anomalyArg == 'mean':
            self._trueAnomaly = meanToTrueAnomaly(angle, self._ecc)
            self._meanAnomaly = angle
        else:
            raise ValueError(f"anomaly must be 'true' or 'mean', not {anomaly}")

        self._epoch = epoch


def radiusAtPeriapsis(semiMajorAxis: float, eccentricity: float) -> float:
    """Computes periapsis from semi-major axis in kilometers and eccentricity."""

    return semiMajorAxis * (1 - eccentricity)


def radiusAtApoapsis(semiMajorAxis: float, eccentricity: float) -> float:
    """Computes apoapsis from semi-major axis in kilometers and eccentricity."""

    return semiMajorAxis * (1 + eccentricity)


def trueToMeanAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly in radians to an eccentric anomaly in radians."""

    eccentricAnomaly = trueToEccentricAnomaly(trueAnomaly, eccentricity)
    return eccentricToMeanAnomaly(eccentricAnomaly, eccentricity)


def trueToEccentricAnomaly(trueAnomaly: float, eccentricity: float) -> float:
    """Converts a true anomaly to an eccentric anomaly in radians."""

    y = sqrt(1 - (eccentricity * eccentricity)) * sin(trueAnomaly)
    return atan3(y, cos(trueAnomaly) + eccentricity)


# todo: come up with an analytical rational for this number
MAXIMUM_ECCENTRICITY = 0.17


def meanToTrueAnomaly(meanAnomaly: float, eccentricity: float) -> float:
    """Computes a mean anomaly in radians to a true anomaly in radians."""

    if eccentricity < MAXIMUM_ECCENTRICITY:
        return _meanToTrueAnomalyFast(meanAnomaly, eccentricity)

    return _meanToTrueAnomalyNewton(meanAnomaly, eccentricity)


def _meanToTrueAnomalyNewton(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to a true anomaly in radians via Newton's method to solve Kepler's
    equation."""

    eccAnom = meanToEccentricAnomaly(meanAnomaly, eccentricity)
    return eccentricToTrueAnomaly(eccAnom, eccentricity)


def _meanToTrueAnomalyFast(meanAnomaly: float, eccentricity: float) -> float:
    """Converts a mean anomaly in radians to a true anomaly in radians via an approximation. The function is valid for
    small eccentricity as error is on the order of (ecc^3)."""

    term1 = 2 * eccentricity * sin(meanAnomaly)
    term2 = 1.25 * eccentricity * eccentricity * sin(2 * meanAnomaly)

    return meanAnomaly + term1 + term2


def meanToEccentricAnomaly(meanAnomaly: float, eccentricity: float) -> float:
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


def eccentricToTrueAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly in radians to a true anomaly in radians."""

    beta = eccentricity / (1 + sqrt(1 - eccentricity * eccentricity))
    return eccentricAnomaly + 2 * atan2(beta * sin(eccentricAnomaly), 1 - beta * cos(eccentricAnomaly))


def eccentricToMeanAnomaly(eccentricAnomaly: float, eccentricity: float) -> float:
    """Converts an eccentric anomaly to a mean anomaly in radians."""

    return eccentricAnomaly - eccentricity * sin(eccentricAnomaly)


def radiusAtAnomaly(sma: float, eccentricity: float, trueAnomaly: float) -> float:
    """Computes the radius of a satellite at a given position in its orbit. Semi-major axis is in kilometers and
    trueAnomaly in radians."""

    numerator = sma * (1 - eccentricity * eccentricity)
    denominator = 1 + eccentricity * cos(trueAnomaly)

    return numerator / denominator


def flightAngleAtAnomaly(eccentricity: float, trueAnomaly: float) -> float:
    """Computes the flight angle between velocity and position vectors at a given true anomaly in radians."""

    numerator = eccentricity * sin(trueAnomaly)
    denominator = 1 + eccentricity * cos(trueAnomaly)

    return atan2(numerator, denominator)


def velocityAtAnomaly(sma: float, radius: float, mu: float) -> float:
    """Compute the magnitude of the velocity at a certain position, where the position is specified by its radius
    in kilometers. Semi-major axis is also in kilometers."""

    tmp = (2 / radius) - (1 / sma)
    return sqrt(mu * tmp)


def nextMeanAnomaly(meanMotion: float, m0: float, epoch0: 'JulianDate', m1: float, time: 'JulianDate') -> 'JulianDate':
    """Finds the time a satellite next achieves a mean anomaly based on the time of a previous anomaly and mean motion.
    Anomalies are in radians and mean motion is in revolutions / day."""

    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion

    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m1 < m0:
        dm += TWOPI

    return time.future((dm / meanMotion) / TWOPI)


def previousMeanAnomaly(meanMotion: float, m0: float, epoch0: 'JulianDate', m1: float, time: 'JulianDate')\
        -> 'JulianDate':
    """Finds the time a satellite previously achieved a mean anomaly based on the time of a previous anomaly and mean
    motion. Anomalies are in radians and mean motion is in revolutions / day."""

    periapsisPass = epoch0.future(-m0 / (TWOPI * meanMotion))
    revolutions = (time - periapsisPass) * meanMotion

    m0 = (revolutions - floor(revolutions)) * TWOPI
    dm = m1 - m0
    if m0 < m1:
        dm -= TWOPI

    return time.future((dm / meanMotion) / TWOPI)


def nextTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: 'JulianDate',
                    t1: float, time: 'JulianDate') -> 'JulianDate':
    """Finds the time a satellite next achieves a true anomaly based on the time of a previous anomaly and mean motion.
    Anomalies are in radians and mean motion is in revolutions / day."""

    m0 = trueToMeanAnomaly(t0, eccentricity)
    m1 = trueToMeanAnomaly(t1, eccentricity)

    return nextMeanAnomaly(meanMotion, m0, epoch0, m1, time)


def previousTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0: 'JulianDate',
                        t1: float, time: 'JulianDate') -> 'JulianDate':
    """Finds the time a satellite previously achieved a true anomaly based on the time of a previous anomaly and mean
    motion. Anomalies are in radians and mean motion is in revolutions / day."""

    m0 = trueToMeanAnomaly(t0, eccentricity)
    m1 = trueToMeanAnomaly(t1, eccentricity)

    return previousMeanAnomaly(meanMotion, m0, epoch0, m1, time)


def nearestTrueAnomaly(meanMotion: float, eccentricity: float, t0: float, epoch0, t1: float):
    """Finds the time of the nearest true anomaly (radians), either forward or backward in time. Mean motion
    must be in revolutions / day."""

    m0 = trueToMeanAnomaly(t0, eccentricity)
    m1 = trueToMeanAnomaly(t1, eccentricity)

    return nearestMeanAnomaly(meanMotion, m0, epoch0, m1)


def nearestMeanAnomaly(meanMotion: float, m0: float, epoch0, m1: float):
    """Finds the time of the nearest mean anomaly (radians), either forward or backward in time. Mean motion
    must be in revolutions / day."""

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

    return epoch0.future(dma / (meanMotion * TWOPI))


def meanAnomalyAtTime(meanMotion: float, m0: float, epoch0: 'JulianDate', time: 'JulianDate') -> float:
    """Computes the mean anomaly at a given time. Mean anomaly is in radians and mean motion is in revolutions / day."""

    dM = meanMotion * (time - epoch0) * SECONDS_PER_DAY
    return (m0 + dM) % TWOPI


def trueAnomalyAtTime(meanMotion: float, eccentricity: float, t0: float, epoch0: 'JulianDate', time: 'JulianDate')\
        -> float:
    """Computes the true anomaly at a given time. True anomaly is in radians and mean motion is in revolutions / day."""

    m0 = trueToMeanAnomaly(t0, eccentricity)
    m1 = meanAnomalyAtTime(meanMotion, m0, epoch0, time)

    return meanToTrueAnomaly(m1, eccentricity)


def computeAnomaly(anomaly: str, position: Vector, velocity: Vector, mu: float = EARTH_MU) -> float:
    """Computes the true anomaly in radians from a set of state vectors. Mu can be any Body's
    standard gravitational parameter, but defaults to Earth's mu."""

    elements = elementsFromState(position, velocity, mu)
    anomalyArg = anomaly.lower()

    if anomalyArg == 'true':
        return elements[-1]
    elif anomalyArg == 'mean':
        return elements[-2]
    else:
        raise ValueError(f"anomaly must be 'true' or 'mean', not {anomaly}")


def smaToMeanMotion(semiMajorAxis: float, mu: float) -> float:
    """Converts a semi-major axis in kilometers to a mean motion in radians / second."""

    smaCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis
    return sqrt(mu / smaCubed)


def meanMotionToSma(meanMotion: float, mu: float) -> float:
    """Converts a mean motion in radians / second to semi-major axis in kilometers."""

    oneThird = 1 / 3
    return (mu ** oneThird) / (meanMotion ** (oneThird * 2))
