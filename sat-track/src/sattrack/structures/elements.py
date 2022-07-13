from math import sqrt, radians, degrees, pi, sin, cos, acos, atan

from pyevspace import EVector, dot, cross, norm

from sattrack.rotation.order import Order
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateMatrixFrom
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import meanToTrue, trueToMean, trueToEccentric
from sattrack.util.constants import EARTH_MU
from sattrack.util.conversions import meanMotionToSma, smaToMeanMotion


class OrbitalElements:
    """This object contains a celestial object's classical orbital elements, namely:
        semi-major axis, eccentricity, inclination, right-ascension of ascending node, argument of periapsis,
        true anomaly and the epoch of the true anomaly."""

    # todo: set epoch to default to None
    def __init__(self, *, sma: float, ecc: float, inc: float, raan: float, aop: float, meanAnomaly: float,
                 epoch: JulianDate = 0):
        """Constructs an instance with the given values, with the epoch optional.
        Parameters:
        sma:    Semi-major axis in meters.
        ecc:    Eccentricity of the orbit.
        inc:    Inclination of the orbit in degrees.
        raan:   Right-ascension of the ascending node, measured in degrees.
        aop:    Argument of Periapsis, measured in degrees.
        meanAnomaly:    Mean anomaly of the orbit at epoch.
        epoch:  Epoch used to distinguish when the mean anomaly occurs."""
        self._sma = sma
        self._ecc = ecc
        self._inc = radians(inc)
        self._raan = radians(raan)
        self._aop = radians(aop)
        self._meanAnomaly = radians(meanAnomaly)
        self._epoch = epoch

    @classmethod
    def fromTle(cls, tle: TwoLineElement, jd: JulianDate = 0):
        """Class method used to instantiate an object from a two-line element object at a given time.
        Parameters:
            tle -- Two-Line element of the object to obtain orbital elements of.
            jd -- Julian Date of the time to compute orbital elements of"""

        if jd == 0:
            jd = tle.epoch()
        dt = jd.difference(tle.epoch())
        inc = radians(tle.inclination())
        n = tle.meanMotionRad()
        a0 = ((EARTH_MU / (n*n)) ** (1/3))  # km todo: replace this with the sma() method?
        dM = (tle.meanMotion() * dt # todo: regroup this for better efficiency
              + tle.meanMotionDot() * dt * dt
              + tle.meanMotionDDot() * dt * dt * dt) * 2 * pi
        M1 = (radians(tle.meanAnomaly()) + dM) % (2 * pi)
        tleEcc = tle.eccentricity()
        n0 = tle.meanMotion()
        n0dot = tle.meanMotionDot() * 2
        aDot = -2 * a0 * n0dot / (3 * n0)
        sma = a0 + aDot * dt
        eDot = -2 * (1 - tleEcc) * n0dot / (3 * n0)
        ecc = tleEcc + eDot * dt
        temp = (a0 ** -3.5) / ((1 - (tleEcc * tleEcc)) ** 2)
        #   Perturbations
        #   Non-spherical Earth
        lanJ2Dot = -2.06474e14 * temp * cos(inc)
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * ((sin(inc)) ** 2))
        #   Third-Body
        lanMoon = -0.00338 * cos(inc) / n0
        lanSun = -0.00154 * cos(inc) / n0
        aopMoon = 0.00169 * (4 - (5 * sin(inc) ** 2)) / n0
        aopSun = 0.00077 * (4 - (5 * sin(inc) **2)) / n0
        ra = (tle.raan() + (lanJ2Dot + lanMoon + lanSun) * dt) % 360
        aop = (tle.argumentOfPeriapsis() + (aopJ2Dot + aopMoon + aopSun) * dt) % 360
        return cls(sma=sma, ecc=ecc, inc=degrees(inc), raan=ra, aop=aop, meanAnomaly=degrees(M1), epoch=jd)

        '''if jd == 0:
            jd = tle.epoch()
        dt = jd.difference(tle.epoch())
        inc = tle.inclination()
        nConvert = tle.meanMotion() * 2 * pi / 86400
        a0 = ((EARTH_MU / (nConvert * nConvert)) ** (1 / 3))  # km
        dM = (tle.meanMotion() * dt
              + tle.meanMotionDot() * dt * dt
              + tle.meanMotionDDot() * dt * dt * dt) * 2 * pi
        M1 = ((pi * (tle.meanAnomaly() / 180) + dM) % (2 * pi))
        tleEcc = tle.eccentricity()
        n0 = tle.meanMotion()
        n0dot = tle.meanMotionDot() * 2
        aDot = -2 * a0 * n0dot / (3 * n0)
        sma = (a0 + aDot * dt)
        eDot = -2 * (1 - tleEcc) * n0dot / (3 * n0)
        ecc = tleEcc + eDot * dt
        temp = (a0 ** -3.5) / ((1 - (tleEcc * tleEcc)) ** 2)
        #   Perturbations
        #   Non-spherical Earth
        lanJ2Dot = -2.06474e14 * temp * cos(radians(inc))
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * ((sin(radians(inc))) ** 2))
        #   Third-Body
        lanMoon = -0.00338 * cos(inc) / n0
        lanSun = -0.00154 * cos(inc) / n0
        aopMoon = 0.00169 * (4 - (5 * sin(inc) ** 2)) / n0
        aopSun = 0.00077 * (4 - (5 * sin(inc) ** 2)) / n0
        ra = tle.raan() + (lanJ2Dot + lanMoon + lanSun) * dt
        aop = tle.argumentOfPeriapsis() + (aopJ2Dot + aopMoon + aopSun) * dt
        return cls(sma=sma, ecc=ecc, inc=inc, raan=ra, aop=aop, meanAnomaly=degrees(M1), epoch=jd)'''

    @classmethod
    def fromState(cls, position: EVector, velocity: EVector, jd: JulianDate = 0):
        """Class method used to instantiate an object from a known state at a given time.
        Parameters:
        position:   Position state of the object at epoch.
        velocity:   Velocity state of the object at epoch.
        jd:         Julian Date of the time to compute orbital elements of."""

        angMom = cross(position, velocity)
        lineOfNodes = norm(cross(EVector(0, 0, 1), angMom))
        eccVec = computeEccentricVector(position, velocity)

        ecc = eccVec.mag()
        inc = acos(angMom[2] / angMom.mag())
        ra = acos(lineOfNodes[0] / lineOfNodes.mag())
        if lineOfNodes[1] < 0:
            ra = (2*pi) - ra
        aop = acos(dot(lineOfNodes, eccVec) / (lineOfNodes.mag() * eccVec.mag()))
        if eccVec[2] < 0:
            aop = (2*pi) - aop
        tAnom = acos(dot(eccVec, position) / (eccVec.mag() * position.mag()))
        if dot(position, velocity) < 0:
            tAnom = (2*pi) - tAnom
        sma = (angMom.mag() ** 2) / ((1 - (ecc * ecc)) * EARTH_MU)
        return cls(sma=sma, ecc=ecc, inc=inc, raan=ra, aop=aop, meanAnomaly=trueToMean(tAnom, ecc), epoch=jd)

    def __str__(self) -> str:
        """Returns a string representation of the orbital elements."""

        return (f'Semi-major axis: {self._sma}, Eccentricity: {self._ecc}, Inclination: {degrees(self._inc)}\n'
                + f'RAAN: {degrees(self._raan)}, Argument of Perigee: {degrees(self._aop)}, '
                  f'Mean Anomaly: {degrees(self._meanAnomaly)}\nEpoch: {str(self._epoch)}')

    def getState(self, jd: JulianDate = None) -> tuple[EVector]:
        """Computes state vectors based on the orbital elements at a given time.
        Parameters:
        jd:     Time to compute the state vectors.
        returns: A tuple containing the position and velocity vectors in m and m/s respectively."""
        if jd is None:
            jd = self._epoch
        tAnom = meanToTrue(meanAnomalyAt(self, jd), self._ecc)
        eAnom = trueToEccentric(tAnom, self._ecc)
        r = self._sma * (1 - self._ecc * cos(eAnom))
        pOrbit = EVector(cos(tAnom), sin(tAnom), 0) * r
        vOrbit = EVector(-sin(eAnom), sqrt(1 - self._ecc * self._ecc) * cos(eAnom), 0) * (
                sqrt(EARTH_MU * self._sma) / r)
        rot = getEulerMatrix(Order.ZXZ, EulerAngles(self._raan, self._inc, self._aop))
        return rotateMatrixFrom(rot, pOrbit), rotateMatrixFrom(rot, vOrbit)

    def setSma(self, sma: float):
        """Sets the semi-major axis of the orbit.
        Parameters:
        sma:    Semi-major axis measured in meters."""
        self._sma = sma

    def getSma(self) -> float:
        """Returns the semi-major axis of the orbit in meters."""
        return self._sma

    def setEcc(self, ecc: float):
        """Sets the eccentricity of the orbit."""
        self._ecc = ecc

    def getEcc(self):
        """Returns the eccentricity of the orbit."""
        return self._ecc

    def setInc(self, inc: float):
        """Sets the inclination of the orbit in radians."""
        self._inc = inc

    def getInc(self):
        """Returns the inclination of the orbit in radians."""
        return self._inc

    def setRaan(self, raan: float):
        """Sets the right-ascension of the ascending node in radians.
            Synonymous with longitude of ascending node."""
        self._raan = raan

    def getRaan(self):
        """Returns the right-ascension of the ascending node in radians.
            Synonymous with longitude of ascending node."""
        return self._raan

    def setAop(self, aop: float):
        """Sets the argument of periapsis measured in radians."""
        self._aop = aop

    def getAop(self):
        """Returns the argument of periapsis measured in radians."""
        return self._aop

    def setMeanAnomaly(self, meanAnom: float):
        """Sets the mean anomaly of the orbit in radians. For the positional methods to work properly
        the epoch should be changed whenever this value is changed."""
        self._meanAnomaly = meanAnom

    def getMeanAnomaly(self):
        """Returns the mean anomaly of the orbit in radians."""
        return self._meanAnomaly

    def setEpoch(self, epoch: JulianDate):
        """Sets the epoch associated with the true anomaly. For the positional methods to work properly,
        the true anomaly should be changed whenever this value is changed. Argument should be a JulianDate object."""
        self._epoch = epoch

    def getEpoch(self):
        """Returns the epoch associated with the objects true anomaly."""
        return self._epoch


def raanProcession(tle: TwoLineElement) -> float:
    """Computes the procession of the right-ascension of the ascending node due to third-body perturbations
     and the non-spherical earth. Values are negative for a prograde (direct) orbit, positive for retrograde.
    Parameters:
    tle:    Two-line element set for an object.
    returns: The rate of procession of the RAAN in radians per day."""

    a = meanMotionToSma(tle.meanMotion())
    dRaanSphere = -2.06474e14 * (a ** -3.5) * cos(tle.inclination()) / (
                (1 - tle.eccentricity() * tle.eccentricity()) ** 2)
    dRaanMoon = -0.00338 * cos(tle.inclination()) / tle.meanMotion()
    dRaanSun = -0.00154 * cos(tle.inclination()) / tle.meanMotion()
    return radians(dRaanSphere + dRaanMoon + dRaanSun)


def aopProcession(tle: TwoLineElement) -> float:
    """Computes the procession of the argument of periapsis due to third-body perturbations and the
    non-spherical Earth. Values are positive for a prograde (direct) orbit, negative for retrograde.
    Parameters:
    tle:    Two-line element set for an object.
    returns: The rate of procession of the AOP in radians per day."""

    a = meanMotionToSma(tle.meanMotion())
    dAopSphere = 1.03237e14 * (a ** -3.5) * (4 - 5 * (sin(tle.inclination()) ** 2)) / \
        ((1 - tle.eccentricity() * tle.eccentricity()) ** 2)
    dAopMoon = 0.00169 * (4 - 5 * (sin(tle.inclination()) ** 2)) / tle.meanMotion()
    dAopSun = 0.00077 * (4 - 5 * (sin(tle.inclination()) ** 2)) / tle.meanMotion()
    return radians(dAopSphere + dAopMoon + dAopSun)


# todo: figure out if the magnitude here is accurate
def computeEccentricVector(position: EVector, velocity: EVector) -> EVector:
    """Compute the eccentric vector of an orbit.
    Parameters:
    position:   Position state of the object.
    velocity:   Velocity state of the object."""
    lhs = position * (velocity.mag2() / EARTH_MU - 1 / position.mag())
    rhs = velocity * (dot(position, velocity) / EARTH_MU)
    return lhs - rhs


def computeVelocity(elements: OrbitalElements) -> float:
    radius = computeRadius(elements)
    return sqrt(EARTH_MU * ((2 / radius) - (1 / elements.getSma())))


def computeRadius(elements: OrbitalElements) -> float:
    tAnom = meanToTrue(elements.getMeanAnomaly(), elements.getEcc())
    ecc = elements.getEcc()
    return elements.getSma() * (1 - ecc * ecc) / (1 + ecc * cos(tAnom))


def computeFlightAngle(elements: OrbitalElements) -> float:
    tAnom = meanToTrue(elements.getMeanAnomaly(), elements.getEcc())
    return atan((elements.getEcc() * sin(tAnom)) / (1 + elements.getEcc() * cos(tAnom)))


def meanAnomalyAt(elements: OrbitalElements, jd: JulianDate) -> float:
    if elements.getEpoch() == 0:
        raise ValueError('Epoch was not set for this instance.')
    n = smaToMeanMotion(elements.getSma())
    dt = jd.difference(elements.getEpoch()) * 86400.0
    mAnom = n * dt + elements.getMeanAnomaly()
    return mAnom % (2*pi)


def timeToMeanAnomaly(elements: OrbitalElements, meanAnom: float) -> float:
    n = smaToMeanMotion(elements.getSma())
    m1 = elements.getMeanAnomaly()
    dM = (meanAnom - m1) % (2*pi)
    if meanAnom < m1:
        dM += (2*pi)
    return dM / n / 86400.0
