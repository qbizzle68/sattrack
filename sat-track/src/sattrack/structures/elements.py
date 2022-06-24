from pyevspace import EVector, dot, cross, norm, vang
from math import sqrt, radians, degrees, pi, sin, cos, acos, atan

from sattrack.rotation.order import Order
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, rotateMatrixFrom
from sattrack.spacetime.juliandate import JulianDate
from sattrack.util.anomalies import meanToTrue, trueToMean, trueToEccentric
from sattrack.util.constants import EARTH_MU
from sattrack.util.conversions import smaToMeanMotion, meanMotionToSma
from sattrack.structures.tle import TwoLineElement


class OrbitalElements:
    """This object contains a celestial object's classical orbital elements, namely:
        semi-major axis, eccentricity, inclination, right-ascension of ascending node, argument of periapsis,
        true anomaly and the epoch of the true anomaly."""

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
        self._inc = inc
        self._raan = raan
        self._aop = aop
        self._meanAnomaly = meanAnomaly
        self._epoch = epoch

    @classmethod
    def fromTLE(cls, tle: TwoLineElement, jd: JulianDate):
        """Class method used to instantiate an object from a two-line element object at a given time.
        Parameters:
        tle:    Two-Line element of the object to obtain orbital elements of.
        jd:     Julian Date of the time to compute orbital elements of."""

        dt = jd.difference(tle.epoch())
        inc = tle.inclination()
        nConvert = tle.meanMotion() * 2 * pi / 86400
        a0 = ((EARTH_MU / (nConvert * nConvert)) ** (1 / 3)) / 1000  # km
        dM = (tle.meanMotion() * dt
              + tle.meanMotionDot() * dt * dt
              + tle.meanMotionDDot() * dt * dt * dt) * 2 * pi
        M1 = (2 * pi * (tle.meanAnomaly() / 360) + dM % (2 * pi))
        tleEcc = tle.eccentricity()
        n0 = tle.meanMotion()
        n0dot = tle.meanMotionDot() * 2
        aDot = -2 * a0 * n0dot / (3 * n0)
        sma = (a0 + aDot * dt) * 1000
        eDot = -2 * (1 - tleEcc) * n0dot / (3 * n0)
        ecc = tleEcc + eDot * dt # todo: do we use this for the next line?
        temp = (a0 ** -3.5) / ((1 - (tleEcc * tleEcc)) ** 2)
        # todo: include the moon and sun perturbations to this
        lanJ2Dot = -2.06474e14 * temp * cos(radians(inc))
        ra = tle.raan() + lanJ2Dot * dt
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * ((sin(radians(inc))) ** 2))
        aop = tle.argumentOfPeriapsis() + aopJ2Dot * dt
        return cls(sma=sma, ecc=ecc, inc=inc, raan=ra, aop=aop, meanAnomaly=M1, epoch=jd)

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
        inc = degrees(acos(angMom[2] / angMom.mag()))
        ra = degrees(acos(lineOfNodes[0] / lineOfNodes.mag()))
        if lineOfNodes[1] < 0:
            ra = 360 - ra   # todo: AOP is significantly off, why?
        aop = degrees(acos(dot(lineOfNodes, eccVec) / lineOfNodes.mag() * eccVec.mag()))
        if eccVec[2] < 0:
            aop = 360 - aop
        tAnom = degrees(acos(dot(eccVec, position) / (eccVec.mag() * position.mag())))
        if dot(position, velocity) < 0:
            tAnom = 360 - tAnom
        sma = (angMom.mag() ** 2) / ((1 - (ecc * ecc)) * EARTH_MU)
        return cls(sma=sma, ecc=ecc, inc=inc, raan=ra, aop=aop, meanAnomaly=trueToMean(tAnom, ecc), epoch=jd)

    def __str__(self) -> str:
        """Returns a string representation of the orbital elements."""

        return (f'Semi-major axis: {self._sma}, Eccentricity: {self._ecc}, Inclination: {self._inc}\n'
                + f'RAAN: {self._raan}, Argument of Perigee: {self._aop}, Mean Anomaly: {self._meanAnomaly}\n'
                + f'Epoch: {str(self._epoch)}')

    def timeToMeanAnomaly(self, meanAnomaly: float) -> float:
        """Computes the time taken for the object to reach a given mean anomaly. To incorporate multiple orbits,
        increase the argument by the orbit count multiple of 360 degrees.
        Parameters:
        meanAnomaly:    The mean anomaly to compute the time until, measured in degrees."""

        n = smaToMeanMotion(self._sma)
        dM = meanAnomaly - self._meanAnomaly
        if dM < 0:
            dM %= 360
        return radians(dM) / n

    def meanAnomalyAt(self, jd: JulianDate) -> float:
        """Computes the mean anomaly at a given time.
        Parameters:
        jd: The time to compute the mean anomaly.
        Returns the mean anomaly in degrees."""

        if self._epoch == 0:
            raise ValueError("Epoch was not set for this instance.")
        n = smaToMeanMotion(self._sma)
        dt = jd.difference(self._epoch)
        mRad = n * dt + radians(self._meanAnomaly)
        return degrees(mRad) % 360

    def getState(self, jd: JulianDate = None) -> tuple[EVector]:
        if not jd:
            jd = self._epoch
        tAnom = meanToTrue(self.meanAnomalyAt(jd), self._ecc)
        eAnom = trueToEccentric(tAnom, self._ecc)
        r = radiusAtTrueAnomaly(self._sma, self._ecc, tAnom)
        pOrbit = EVector(cos(radians(tAnom)), sin(radians(tAnom)), 0) * r
        vOrbit = EVector(-sin(radians(eAnom)), sqrt(1 - self._ecc * self._ecc) * cos(radians(eAnom)), 0) * (
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
        """Sets the inclination of the orbit in degrees."""
        self._inc = inc

    def getInc(self):
        """Returns the inclination of the orbit in degrees."""
        return self._inc

    def setRaan(self, raan: float):
        """Sets the right-ascension of the ascending node in degrees.
            Synonymous with longitude of ascending node."""
        self._raan = raan

    def getRaan(self):
        """Returns the right-ascension of the ascending node in degrees.
            Synonymous with longitude of ascending node."""
        return self._raan

    def setAop(self, aop: float):
        """Sets the argument of periapsis measured in degrees."""
        self._aop = aop

    def getAop(self):
        """Returns the argument of periapsis measured in degrees."""
        return self._aop

    def setMeanAnomaly(self, meanAnom: float):
        """Sets the mean anomaly of the orbit in degrees. For the positional methods to work properly
        the epoch should be changed whenever this value is changed."""
        self._meanAnomaly = meanAnom

    def getMeanAnomaly(self):
        """Returns the mean anomaly of the orbit in degrees."""
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
    returns: The rate of procession of the RAAN in degrees per day."""

    a = meanMotionToSma(tle.meanMotion()) / 1000.0
    dRaanSphere = -2.06474e14 * (a ** -3.5) * cos(radians(tle.inclination())) / ((1 - tle.eccentricity() * tle.eccentricity()) ** 2)
    dRaanMoon = -0.00338 * cos(radians(tle.inclination())) / tle.meanMotion()
    dRaanSun = -0.00154 * cos(radians(tle.inclination())) / tle.meanMotion()
    return dRaanSphere + dRaanMoon + dRaanSun


def aopProcession(tle: TwoLineElement) -> float:
    """Computes the procession of the argument of periapsis due to third-body perturbations and the
    non-spherical Earth. Values are positive for a prograde (direct) orbit, negative for retrograde.
    Parameters:
    tle:    Two-line element set for an object.
    returns: The rate of procession of the AOP in degrees per day."""

    a = meanMotionToSma(tle.meanMotion()) / 1000.0
    dAopSphere = 1.03237e14 * (a ** -3.5) * (4 - 5 * (sin(radians(tle.inclination())) ** 2)) / ((1 - tle.eccentricity() * tle.eccentricity()) ** 2)
    dAopMoon = 0.00169 * (4 - 5 * (sin(radians(tle.inclination())) ** 2)) / tle.meanMotion()
    dAopSun = 0.00077 * (4 - 5 * (sin(radians(tle.inclination())) ** 2)) / tle.meanMotion()
    return dAopSphere + dAopMoon + dAopSun


def computeEccentricVector(position: EVector, velocity: EVector) -> EVector:
    """Compute the eccentric vector of an orbit.
    Parameters:
    position:   Position state of the object.
    velocity:   Velocity state of the object."""
    rtn = velocity * dot(position, velocity)
    rtn = position * ((velocity.mag() ** 2) - (EARTH_MU / position.mag())) - rtn
    return rtn / EARTH_MU


def radiusFromElements(elements: OrbitalElements) -> float:
    return radiusAtTrueAnomaly(elements.getSma(), elements.getEcc(),
                               meanToTrue(elements.getMeanAnomaly(), elements.getEcc()))


def radiusAtTrueAnomaly(sma: float, ecc: float, trueAnomaly: float) -> float:
    return sma * (1 - ecc * ecc) / (1 + ecc * cos(radians(trueAnomaly)))


def flightAngleFromElements(elements: OrbitalElements) -> float:
    return flightAngleAtTrueAnomaly(elements.getEcc(), meanToTrue(elements.getMeanAnomaly(), elements.getEcc()))


def flightAngleAtTrueAnomaly(ecc: float, trueAnom: float) -> float:
    taRad = radians(trueAnom)
    return atan((ecc * sin(taRad)) / (1 + ecc * cos(taRad)))


def velocityFromElements(elements: OrbitalElements) -> float:
    return velocityAtTrueAnomaly(elements.getSma(), elements.getEcc(),
                                 meanToTrue(elements.getMeanAnomaly(), elements.getEcc()))


def velocityAtTrueAnomaly(sma: float, ecc: float, trueAnom: float) -> float:
    r = radiusAtTrueAnomaly(sma, ecc, trueAnom)
    return sqrt(EARTH_MU * ((2 / r) - (1 / sma)))


'''how can we compute this without a velocity vector'''
'''idea: use another vector in the orbital plane (second position or line of nodes?) to generate angular momentum,
    use this vector to create a quaternion or a rotation matrix based on the angle-axis rotation, and rotate the 
    position vector by (90 - flight angle) ... need true anomaly for flight angle'''
'''new idea: generate an angular momentum direction with some other vector, generate line of nodes with this vector,
    find the angle from this vector, then account for aop to find true anomaly. need to figure out how to differentiate
    from angles between 0-180 and 180-360.'''


# todo: why does this return NaN for 0 and 180
def trueAnomalyFromState(position: EVector, velocity: EVector) -> float:
    eccVec = computeEccentricVector(position, velocity)
    ang = vang(position, eccVec)
    return 360 - ang if norm(cross(position, eccVec)) == norm(cross(position, velocity)) else ang


def trueAnomalyFromMomentum(position: EVector, momentum: EVector, aop: float) -> float:
    lineOfNodes = cross(EVector(0, 0, 1), momentum)
    ang = vang(position, lineOfNodes)
    angAdjust = 360 - ang if norm(cross(position, lineOfNodes)) == norm(momentum) else ang
    return (angAdjust - aop) % 360
