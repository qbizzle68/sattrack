from pyevspace import EVector, dot, cross, norm, vang

from rotation.order import Order
from rotation.rotation import getEulerMatrix, EulerAngles, rotateMatrixFrom
from spacetime.spacetime import JulianDate
from util.constants import EARTH_MU, CJ2
from util.conversions import atan2
from math import sqrt, radians, degrees, pi, sin, cos, acos, atan

from structures.tle import TwoLineElement


class OrbitalElements:
    """This object contains a celestial object's classical orbital elements, namely:
        semi-major axis, eccentricity, inclination, right-ascension of ascending node, argument of periapsis,
        true anomaly and the epoch of the true anomaly.
    Notice the true anomaly is used to represent the angular position of the orbit as opposed to the mean anomaly.
    The epoch is not required in instantiate an object, but the class supports methods to move the object in time,
    and behavior of the methods is undefined if the epoch is not set to the matching true anomaly."""

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
        ecc = tleEcc + eDot * dt
        temp = (a0 ** -3.5) * ((1 - (tleEcc * tleEcc)) ** -2)
        lanJ2Dot = CJ2 * temp * cos(radians(inc))
        ra = tle.raan() + lanJ2Dot * dt
        aopJ2Dot = (CJ2 / 2) * temp * (4 - 5 * ((sin(radians(inc))) ** 2))
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
        eccVec = cls.computeEccentricVector(position, velocity)

        ecc = eccVec.mag()
        inc = degrees(acos(angMom[2] / angMom.mag()))
        ra = degrees(acos(lineOfNodes[0] / lineOfNodes.mag()))
        if lineOfNodes[1] < 0:
            ra = 360 - ra
        aop = degrees(acos(dot(lineOfNodes, eccVec) / lineOfNodes.mag() * eccVec.mag()))
        if eccVec[2] < 0:
            aop = 360 - aop
        tAnom = degrees(acos(dot(eccVec, position) / (eccVec.mag() * position.mag())))
        if dot(position, velocity) < 0:
            tAnom = 360 - tAnom
        sma = (angMom.mag() ** 2) / ((1 - (ecc * ecc)) * EARTH_MU)
        return cls(sma=sma, ecc=ecc, inc=inc, raan=ra, aop=aop, meanAnomaly=trueToMean(tAnom, ecc), epoch=jd)

    @classmethod
    def computeEccentricVector(cls, position: EVector, velocity: EVector):
        """Compute the eccentric vector of an orbit.
        Parameters:
        position:   Position state of the object.
        velocity:   Velocity state of the object."""
        rtn = velocity * dot(position, velocity)
        rtn = position * ((velocity.mag() ** 2) - (EARTH_MU / position.mag())) - rtn
        return rtn / EARTH_MU

    def __str__(self) -> str:
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

    '''def moveTo(self, jd: JulianDate):
        """Moves a satellites angular position to a specified time using the mean motion of the object. Time is
        determined by a future or past Julian Date.s
        Parameters:
        jd: The Julian Date you wish to move the object to."""
        self.moveBy(jd.difference(self._epoch))

    def moveBy(self, dt: float):
        """Moves a satellites angular position by a number of solar days, using the mean motion of the object.
        A negative value indicates a past position.
        Parameters:
        dt: The number of solar days to move the object to."""
        dtp = 86400 * dt
        M1rad = radians(trueToMean(self._trueAnomaly, self._ecc)) + dtp * sqrt(EARTH_MU / (self._sma ** 3))
        self._trueAnomaly = meanToTrue(degrees(M1rad), self._ecc)
        self._epoch = self._epoch.future(dt)

    def timeUntil(self, trueAnomaly: float) -> float:
        """Returns the time until the object achieves a specified true anomaly.
            If the trueAnomaly parameter is equal to the current true anomaly, the returned value will be 0,
            NOT the orbital period. I.e. the angle needed to travel will be 0, not 360 degrees.
        Parameters:
        trueAnomaly:    The time computed will be the time it takes to travel to this true anomaly."""
        M0 = trueToMean(self._trueAnomaly, self._ecc)
        M1 = trueToMean(trueAnomaly, self._ecc)
        da = M1 - M0 if M1 >= M0 else M1 + 360.0 - M0
        return radians(da) / sqrt(EARTH_MU / (self._sma ** 3))'''

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


def computeEccentricVector(position: EVector, velocity: EVector) -> EVector:
    """Compute the eccentric vector of an orbit.
    Parameters:
    position:   Position state of the object.
    velocity:   Velocity state of the object."""
    rtn = velocity * dot(position, velocity)
    rtn = position * ((velocity.mag() ** 2) - (EARTH_MU / position.mag())) - rtn
    return rtn / EARTH_MU


def meanMotionToSma(meanMotion: float) -> float:
    """Converts mean motion from a two-line element to semi-major axis using
    Kepler's 2nd law.
    Parameters:
    meanMotion: Mean motion measured in revolutions per day.
    Returns the semi-major axis in meters."""
    mMotionRad = meanMotion * 2 * pi / 86400.0
    return (EARTH_MU ** (1.0 / 3.0)) / (mMotionRad ** (2.0 / 3.0))


def smaToMeanMotion(sma: float) -> float:
    """Converts semi-major axis to mean motion.
    Parameters:
    sma: Semi-major axis measured in meters.
    Returns the mean motion in radians per second."""
    return sqrt(EARTH_MU / (sma ** 3))


def meanToTrue(meanAnomaly: float, eccentricity: float) -> float:
    return eccentricToTrue(
        meanToEccentric(meanAnomaly, eccentricity),
        eccentricity
    )


def meanToEccentric(meanAnomaly: float, eccentricity: float) -> float:
    e = __m2ENewtonRaphson(meanAnomaly, meanAnomaly, eccentricity) % 360.0
    return e if e >= 0 else e + 360.0


def eccentricToTrue(eccAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(eccAnom))
    return degrees(atan2(y, cos(radians(eccAnom)) - ecc))


def trueToMean(trueAnom: float, ecc: float) -> float:
    return eccentricToMean(
        trueToEccentric(trueAnom, ecc), ecc
    )


def trueToEccentric(trueAnom: float, ecc: float) -> float:
    y = sqrt(1 - (ecc * ecc)) * sin(radians(trueAnom))
    return degrees(atan2(y, cos(radians(trueAnom)) + ecc))


def eccentricToMean(eccAnom: float, ecc: float) -> float:
    m = (degrees(radians(eccAnom) - ecc * sin(radians(eccAnom)))) % 360.0
    return m if m >= 0 else m + 360.0


def __m2ENewtonRaphson(M: float, Ej: float, ecc: float) -> float:
    M = radians(M)
    Ej = radians(Ej)
    while True:
        num = Ej - (ecc * sin(Ej)) - M
        den = 1 - ecc * cos(Ej)
        Ej1 = Ej - (num / den)
        if abs(Ej1 - Ej) <= 1e-7:
            return degrees(Ej1)
        else:
            Ej = Ej1
    '''Ej1 = Ej - ((Ej - ecc * sin(radians(Ej)) - M) / (1 - ecc * cos(radians(Ej))))
    return Ej1 if abs(Ej1 - Ej) <= 1e-7 else __m2ENewtonRaphson(M, Ej1, ecc)'''

def __m2E(M, Ej, ecc):
    return Ej - ((Ej - ecc * sin(radians(Ej)) - M) / (1 - ecc * cos(radians(Ej))))

def __m2ENewtonRaphson2(M: float, ecc: float) -> float:
    Ei1 = M
    for i in range(0):
        Ei = Ei1
        Eo = __m2E(M, Ei, ecc)
        Eoo = __m2E(M, Eo, ecc)
        Ei1 = 0.5 * (Eo + Eoo)
    return __m2ENewtonRaphson(M, Ei1, ecc)


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
    use this vector to create a quaternion or a rotation matrix based on the angle-axis rotation, and rotate the position
    vector by (90 - flight angle) ... need true anomaly for flight angle'''
'''new idea: generate an angular momentum direction with some other vector, generate line of nodes with this vector,
    find the angle from this vector, then account for aop to find true anomaly. need to figure out how to differentiate
    from angles between 0-180 and 180-360.'''
#todo: why does this return NaN for 0 and 180
def trueAnomalyFromState(position: EVector, velocity: EVector) -> float:
    eccVec = OrbitalElements.computeEccentricVector(position, velocity)
    ang = vang(position, eccVec)
    return 360 - ang if norm(cross(position, eccVec)) == norm(cross(position, velocity)) else ang

def trueAnomalyFromMomentum(position: EVector, momentum: EVector, aop: float) -> float:
    lineOfNodes = cross(EVector(0, 0, 1), momentum)
    ang = vang(position, lineOfNodes)
    angAdjust = 360 - ang if norm(cross(position, lineOfNodes)) == norm(momentum) else ang
    return (angAdjust - aop) % 360
