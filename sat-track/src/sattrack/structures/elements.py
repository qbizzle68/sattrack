from math import sqrt, radians, degrees, pi, sin, cos, acos, atan

from pyevspace import EVector, dot, cross, norm

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, ReferenceFrame
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import meanToTrue, trueToMean, trueToEccentric
from sattrack.util.constants import EARTH_MU, TWOPI
from sattrack.util.conversions import meanMotionToSma, smaToMeanMotion


class OrbitalElements:
    """A container class that computes and contains orbital elements of a satellite, either imported manually, or
    derived from an approximation of TLE values at a given time. The mean anomaly is the anomaly used for positioning,
    and an epoch can be set to compute times for future/past anomalies. The initialized angles are in degrees, for ease
    of use, but the angles are converted and implemented as radians."""

    def __init__(self, *, sma: float, ecc: float, inc: float, raan: float, aop: float, meanAnomaly: float,
                 epoch: JulianDate = None):
        """
        Constructs an instance with the given values. Not including an epoch will not allow future
        anomalies to be computed. For ease of user interface the angle parameters are measured in degrees,
        but are converted into radians to be consistent with the rest of the APIs. Due to the number of
        parameters, they are keyword only, NOT positional.

        Args:
            sma: Semi-major axis in kilometers.
            ecc: Eccentricity of the orbit.
            inc: Inclination of the orbit in degrees.
            raan: Right-ascension of ascending node in degrees.
            aop: Argument of periapsis in degrees.
            meanAnomaly: Mean anomaly at epoch in degrees.
            epoch: Epoch of meanAnomaly (Default = None)
        """

        self._sma = sma
        self._ecc = ecc
        self._inc = radians(inc)
        self._raan = radians(raan)
        self._aop = radians(aop)
        self._meanAnomaly = radians(meanAnomaly)
        self._epoch = epoch
        self._tle = None

    @classmethod
    def fromTle(cls, tle: TwoLineElement, time: JulianDate = None):
        """
        Class method used to instantiate an object from a two-line element object at a given time.

        Args:
            tle: Two-Line element of the object to obtain orbital elements of.
            time: Julian Date of the time to compute orbital elements of.
        """

        rtn = cls(**cls.__tleToElements(tle, time))
        rtn._tle = tle
        return rtn

    @classmethod
    def fromState(cls, position: EVector, velocity: EVector, time: JulianDate = None):
        """
        Class method used to instantiate an object from a known state at a given time.

        Args:
            position: Position state of the object at epoch.
            velocity: Velocity state of the object at epoch.
            time: Julian Date of the time to compute orbital elements of.
        """

        elements = cls.__stateToElements(position, velocity, time)
        return cls(**elements)

    def __str__(self) -> str:
        """Returns a string representation of the orbital elements."""

        return (f'Semi-major axis: {self._sma}, Eccentricity: {self._ecc}, Inclination: {degrees(self._inc)}\n'
                + f'RAAN: {degrees(self._raan)}, Argument of Perigee: {degrees(self._aop)}, '
                  f'Mean Anomaly: {degrees(self._meanAnomaly)}\nEpoch: {str(self._epoch)}')

    @classmethod
    def __tleToElements(cls, tle: TwoLineElement, time: JulianDate) -> dict:
        """Logic for computing elements from TLE."""
        dt = time.difference(tle.getEpoch())
        args = {'inc': tle.getInc()}
        inc = radians(tle.getInc())
        n0 = tle.getMeanMotion()
        dM = (dt * (n0 + dt * (tle.getMeanMotionDot() + dt * tle.getMeanMotionDDot()))) * 360.0
        args['meanAnomaly'] = (tle.getMeanAnomaly() + dM) % 360.0
        ecc0 = tle.getEcc()
        n0dot = tle.getMeanMotionDot() * 2  # todo: ensure this is multiplied by 2 (if not move it up before dM)
        a0 = tle.getSma()
        aDot = -2 * a0 * n0dot / (3 * n0)
        args['sma'] = a0 + aDot * dt
        eDot = -2 * (1 - ecc0) * n0dot / (3 * n0)
        args['ecc'] = ecc0 + eDot * dt
        temp = (a0 ** -3.5) / ((1 - (ecc0 * ecc0)) ** 2)

        #   Perturbations
        #   Non-spherical Earth
        lanJ2Dot = -2.06474e14 * temp * cos(inc)
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * ((sin(inc)) ** 2))
        #   Third-Body
        lanMoon = -0.00338 * cos(inc) / n0
        lanSun = -0.00154 * cos(inc) / n0
        aopMoon = 0.00169 * (4 - (5 * sin(inc) ** 2)) / n0
        aopSun = 0.00077 * (4 - (5 * sin(inc) ** 2)) / n0

        args['raan'] = (tle.getRaan() + (lanJ2Dot + lanMoon + lanSun) * dt) % 360.0
        args['aop'] = (tle.getAop() + (aopJ2Dot + aopMoon + aopSun) * dt) % 360.0
        args['epoch'] = time
        return args

    @classmethod
    def __stateToElements(cls, position: EVector, velocity: EVector, time: JulianDate) -> dict:
        """Logic for computing elements from a known state."""
        angMom = cross(position, velocity)
        lineOfNodes = norm(cross(EVector(0, 0, 1), angMom))
        eccVec = computeEccentricVector(position, velocity)

        # todo: test if removing the lineOfNodes.mag() messes anything up (value should be 1)
        args = {'ecc': eccVec.mag(), 'inc': degrees(acos(angMom[2] / angMom.mag()))}
        raan = degrees(acos(lineOfNodes[0] / lineOfNodes.mag()))
        args['raan'] = 360.0 - raan if lineOfNodes[1] < 0 else raan
        aop = degrees(acos(dot(lineOfNodes, eccVec) / (lineOfNodes.mag() * eccVec.mag())))
        args['aop'] = 360.0 - aop if eccVec[2] < 0 else aop
        tAnom = acos(dot(eccVec, position) / (eccVec.mag() * position.mag()))
        if dot(position, velocity) < 0:
            tAnom = 2 * pi - tAnom
        args['meanAnomaly'] = degrees(trueToMean(tAnom, args['ecc']))
        args['sma'] = (angMom.mag() ** 2) / ((1 - (args['ecc'] * args['ecc'])) * EARTH_MU)
        args['epoch'] = time
        return args

    def __setFromDict(self, elements: dict):
        """Sets the internal values from the computed dictionary."""
        self._sma = elements['sma']
        self._ecc = elements['ecc']
        self._inc = radians(elements['inc'])
        self._raan = radians(elements['raan'])
        self._aop = radians(elements['aop'])
        self._meanAnomaly = radians(elements['meanAnomaly'])
        self._epoch = elements['epoch']

    def update(self, time: JulianDate):
        """
        Updates elements to their values at the given time. For an object created with a TLE, this will
        reevaluate the values as if initializing with the given time. For an object created from a known state
        or with individual values, only the mean anomaly can be adjusted, as well as the epoch.

        Args:
            time: The time to update the orbital elements to.
        """

        if self._tle is not None:
            elements = self.__tleToElements(self._tle, time)
            self.__setFromDict(elements)
        else:
            # todo: adjust mean anomaly
            pass  # cannot recompute elements, so adjust mean anomaly

    def updateFromState(self, position: EVector, velocity: EVector, time: JulianDate = None):
        """
        Updates current object's elements from another known state. This method is to keep
        from creating a new object when a more up-to-date state is known.

        Args:
            position: Position state vector in kilometers.
            velocity: Velocity state vector in kilometers / second.
            time: Epoch of the state (Default = None).
        """

        elements = self.__stateToElements(position, velocity, time)
        self.__setFromDict(elements)

    def getState(self, time: JulianDate = None) -> tuple[EVector]:
        """
        Computes state vectors based on the orbital elements at a given time.

        Args:
            time: Time to compute the state vectors.

        Returns:
            A tuple containing the position and velocity vectors in kilometers and kilometers / second
            respectively.
        """

        if time is None:
            tAnom = meanToTrue(self._meanAnomaly, self._ecc)
            elements = {'sma': self._sma, 'inc': self._inc, 'ecc': self._ecc, 'raan': self._raan, 'aop': self._aop}
        else:
            elements = self.__tleToElements(self._tle, time)
            elements['inc'] = radians(elements['inc'])
            elements['raan'] = radians(elements['raan'])
            elements['aop'] = radians(elements['aop'])
            tAnom = meanToTrue(radians(elements['meanAnomaly']), elements['ecc'])

        eAnom = trueToEccentric(tAnom, elements['ecc'])
        r = elements['sma'] * (1 - elements['ecc'] * cos(eAnom))
        pOrbit = EVector(cos(tAnom), sin(tAnom), 0) * r
        vOrbit = EVector(-sin(eAnom), sqrt(1 - elements['ecc'] * elements['ecc']) * cos(eAnom), 0) * (
                sqrt(EARTH_MU * elements['sma']) / r)
        rot = getEulerMatrix(ZXZ, EulerAngles(elements['raan'], elements['inc'], elements['aop']))
        return rot @ pOrbit, rot @ vOrbit

    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        """
        Computes the perifocal reference frame of the satellite with the instance's elements. If a time
        is provided and a TLE was used to initialize, the reference frame will be computed with up-to-date
        elements computed from the TLE.

        Args:
            time: Time to adjust elements used to compute the reference frame (Default = None).

        Returns:
            A ReferenceFrame object corresponding to the elements.
        """

        if self._tle and time is not None:
            elements = self.__tleToElements(self._tle, time)
            return ReferenceFrame(ZXZ, EulerAngles(
                radians(elements['raan']),
                radians(elements['inc']),
                radians(elements['aop'])
            ))
        else:
            return ReferenceFrame(ZXZ, EulerAngles(
                self._raan,
                self._inc,
                self._aop
            ))

    def meanAnomalyAt(self, time: JulianDate) -> float:
        """
        Computes the most accurate mean anomaly at a given epoch. If object was instantiated with
        a TLE, the mean motion derivative are used in a simplified approximation for perturbations.
        Other instances are propagated forward with constant mean motion.

        Args:
            time: Epoch to find the mean anomaly.

        Returns:
            The mean anomaly in radians.

        Raises ValueError: If an epoch has not been set.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        if self._tle is not None:
            dt = time.difference(self._tle.getEpoch())
            dM = (dt * (self._tle.getMeanMotion() + dt * (self._tle.getMeanMotionDot() + dt *
                                                          self._tle.getMeanMotionDDot()))) * TWOPI
            return (radians(self._tle.getMeanAnomaly()) + dM) % TWOPI
        #   non-TLE computation
        dt = time.difference(self._epoch)
        dM = smaToMeanMotion(self._sma) * dt * 86400.0
        return (self._meanAnomaly + dM) % TWOPI

    def trueAnomalyAt(self, time: JulianDate) -> float:
        """
        Computes the most accurate true anomaly at a given epoch. If the object was instantiated with
        a TLE, the mean motion derivatives are used in a simplified approximation for perturbations.
        Other instances are propagated forward with constant mean motion.

        Args:
            time: Epoch to find the true anomaly.

        Returns:
            The true anomaly in radians.

        Raises ValueError: If an epoch has not been set.
        """

        return meanToTrue(self.meanAnomalyAt(time), self._ecc)

    def getEccentricVector(self) -> EVector:
        """
        Computes the eccentric vector of the current elements.

        Raises ValueError: If an epoch was not set for the class.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return computeEccentricVector(*self.getState(self._epoch))

    def setSma(self, sma: float):
        """Sets the semi-major axis of the orbit in kilometers."""
        self._sma = sma

    def getSma(self) -> float:
        """Returns the semi-major axis of the orbit in kilometers."""
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
        """Sets the right-ascension of the ascending node in radians. Synonymous with longitude of ascending node."""
        self._raan = raan

    def getRaan(self):
        """Returns the right-ascension of the ascending node in radians. Synonymous with longitude of ascending node."""
        return self._raan

    def setAop(self, aop: float):
        """Sets the argument of periapsis measured in radians."""
        self._aop = aop

    def getAop(self):
        """Returns the argument of periapsis measured in radians."""
        return self._aop

    def setMeanAnomaly(self, meanAnom: float):
        """Sets the mean anomaly of the orbit in radians. For the positional methods to work properly the epoch should
        be changed whenever this value is changed."""
        self._meanAnomaly = meanAnom

    def getMeanAnomaly(self):
        """Returns the mean anomaly of the orbit in radians."""
        return self._meanAnomaly

    def setEpoch(self, epoch: JulianDate):
        """Sets the epoch associated with the mean anomaly. For the positional methods to work properly, the mean
        anomaly should be changed whenever this value is changed."""
        self._epoch = epoch

    def getEpoch(self):
        """Returns the epoch associated with the objects true anomaly."""
        return self._epoch


def raanProcessionRate(tle: TwoLineElement) -> float:
    """
    Computes the procession rate of the right-ascension of the ascending node due to third-body non-spherical Earth
    perturbations. Values are negative for a prograde (direct) orbit, positive for retrograde.

    Args:
        tle: Two-line elements set for a satellite.

    Returns:
        The rate of procession of the RAAN in radians per day.
    """

    a = meanMotionToSma(tle.getMeanMotion())
    dRaanSphere = -2.06474e14 * (a ** -3.5) * cos(tle.getInc()) / (
            (1 - tle.getEcc() * tle.getEcc()) ** 2)
    dRaanMoon = -0.00338 * cos(tle.getInc()) / tle.getMeanMotion()
    dRaanSun = -0.00154 * cos(tle.getInc()) / tle.getMeanMotion()
    return radians(dRaanSphere + dRaanMoon + dRaanSun)


def aopProcessionRate(tle: TwoLineElement) -> float:
    """
    Computes the procession rate of the argument of periapsis due to third-body and non-spherical Earth perturbations.
    Values are positive for a prograde (direct) orbit, negative for retrograde.

    Args:
        tle: Two-line element set for a satellite.

    Returns:
        The rate of procession of the argument of periapsis in radians per day.
    """

    a = meanMotionToSma(tle.getMeanMotion())
    dAopSphere = 1.03237e14 * (a ** -3.5) * (4 - 5 * (sin(tle.getInc()) ** 2)) / \
        ((1 - tle.getEcc() * tle.getEcc()) ** 2)
    dAopMoon = 0.00169 * (4 - 5 * (sin(tle.getInc()) ** 2)) / tle.getMeanMotion()
    dAopSun = 0.00077 * (4 - 5 * (sin(tle.getInc()) ** 2)) / tle.getMeanMotion()
    return radians(dAopSphere + dAopMoon + dAopSun)


def computeEccentricVector(position: EVector, velocity: EVector) -> EVector:
    """
    Computes the eccentric vector of an orbit, which points in the positive x-axis of it's perifocal reference frame.

    Args:
        position: Position state vector of the satellite.
        velocity: Velocity state vector of the satellite.

    Returns:
        The eccentric vector of the satellite's orbit.
    """

    lhs = position * (velocity.mag2() / EARTH_MU - 1 / position.mag())
    rhs = velocity * (dot(position, velocity) / EARTH_MU)
    return lhs - rhs


def computeVelocity(sma: float, ecc: float, trueAnom: float) -> float:
    """
    Computes the velocity of a satellite at a given true anomaly.

    Args:
        sma: Semi-major axis in kilometers.
        ecc: Eccentricity of the orbit.
        trueAnom: True anomaly in radians.

    Returns:
        The magnitude of the velocity of the satellite in kilometers / second.
    """

    radius = computeRadius(sma, ecc, trueAnom)
    return sqrt(EARTH_MU * ((2 / radius) - (1 / sma)))


def computeRadius(sma: float, ecc: float, trueAnom: float) -> float:
    """
    Computes thr radius of a satellite at a given true anomaly.

    Args:
        sma: Semi-major axis in kilometers
        ecc: Eccentricity of the orbit.
        trueAnom: True anomaly in radians.

    Returns:
        The radius of the satellite position in kilometers.
    """

    return sma * (1 - ecc * ecc) / (1 + ecc * cos(trueAnom))


def computeFlightAngle(ecc: float, trueAnom: float) -> float:
    """
    Computes the flight angle of a satellite at a given true anomaly.

    Args:
        ecc: Eccentricity of the orbit.
        trueAnom: True anomaly in radians.

    Returns:
        The flight angle of the satellite in radians.
    """

    return atan((ecc * sin(trueAnom)) / (1 + ecc * cos(trueAnom)))
