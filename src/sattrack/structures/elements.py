import math as _math

from pyevspace import EVector, dot, cross, norm

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import getEulerMatrix, EulerAngles, ReferenceFrame
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.body import Body, EARTH_BODY
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import meanToTrue, trueToMean, trueToEccentric, timeToNextMeanAnomaly, \
    timeToPrevMeanAnomaly, meanAnomalyAt
from sattrack.util.constants import TWOPI
from sattrack.util.conversions import meanMotionToSma, smaToMeanMotion


class OrbitalElements:
    """A container class that computes and contains orbital elements of a satellite, either imported manually, or
    derived from an approximation of TLE values at a given time. The mean anomaly is the anomaly used for positioning,
    and an epoch can be set to compute times for future/past anomalies. The initialized angles are in degrees, for ease
    of use, but the angles are converted and implemented as radians."""

    def __init__(self, *, sma: float, ecc: float, inc: float, raan: float, aop: float, meanAnomaly: float,
                 name: str = "", epoch: JulianDate = None, body: Body = EARTH_BODY):
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
            name: Name of the object.
            epoch: Epoch of meanAnomaly (Default = None).
            body: Body of the orbiting satellite (Default = EARTH_BODY).
        """

        self._sma = sma
        self._ecc = ecc
        self._inc = _math.radians(inc)
        self._raan = _math.radians(raan)
        self._aop = _math.radians(aop)
        self._meanAnomaly = _math.radians(meanAnomaly)
        self._name = name
        self._epoch = epoch
        self._tle = None
        self._body = body

    @classmethod
    def fromTle(cls, tle: TwoLineElement, time: JulianDate = None):
        """
        Class method used to instantiate an object from a two-line element object at a given time.

        Args:
            tle: Two-Line element of the object to obtain orbital elements of.
            time: Julian Date of the time to compute orbital elements of.
        """

        rtn = cls(**cls.__tleToElements(tle, time), body=EARTH_BODY)
        rtn._tle = tle
        rtn._name = tle.getName()
        return rtn

    @classmethod
    def fromState(cls, position: EVector, velocity: EVector, *, time: JulianDate = None, body: Body = EARTH_BODY):
        """
        Class method used to instantiate an object from a known state at a given time.

        Args:
            position: Position state of the object at epoch.
            velocity: Velocity state of the object at epoch.
            time: Julian Date of the time to compute orbital elements of.
            body: Body of the orbiting satellite (Default = EARTH_BODY).
        """

        elements = cls.__stateToElements(position, velocity, time, body)
        return cls(**elements, body=body)

    def __str__(self) -> str:
        """Returns a string representation of the orbital elements."""

        return f'Name: {self._name}\nSemi-major axis: {self._sma}, Eccentricity: {self._ecc}, Inclination: ' \
            f'{_math.degrees(self._inc)}\nRAAN: {_math.degrees(self._raan)}, Argument of Perigee: {_math.degrees(self._aop)}, ' \
            f'Mean Anomaly: {_math.degrees(self._meanAnomaly)}\nEpoch: {str(self._epoch)}'

    @classmethod
    def __tleToElements(cls, tle: TwoLineElement, time: JulianDate) -> dict:
        """Logic for computing elements from TLE."""
        # dt = time.difference(tle.getEpoch())
        dt = time - tle.getEpoch()
        args = {'inc': tle.getInc()}
        inc = _math.radians(tle.getInc())
        n0 = tle.getMeanMotion()
        dM = (dt * (n0 + dt * (tle.getMeanMotionDot() + dt * tle.getMeanMotionDDot()))) * 360.0
        args['meanAnomaly'] = (tle.getMeanAnomaly() + dM) % 360.0
        ecc0 = tle.getEcc()
        n0dot = tle.getMeanMotionDot() * 2
        a0 = tle.getSma()
        aDot = -2 * a0 * n0dot / (3 * n0)
        args['sma'] = a0 + aDot * dt
        eDot = -2 * (1 - ecc0) * n0dot / (3 * n0)
        args['ecc'] = ecc0 + eDot * dt
        temp = (a0 ** -3.5) / ((1 - (ecc0 * ecc0)) ** 2)

        #   Perturbations
        #   Non-spherical Earth
        lanJ2Dot = -2.06474e14 * temp * _math.cos(inc)
        aopJ2Dot = 1.03237e14 * temp * (4 - 5 * ((_math.sin(inc)) ** 2))
        #   Third-Body
        lanMoon = -0.00338 * _math.cos(inc) / n0
        lanSun = -0.00154 * _math.cos(inc) / n0
        aopMoon = 0.00169 * (4 - (5 * _math.sin(inc) ** 2)) / n0
        aopSun = 0.00077 * (4 - (5 * _math.sin(inc) ** 2)) / n0

        args['raan'] = (tle.getRaan() + (lanJ2Dot + lanMoon + lanSun) * dt) % 360.0
        args['aop'] = (tle.getAop() + (aopJ2Dot + aopMoon + aopSun) * dt) % 360.0
        args['epoch'] = time
        return args

    @classmethod
    def __stateToElements(cls, position: EVector, velocity: EVector, time: JulianDate = None,
                          body: Body = EARTH_BODY) -> dict:
        """Logic for computing elements from a known state."""
        angMom = cross(position, velocity)
        lineOfNodes = norm(cross(EVector(0, 0, 1), angMom))  # lineOfNodes.mag() not used because it's normalized
        eccVec = computeEccentricVector(position, velocity)

        args = {'ecc': eccVec.mag(), 'inc': _math.degrees(_math.acos(angMom[2] / angMom.mag()))}
        raan = _math.degrees(_math.acos(lineOfNodes[0]))
        args['raan'] = 360.0 - raan if lineOfNodes[1] < 0 else raan
        aop = _math.degrees(_math.acos(dot(lineOfNodes, eccVec) / eccVec.mag()))
        args['aop'] = 360.0 - aop if eccVec[2] < 0 else aop
        tAnom = _math.acos(dot(eccVec, position) / (eccVec.mag() * position.mag()))
        if dot(position, velocity) < 0:
            tAnom = 2 * _math.pi - tAnom
        args['meanAnomaly'] = _math.degrees(trueToMean(tAnom, args['ecc']))
        args['sma'] = (angMom.mag() ** 2) / ((1 - (args['ecc'] * args['ecc'])) * body.getMu())
        args['epoch'] = time
        return args

    def __setFromDict(self, elements: dict):
        """Sets the internal values from the computed dictionary."""
        self._sma = elements['sma']
        self._ecc = elements['ecc']
        self._inc = _math.radians(elements['inc'])
        self._raan = _math.radians(elements['raan'])
        self._aop = _math.radians(elements['aop'])
        self._meanAnomaly = _math.radians(elements['meanAnomaly'])
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
            if self._epoch is None:
                raise ValueError('Epoch was not set for this instance.')
            self._meanAnomaly = meanAnomalyAt(smaToMeanMotion(self._sma), self._meanAnomaly, self._epoch, time)

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
            elements['inc'] = _math.radians(elements['inc'])
            elements['raan'] = _math.radians(elements['raan'])
            elements['aop'] = _math.radians(elements['aop'])
            tAnom = meanToTrue(_math.radians(elements['meanAnomaly']), elements['ecc'])

        eAnom = trueToEccentric(tAnom, elements['ecc'])
        r = elements['sma'] * (1 - elements['ecc'] * _math.cos(eAnom))
        pOrbit = EVector(_math.cos(tAnom), _math.sin(tAnom), 0) * r
        vOrbit = EVector(-_math.sin(eAnom), _math.sqrt(1 - elements['ecc'] * elements['ecc']) * _math.cos(eAnom), 0) * (
                _math.sqrt(self._body.getMu() * elements['sma']) / r)
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
                _math.radians(elements['raan']),
                _math.radians(elements['inc']),
                _math.radians(elements['aop'])
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
            # dt = time.difference(self._tle.getEpoch())
            dt = time - self._tle.getEpoch()
            dM = (dt * (self._tle.getMeanMotion() + dt * (self._tle.getMeanMotionDot() + dt *
                                                          self._tle.getMeanMotionDDot()))) * TWOPI
            return (_math.radians(self._tle.getMeanAnomaly()) + dM) % TWOPI
        #   non-TLE computation
        # dt = time.difference(self._epoch)
        dt = time - self._epoch
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

    def timeToNextMeanAnomaly(self, mAnom: float, time: JulianDate) -> JulianDate:
        """
        Computes the next time a satellite passes through the mean anomaly after the given time.

        Args:
            mAnom: Mean anomaly in radians
            time: Relative time to find the next anomaly.

        Returns:
            The next time the satellite's position is mAnom.

        Raises ValueError: If an epoch was not set.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return timeToNextMeanAnomaly(
            smaToMeanMotion(self._sma),
            self._meanAnomaly,
            self._epoch,
            mAnom,
            time
        )

    def timeToPrevMeanAnomaly(self, mAnom: float, time: JulianDate) -> JulianDate:
        """
        Computes the previous time the satellite passed through the mean anomaly before the given time.

        Args:
            mAnom: Mean anomaly in radians.
            time: Relative time to find the previous anomaly.

        Returns:
            The previous time the satellite's position was mAnom.

        Raises ValueError: If an epoch was not set.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return timeToPrevMeanAnomaly(
            smaToMeanMotion(self._sma),
            self._meanAnomaly,
            self._epoch,
            mAnom,
            time
        )

    def timeToNextTrueAnomaly(self, tAnom: float, time: JulianDate) -> JulianDate:
        """
        Computes the next time the satellite passes through the true anomaly after the time given.

        Args:
            tAnom: True anomaly in radians.
            time: Relative time to find the next anomaly

        Returns:
            The next time the satellite's position is tAnom.

        Raises ValueError: If an epoch was not set.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return timeToNextMeanAnomaly(
            smaToMeanMotion(self._sma),
            self._meanAnomaly,
            self._epoch,
            trueToMean(tAnom, self._ecc),
            time
        )

    def timeToPrevTrueAnomaly(self, tAnom: float, time: JulianDate) -> JulianDate:
        """
        Computes the previous time the satellite passes through the true anomaly before the time given.

        Args:
            tAnom: True anomaly in radians.
            time: Relative time to find the previous anomaly.

        Returns:
            The previous time the satellite position is tAnom.

        Raises ValueError: If an epoch was not set.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return timeToPrevMeanAnomaly(
            smaToMeanMotion(self._sma),
            self._meanAnomaly,
            self._epoch,
            trueToMean(tAnom, self._ecc),
            time
        )

    def getEccentricVector(self) -> EVector:
        """
        Computes the eccentric vector of the current elements.

        Raises ValueError: If an epoch was not set for the class.
        """

        if self._epoch is None:
            raise ValueError('Epoch was not set for this instance.')
        return computeEccentricVector(*self.getState(self._epoch))

    @property
    def sma(self) -> float:
        """Returns the semi-major axis of the orbit in kilometers."""
        return self._sma

    @sma.setter
    def sma(self, value):
        self._sma = value

    @property
    def ecc(self):
        """Returns the eccentricity of the orbit."""
        return self._ecc

    @ecc.setter
    def ecc(self, value):
        self._ecc = value

    @property
    def inc(self):
        """Returns the inclination of the orbit in radians."""
        return self._inc

    @inc.setter
    def inc(self, value):
        self._inc = value

    @property
    def raan(self):
        """Returns the right-ascension of the ascending node in radians. Synonymous with longitude of ascending node."""
        return self._raan

    @raan.setter
    def raan(self, value):
        self._raan = value

    @property
    def aop(self):
        """Returns the argument of periapsis measured in radians."""
        return self._aop

    @aop.setter
    def aop(self, value):
        self._aop = value

    @property
    def meanAnomaly(self):
        """Returns the mean anomaly of the orbit in radians."""
        return self._meanAnomaly

    @meanAnomaly.setter
    def meanAnomaly(self, value):
        self._meanAnomaly = value

    @property
    def name(self) -> str:
        """Returns the name of the object whose elements these represent."""
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def epoch(self):
        """Returns the epoch associated with the objects true anomaly."""
        return self._epoch

    @epoch.setter
    def epoch(self, value):
        self._epoch = value


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
    dRaanSphere = -2.06474e14 * (a ** -3.5) * _math.cos(tle.getInc()) / (
            (1 - tle.getEcc() * tle.getEcc()) ** 2)
    dRaanMoon = -0.00338 * _math.cos(tle.getInc()) / tle.getMeanMotion()
    dRaanSun = -0.00154 * _math.cos(tle.getInc()) / tle.getMeanMotion()
    return _math.radians(dRaanSphere + dRaanMoon + dRaanSun)


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
    dAopSphere = 1.03237e14 * (a ** -3.5) * (4 - 5 * (_math.sin(tle.getInc()) ** 2)) / \
                 ((1 - tle.getEcc() * tle.getEcc()) ** 2)
    dAopMoon = 0.00169 * (4 - 5 * (_math.sin(tle.getInc()) ** 2)) / tle.getMeanMotion()
    dAopSun = 0.00077 * (4 - 5 * (_math.sin(tle.getInc()) ** 2)) / tle.getMeanMotion()
    return _math.radians(dAopSphere + dAopMoon + dAopSun)


def computeEccentricVector(position: EVector, velocity: EVector, body: Body = EARTH_BODY) -> EVector:
    """
    Computes the eccentric vector of an orbit, which points in the positive x-axis of it's perifocal reference frame.

    Args:
        position: Position state vector of the satellite.
        velocity: Velocity state vector of the satellite.
        body: Body of the orbiting satellite (Default = EARTH_BODY).

    Returns:
        The eccentric vector of the satellite's orbit.
    """

    lhs = position * (velocity.mag2() / body.getMu() - 1 / position.mag())
    rhs = velocity * (dot(position, velocity) / body.getMu())
    return lhs - rhs


def computeVelocity(sma: float, ecc: float, trueAnom: float, body: Body = EARTH_BODY) -> float:
    """
    Computes the velocity of a satellite at a given true anomaly.

    Args:
        sma: Semi-major axis in kilometers.
        ecc: Eccentricity of the orbit.
        trueAnom: True anomaly in radians.
        body: Body of the orbiting satellite (Default = EARTH_BODY).

    Returns:
        The magnitude of the velocity of the satellite in kilometers / second.
    """

    radius = computeRadius(sma, ecc, trueAnom)
    return _math.sqrt(body.getMu() * ((2 / radius) - (1 / sma)))


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

    return sma * (1 - ecc * ecc) / (1 + ecc * _math.cos(trueAnom))


def computeFlightAngle(ecc: float, trueAnom: float) -> float:
    """
    Computes the flight angle of a satellite at a given true anomaly.

    Args:
        ecc: Eccentricity of the orbit.
        trueAnom: True anomaly in radians.

    Returns:
        The flight angle of the satellite in radians.
    """

    return _math.atan((ecc * _math.sin(trueAnom)) / (1 + ecc * _math.cos(trueAnom)))
