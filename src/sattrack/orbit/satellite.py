from abc import abstractmethod, ABC
from copy import deepcopy
from typing import TYPE_CHECKING

from pyevspace import Vector, ZXZ, Z_AXIS, Angles, getMatrixEuler, rotateMatrixFrom, getMatrixAxis, \
    ReferenceFrame

from sattrack.orbit.sgp4 import getState
from sattrack.bodies.earth import Earth
from sattrack.orbit.elements import elementsFromTle, Elements, radiusAtPeriapsis, radiusAtApoapsis, meanToTrueAnomaly, \
    radiusAtAnomaly, flightAngleAtAnomaly, velocityAtAnomaly, nextMeanAnomaly, previousMeanAnomaly, nextTrueAnomaly, \
    previousTrueAnomaly, nearestTrueAnomaly, nearestMeanAnomaly, meanAnomalyAtTime, trueAnomalyAtTime, smaToMeanMotion
from sattrack.bodies.topocentric import toTopocentricState
from sattrack.core.juliandate import now
from sattrack.util.constants import TWOPI, SECONDS_PER_DAY, SIDEREAL_PER_SOLAR

if TYPE_CHECKING:
    from sattrack.orbit.sgp4 import TwoLineElement
    from sattrack.bodies.body import Body
    from sattrack.core.coordinates import GeoPosition
    from sattrack.core.juliandate import JulianDate


class Orbitable(ABC):
    """Abstract object to represent an orbitable object around a parent body. This provides the similar attributes and
    methods of the Orbit and Satellite classes. An Orbit object is an Orbitable created from an Elements object, and
    a Satellite object is an Orbitable created with a TwoLineElement object."""
    __slots__ = '_name', '_body'

    def __init__(self, name: str, body: 'Body' = Earth):
        """Initialize the orbitable with the common requirements name and parent body, which defaults to EARTH_BODY."""

        self._name = name
        self._body = body

    @property
    def name(self):
        """Returns the name of the orbitable object."""
        return self._name

    @property
    def body(self):
        """Returns the parent body of the orbitable object."""
        return self._body

    @abstractmethod
    def anomalyAtTime(self, time: 'JulianDate', anomalyType: str = 'true') -> float:
        pass

    @abstractmethod
    def timeToNextAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        pass

    @abstractmethod
    def timeToPreviousAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        pass

    @abstractmethod
    def timeToNearestAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        pass

    @abstractmethod
    def getState(self, time: 'JulianDate') -> (Vector, Vector):
        pass

    def getTopocentricState(self, geo: 'GeoPosition', time: 'JulianDate') -> (Vector, Vector):
        """Computes the state vectors and rotates them to a topocentric reference frame."""

        state = self.getState(time)
        return toTopocentricState(*state, geo, time)

    @abstractmethod
    def getElements(self, time: 'JulianDate') -> Elements:
        pass

    @abstractmethod
    def getReferenceFrame(self, time: 'JulianDate' = None) -> ReferenceFrame:
        pass

    @abstractmethod
    def getPeriapsis(self, time: 'JulianDate' = None) -> float:
        pass

    @abstractmethod
    def getApoapsis(self, time: 'JulianDate' = None) -> float:
        pass

    @abstractmethod
    def isGeosynchronous(self) -> bool:
        """Some package utilities may need to distinguish between geo- and non-geosynchronous satellites."""
        pass


class Orbit(Orbitable):
    """A class derived from the Orbitable abstract class that describes an orbitable object defined by a set of
    orbital elements around a parent body."""

    __slots__ = '_elements', '_periapsis', '_apoapsis'

    def __init__(self, elements: Elements, name: str = '', body: 'Body' = Earth):
        """Initializes an orbit from a set of orbital elements. Future and past positions can't be computed if the
        epoch attribute was not set when instantiating and elements parameter."""

        super().__init__(name, body)
        self._elements = elements
        self._periapsis = radiusAtPeriapsis(self._elements.sma, self._elements.ecc)
        self._apoapsis = radiusAtApoapsis(self._elements.sma, self._elements.ecc)

    def anomalyAtTime(self, time: 'JulianDate', anomalyType: str = 'true') -> float:
        """Finds the anomaly in radians at a given time. Anomaly type must be 'true' or 'mean'."""

        meanMotion = smaToMeanMotion(self._elements.sma, self._body.mu)
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            t0 = meanToTrueAnomaly(self._elements.meanAnomaly, self._elements.ecc)
            return trueAnomalyAtTime(meanMotion, self._elements.ecc, t0, self._elements.epoch, time)
        elif anomalyArg == 'mean':
            return meanAnomalyAtTime(meanMotion, self._elements.meanAnomaly, self._elements.epoch, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToNextAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the next time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'."""

        meanMotion = smaToMeanMotion(self._elements.sma, self._body.mu) * SECONDS_PER_DAY / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            return nextTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly, self._elements.epoch,
                                   anomaly, time)
        elif anomalyArg == 'mean':
            return nextMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToPreviousAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the previous time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'."""

        meanMotion = smaToMeanMotion(self._elements.sma, self._body.mu) * SECONDS_PER_DAY / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            return previousTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly,
                                       self._elements.epoch, anomaly, time)
        elif anomalyArg == 'mean':
            return previousMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToNearestAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the nearest time the orbitable achieves anomaly in radians. Anomaly types are 'true' and 'mean'. The
        time parameter here is not used, and can be None."""

        meanMotion = smaToMeanMotion(self._elements.sma, self._body.mu) * 86400 / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            return nearestTrueAnomaly(meanMotion, self._elements.ecc, self._elements.trueAnomaly, self._elements.epoch,
                                      anomaly)
        elif anomalyArg == 'mean':
            return nearestMeanAnomaly(meanMotion, self._elements.meanAnomaly, self._elements.epoch, anomaly)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def getState(self, time: 'JulianDate') -> (Vector, Vector):
        """Returns the state vectors of the orbitable at a given time in kilometers and kilometers / second."""

        radius = radiusAtAnomaly(self._elements.sma, self._elements.ecc, self._elements.trueAnomaly)
        flightAngle = flightAngleAtAnomaly(self._elements.ecc, self._elements.trueAnomaly)
        velocityMag = velocityAtAnomaly(self._elements.sma, radius, self._body.mu)

        angs = Angles(self._elements.raan, self._elements.inc, self._elements.aop + self._elements.trueAnomaly)
        rotationMatrix = getMatrixEuler(ZXZ, angs)
        position = rotateMatrixFrom(rotationMatrix, Vector.e1) * radius

        # since velocity is the flight angle away from 90 degrees from the position vector, rotating the first reference
        # frame by the flight angle and taking the y-axis gives the velocity direction
        rotationMatrix = rotationMatrix @ getMatrixAxis(Z_AXIS, -flightAngle)
        velocity = rotateMatrixFrom(rotationMatrix, Vector.e2) * velocityMag

        return position, velocity

    def getElements(self, time: 'JulianDate') -> Elements:
        """Returns the orbital elements of the orbitable adjusted by the given time."""

        elements = deepcopy(self._elements)
        meanMotion = smaToMeanMotion(self._elements.sma, self._body.mu)
        meanAnomaly = meanAnomalyAtTime(meanMotion, self._elements.meanAnomaly, self._elements.epoch, time)
        elements.updateAnomaly('mean', meanAnomaly, time)

        return elements

    def getReferenceFrame(self, time: 'JulianDate' = None) -> ReferenceFrame:
        """Computes the perifocal reference frame of the orbitable at a given time. The y-axis is the eccentric vector,
        which points towards periapsis, the z-axis points toward the angular momentum vector and the x-axis completes
        the right-handed coordinate system."""

        angles = Angles(self._elements.raan, self._elements.inc, self._elements.aop)
        return ReferenceFrame(ZXZ, angles)

    def getPeriapsis(self, time: 'JulianDate' = None) -> float:
        """Returns the periapsis, the lowest altitude, of the orbit."""

        return self._periapsis

    def getApoapsis(self, time: 'JulianDate' = None) -> float:
        """Returns the apoapsis, the highest altitude, of the orbit."""

        return self._apoapsis

    def isGeosynchronous(self) -> bool:
        tmp = 1 / SIDEREAL_PER_SOLAR
        period = tmp * TWOPI / SECONDS_PER_DAY
        delta = period * 0.01

        elements = self._elements
        meanMotion = smaToMeanMotion(elements.sma, self.body.mu)

        return period - delta <= meanMotion <= period + delta

    # todo: test this
    def impulse(self, time: 'JulianDate', prograde: float, normal: float, radial: float):
        """Adjusts the orbit of an orbitable based on an impulse taking place at a given instance. Prograde is positive
        in the direction of motion, normal is positive in the direction of the angular momentum vector and radial is
        positive away from the parent body. All delta-v values are in kilometers / second."""

        position, velocity = self.getState(time)
        # x - radial out, y - prograde, z - normal
        flightAngle = flightAngleAtAnomaly(self._elements.ecc, self._elements.trueAnomaly)
        angs = Angles(self._elements.raan, self._elements.inc, self._elements.aop +
                      self._elements.trueAnomaly - flightAngle)
        referenceFrame = ReferenceFrame(ZXZ, angs)
        progradeVector = referenceFrame.rotateFrom(Vector.e2)
        radialVector = referenceFrame.rotateFrom(Vector.e1)
        normalVector = referenceFrame.rotateFrom(Vector.e3)

        newVelocity = velocity + progradeVector * prograde + radialVector * radial + normalVector * normal
        self._elements = Elements.fromState(position, newVelocity, self._elements.epoch)
        self._periapsis = radiusAtPeriapsis(self._elements.sma, self._elements.ecc)
        self._apoapsis = radiusAtApoapsis(self._elements.sma, self._elements.ecc)

    @property
    def elements(self):
        # Return original elements used to instantiate the object.
        return self._elements


class Satellite(Orbitable):
    """Derived from the Orbitable class, a Satellite object is an orbitable described by a TwoLineElement object around
    a parent body."""
    __slots__ = '_tle'

    def __init__(self, tle: 'TwoLineElement'):
        """Initializes the satellite object with a TLE."""

        super().__init__(tle.name, Earth)
        self._tle = tle

    @property
    def tle(self):
        return self._tle

    def anomalyAtTime(self, time: 'JulianDate', anomalyType: str = 'true') -> float:
        """Finds the anomaly in radians at a given time. Anomaly type must be 'true' or 'mean'."""

        _, _, _, eccentricity, _, meanAnomaly = elementsFromTle(self._tle, time)

        if anomalyType == 'true':
            return meanToTrueAnomaly(meanAnomaly, eccentricity)
        elif anomalyType == 'mean':
            return meanAnomaly
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToNextAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the next time the satellite achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""

        _, _, _, eccentricity, sma, meanAnomaly = elementsFromTle(self._tle, time)
        meanMotion = smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            trueAnomaly = meanToTrueAnomaly(meanAnomaly, eccentricity)
            return nextTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly, time)
        elif anomalyArg == 'mean':
            return nextMeanAnomaly(meanMotion, meanAnomaly, time, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToPreviousAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the previous time the orbitable achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""

        _, _, _, eccentricity, sma, meanAnomaly = elementsFromTle(self._tle, time)
        meanMotion = smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            trueAnomaly = meanToTrueAnomaly(meanAnomaly, eccentricity)
            return previousTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly, time)
        elif anomalyArg == 'mean':
            return previousMeanAnomaly(meanMotion, meanAnomaly, time, anomaly, time)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def timeToNearestAnomaly(self, anomaly: float, time: 'JulianDate', anomalyType: str = 'true') -> 'JulianDate':
        """Find the nearest time the orbitable achieves anomaly in radians. Anomaly types are 'true' or 'mean'."""

        _, _, _, eccentricity, sma, meanAnomaly = elementsFromTle(self._tle, time)
        meanMotion = smaToMeanMotion(sma, self._body.mu) * 86400 / TWOPI
        anomalyArg = anomalyType.lower()

        if anomalyArg == 'true':
            trueAnomaly = meanToTrueAnomaly(meanAnomaly, eccentricity)
            return nearestTrueAnomaly(meanMotion, eccentricity, trueAnomaly, time, anomaly)
        elif anomalyArg == 'mean':
            return nearestMeanAnomaly(meanMotion, meanAnomaly, time, anomaly)
        else:
            raise ValueError("anomalyType must be 'mean' or 'true'", anomalyType)

    def getState(self, time: 'JulianDate') -> (Vector, Vector):
        """Computes the position and velocity state vectors in kilometers and kilometers / second from the SGP4
        model."""

        return getState(self._tle, time)

    def getElements(self, time: 'JulianDate') -> Elements:
        """Returns a set of stable orbital elements at a specified time. For oscillating elements see
        Elements.fromState()."""

        return Elements.fromTle(self._tle, time)

    def getReferenceFrame(self, time: 'JulianDate' = None) -> ReferenceFrame:
        """Generates a pyevspace.ReferenceFrame object representing the perifocal reference frame of the satellite's
        orbit."""

        raan, inc, aop, *_ = elementsFromTle(self._tle, time)
        return ReferenceFrame(ZXZ, Angles(raan, inc, aop))

    def getPeriapsis(self, time: 'JulianDate' = None) -> float:
        """Computes the perigee (the smallest altitude) of the satellites orbit relative to the Earth's equitorial
        radius."""

        # todo: do we use the apogee and perigee from the SGP4 or do we compute it from oscillating elements w.r.t. time
        return self._tle.perigee

    def getApoapsis(self, time: 'JulianDate' = None) -> float:
        """Computes the apogee (the highest altitude) of the satellites orbit relative to the Earth's equitorial
        radius."""

        return self._tle.apogee

    def isGeosynchronous(self) -> bool:
        tmp = 1 / SIDEREAL_PER_SOLAR
        period = tmp * TWOPI / SECONDS_PER_DAY
        delta = period * 0.01

        elements = self.getElements(now())
        meanMotion = smaToMeanMotion(elements.sma, self.body.mu)

        return period - delta <= meanMotion <= period + delta
