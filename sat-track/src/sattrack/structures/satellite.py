from math import radians, pi

from pyevspace import EVector

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import ReferenceFrame, EulerAngles
from sattrack.structures.sgp4 import SGP4_Propagator
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.body import Body, EARTH_BODY
from sattrack.structures.elements import OrbitalElements
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import trueToMean, timeToNextMeanAnomaly, timeToPrevMeanAnomaly
from sattrack.util.conversions import smaToMeanMotion


class Satellite:

    def __init__(self, obj: TwoLineElement | OrbitalElements, body: Body = EARTH_BODY):
        """
        Initializes the satellite from a TLE or a set of orbital elements.
        Args:
            obj: A TwoLineElement or a set of orbital elements.
            body: Parent body of the satellite (Default = EARTH_BODY).
        """

        # todo: add name to the sat (need a name from elements)
        self._tle = None
        self._propagator = None
        self._epoch = None
        self._elements = None
        self._periapsisPassage = None
        self.update(obj)
        self._body = body

    def __str__(self) -> str:
        # todo: create a string with name and either tle or elements, plus periapsis passage and body
        pass

    # todo: add a stack to store states, with bool parameter that defaults to false
    def getState(self, jd: JulianDate) -> tuple[EVector]:
        """
        Computes a set of state vectors for the satellite at a given time.

        Args:
            jd: Time to compute the state.

        Returns:
            A tuple of EVectors that contain the position and velocity in kilometers and kilometers / second
            respectively.
        """

        if self._tle:
            rawState = self._propagator.getState(self._tle, jd)
            return EVector(rawState[0][0], rawState[0][1], rawState[0][2]), EVector(rawState[1][0], rawState[1][1],
                                                                                    rawState[1][2])
        else:
            return self._elements.getState(jd)

    def update(self, obj: TwoLineElement | OrbitalElements):
        """
        Updates the TLE or orbital elements of the satellite.

        Args:
            obj: TwoLineElement or OrbitalElements that describe the satellite.
        """

        # todo: is computing periapsis passage still ideal?
        if type(obj) == TwoLineElement:
            self._tle = obj
            self._propagator = SGP4_Propagator(obj)
            self._epoch = obj.getEpoch()
            self._elements = None
            dt = -radians(obj.getMeanAnomaly()) / (2 * pi * obj.getMeanMotion())
            self._periapsisPassage = self._epoch.future(dt)
        elif type(obj) == OrbitalElements:
            self._tle = None
            self._propagator = None
            self._epoch = obj.getEpoch()
            self._elements = obj
            if obj.getEpoch() is not None:
                dt = -radians(obj.getMeanAnomaly()) / smaToMeanMotion(obj.getSma())
                self._periapsisPassage = self._epoch.future(dt)
            else:
                self._periapsisPassage = None

    def getTle(self) -> TwoLineElement:
        """
        Returns the TLE of the satellite.

        Returns:
            The TLE of the satellite.

        Raises ValueError: If a TLE was not set for the satellite.
        """

        if not self._tle:
            raise ValueError('No TLE was set.')
        return self._tle

    def getElements(self) -> OrbitalElements:
        """
        Returns either the orbital elements used to initialize the satellite, or computes a set of elements from the
        stored TLE.

        Returns:
            A set of orbital elements for the satellite at epoch.
        """

        return self._elements if self._elements else OrbitalElements.fromTle(self._tle)

    def getEpoch(self) -> JulianDate:
        """
        Returns the epoch of the satellite if one was set.

        Returns:
            The epoch of the satellite.

        Raises ValueError: If no epoch was set.
        """

        if not self._epoch:
            raise ValueError('No epoch was set.')
        return self._epoch

    def setBody(self, body: Body):
        """Sets the parent body of the satellite."""
        self._body = body

    def getBody(self) -> Body:
        """Returns the parent body of the satellite."""
        return self._body

    # todo: do we need to include this?
    def getPeriapsisPassage(self) -> float:
        """Returns the time of periapsis passage for the satellite."""
        return self._periapsisPassage

    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        """Computes the rotation object for the perifocal reference frame."""
        if time is None:
            return ReferenceFrame(
                ZXZ,
                EulerAngles(
                    radians(self._tle.getRaan()),
                    radians(self._tle.getInc()),
                    radians(self._tle.getAop())
                )
            )
        else:
            elements = OrbitalElements.fromTle(self._tle, time)
            return ReferenceFrame(
                ZXZ,
                EulerAngles(
                    elements.getRaan(),
                    elements.getInc(),
                    elements.getAop()
                )
            )

    def getEccentricVector(self, time: JulianDate = None) -> EVector:
        """
        Computes the eccentric vector of the satellite.

        Args:
            time: Can be used to compute up-to-date orbital elements for a more accurate computation.

        Returns:
            The eccentric vector for the satellite.
        """

        # todo: figure out why we cant compute this from the state.
        if time is None:
            rot = ReferenceFrame(ZXZ, EulerAngles(
                    radians(self._tle.getRaan()),
                    radians(self._tle.getInc()),
                    radians(self._tle.getAop())
                ))
            return rot.RotateFrom(EVector.e1) * self._tle.getEcc()
        else:
            elements = OrbitalElements.fromTle(self._tle, time)
            rot = ReferenceFrame(ZXZ, EulerAngles(
                    elements.getRaan(),
                    elements.getInc(),
                    elements.getAop()
                ))
            return rot.RotateTo(EVector.e1) * elements.getEcc()

    def timeToNextMeanAnomaly(self, mAnom: float, time: JulianDate) -> JulianDate:
        """
        Computes the next time the satellite passes through the mean anomaly after the given time.

        Args:
            mAnom: Mean anomaly in radians.
            time: Relative time to find the next anomaly.

        Returns:
            Time to the next mean anomaly occurrence after the time given.
        """

        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        return timeToNextMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            mAnom,
            time
        )

    def timeToPrevMeanAnomaly(self, mAnom: float, time: JulianDate) -> JulianDate:
        """
        Compute the previous time the satellite passes through the mean anomaly before the given time.

        Args:
            mAnom: Mean anomaly in radians.
            time: Relative time to find the previous anomaly.

        Returns:
            Time to the previous mean anomaly occurrence before the time given.
        """

        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        return timeToPrevMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            mAnom,
            time
        )

    def timeToNextTrueAnomaly(self, tAnom: float, time: JulianDate) -> JulianDate:
        """
        Compute the next time the satellite passes through the true anomaly after the given time.

        Args:
            tAnom: True anomaly in radians.
            time: Relative time to find the next anomaly.

        Returns:
            Time to the next true anomaly occurrence after the time given.
        """

        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        return timeToNextMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            trueToMean(tAnom, elements.getEcc()),
            time
        )

    def timeToPrevTrueAnomaly(self, tAnom: float, time: JulianDate) -> JulianDate:
        """
        Compute the previous time the satellite passes through the true anomaly before the given time.

        Args:
            tAnom: True anomaly in radians.
            time: Relative time to find the previous anomaly.

        Returns:
            Time to the previous true anomaly occurrence before the time given.
        """

        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        return timeToPrevMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            trueToMean(tAnom, elements.getEcc()),
            time
        )
