from math import radians, pi

from pyevspace import EVector
from sattrack.structures.sgp4 import SGP4_Propagator

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import ReferenceFrame, EulerAngles
from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.body import Body, EARTH_BODY
from sattrack.structures.elements import OrbitalElements
from sattrack.structures.tle import TwoLineElement
from sattrack.util.anomalies import trueToMean, timeToNextMeanAnomaly, timeToPrevMeanAnomaly, meanToTrue
from sattrack.util.conversions import smaToMeanMotion


'''class Satellite2:
    """make with tle, state or elements. tle can do all, state can create elements, elements can be used
    to do basics (mostly find anomalies not accounting for perturbations."""

    def __init__(self, obj: TwoLineElement | tuple[EVector] | OrbitalElements, epoch: JulianDate = None):
        if type(obj) == TwoLineElement:
            self._tle = obj
            self._elements = None
            self._epoch = obj.epoch()
            dt = -radians(obj.meanAnomaly()) / (2 * pi * obj.meanMotion())
            self._periapsisPassage = self._epoch.future(dt)
        elif type(obj) == OrbitalElements:
            self._tle = None
            self._elements = obj
            self._epoch = obj.getEpoch()
            dt = -radians(obj.getMeanAnomaly()) / smaToMeanMotion(obj.getSma())
            self._periapsisPassage = self._epoch.future(dt)
        else:
            pass

    @classmethod
    def fromTle(cls, tle: TwoLineElement):
        pass'''


class Satellite:

    def __init__(self, obj: TwoLineElement | OrbitalElements, body: Body = EARTH_BODY):
        '''if type(obj) == TwoLineElement:
            self._tle = obj
            self._propagator = SGP4_Propagator(obj)
            self._epoch = obj.epoch()
            self._elements = None
            dt = -radians(obj.meanAnomaly()) / (2 * pi * obj.meanMotion())
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
                self._periapsisPassage = None'''
        # todo: add name to the sat (need a name from elements)
        self._tle = None
        self._propagator = None
        self._epoch = None
        self._elements = None
        self._periapsisPassage = None
        self.update(obj)
        self._body = body

    def __str__(self) -> str:
        pass

    # todo: add a stack to store states, with bool parameter that defaults to false
    def getState(self, jd: JulianDate) -> tuple[EVector]:
        if self._tle:
            rawState = self._propagator.getState(self._tle, jd)
            return EVector(rawState[0][0], rawState[0][1], rawState[0][2]), EVector(rawState[1][0], rawState[1][1],
                                                                                    rawState[1][2])
        else:
            return self._elements.getState(jd)

    def update(self, obj: TwoLineElement | OrbitalElements):
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

    def tle(self) -> TwoLineElement:
        if not self._tle:
            raise ValueError('No TLE was set.')
        return self._tle

    def elements(self) -> OrbitalElements:
        return self._elements if self._elements else OrbitalElements.fromTle(self._tle)

    def epoch(self) -> JulianDate:
        if not self._epoch:
            raise ValueError('No epoch was set.')
        return self._epoch

    def setBody(self, body: Body):
        self._body = body

    def getBody(self) -> Body:
        return self._body

    def getPeriapsisPassage(self) -> float:
        return self._periapsisPassage

    def getReferenceFrame(self, time: JulianDate = None) -> ReferenceFrame:
        """Computes the rotation object for the perifocal reference frame """
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
        """Computes the next time the satellite passes through the mean anomaly
        after the time given.
        Parameters:
            mAnom -- mean anomaly in radians
            time -- relative time to find the next anomaly
        Returns the time the satellite's position is mAnom.
        """

        '''twoPi = 2 * pi
        if self._tle:
            n = self._tle.meanMotion()  # rev/day
        else:
            n = smaToMeanMotion(self._elements.getSma()) * 86400 / twoPi
        # number of revolutions since self._periapsisPassage
        revs = time.difference(self._periapsisPassage) * n
        m0 = (revs - floor(revs)) * twoPi
        if mAnom < m0:
            dm = mAnom + twoPi - m0
        else:
            dm = mAnom - m0
        return time.future((dm / n) / twoPi)'''
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
        """Computes the previous time the satellite passed through the mean anomaly
        before the given time.
        Parameters:
            mAnom -- mean anomaly in radians
            time -- relative time to find the previous anomaly
        Returns the time the satellite's position is mAnom.
        """

        '''twoPi = 2 * pi
        if self._tle:
            n = self._tle.meanMotion()  # rev/day
        else:
            n = smaToMeanMotion(self._elements.getSma()) * 86400 / twoPi
        revs = time.difference(self._periapsisPassage) * n
        m0 = (revs - floor(revs)) * twoPi
        if m0 < mAnom:
            dm = mAnom - twoPi - m0
        else:
            dm = mAnom - m0
        return time.future((dm / n) / twoPi)'''
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
        """Computes the next time the satellite passes through the true anomaly
        after the time given.
        Parameters:
            tAnom -- true anomaly in radians
            time -- relative time to find the next anomaly
        Returns the time the satellite's position is tAnom.
        """

        '''if self._tle:
            ecc = self._tle.eccentricity()
        else:
            ecc = self._elements.getEcc()
        return self.timeToNextMeanAnomaly(trueToMean(tAnom, ecc), time)'''
        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        mAnom = timeToNextMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            trueToMean(tAnom, elements.getEcc()),
            time
        )
        return meanToTrue(mAnom, elements.getEcc())

    def timeToPrevTrueAnomaly(self, tAnom: float, time: JulianDate) -> JulianDate:
        """Computes the previous time the satellite passed through the true anomaly
        before the given time.
        Parameters:
            tAnom -- true anomaly in radians
            time -- relative time to find the previous anomaly
        Returns the time the satellite's position is tAnom.
        """

        '''if self._tle:
            ecc = self._tle.eccentricity()
        else:
            ecc = self._elements.getEcc()
        return self.timeToPrevMeanAnomaly(trueToMean(tAnom, ecc), time)'''
        if self._tle:
            elements = OrbitalElements.fromTle(self._tle, time)
        else:
            elements = self._elements
        mAnom = timeToPrevMeanAnomaly(
            smaToMeanMotion(elements.getSma()),
            elements.getMeanAnomaly(),
            elements.getEpoch(),
            trueToMean(tAnom, elements.getEcc()),
            time
        )
        return meanToTrue(mAnom, elements.getEcc())
