from pyevspace import EVector
from sattrack.structures._sgp4 import _SGP4_Propagator

from sattrack.spacetime.juliandate import JulianDate
from sattrack.structures.body import Body, EARTH_BODY
from sattrack.structures.elements import OrbitalElements
from sattrack.structures.tle import TwoLineElement


class Satellite:

    def __init__(self, obj: TwoLineElement | OrbitalElements, body: Body = EARTH_BODY):
        if type(obj) == TwoLineElement:
            self._tle = obj
            self._propagator = _SGP4_Propagator(obj)
            self._epoch = obj.epoch()
            self._elements = None
        elif type(obj) == OrbitalElements:
            self._tle = None
            self._propagator = None
            self._epoch = obj.getEpoch()
            self._elements = obj
        self._body = body

    def __str__(self) -> str:
        pass

    def getState(self, jd: JulianDate) -> tuple[EVector]:
        if self._tle:
            rawState = self._propagator.getState(self._tle, jd)
            return EVector(rawState[0][0], rawState[0][1], rawState[0][2]), EVector(rawState[1][0], rawState[1][1], rawState[1][2])
        else:
            return self._elements.getState(jd)

    def update(self, obj: TwoLineElement | OrbitalElements):
        if type(obj) == TwoLineElement:
            self._tle = obj
            self._propagator = _SGP4_Propagator(obj)
            self._epoch = obj.epoch()
            self._elements = None
        elif type(obj) == OrbitalElements:
            self._tle = None
            self._propagator = None
            self._epoch = obj.getEpoch()
            self._elements = obj

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
