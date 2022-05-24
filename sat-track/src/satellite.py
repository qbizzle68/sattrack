from _sgp4 import _SGP4_Propagator
from elements import OrbitalElements
from body import Body, EARTH_BODY
from spacetime import JulianDate
from tle import TwoLineElement


class SimpleSatellite(OrbitalElements):
    """Simple satellite that has a constant mean motion and whose orbit can be described by a set of classical
    orbital elements. The orbit, therefore, of this type of satellite is modeled by a conic section, and is usually
    over idealized, but can be used to more simply compute orbital characteristics/positions."""

    def __init__(self, name: str, sma: float, ecc: float, inc: float, raan: float, aop: float, trueAnomaly: float,
                 *, epoch: JulianDate = 0, body: Body = EARTH_BODY):
        self._name = name
        super().__init__(sma, ecc, inc, raan, aop, trueAnomaly, epoch)
        self._body = body

    def __str__(self) -> str:
        return f"Name: {self._name}\n{super().__str__()}"

    def getBody(self) -> Body:
        """Returns the body the satellite orbits."""
        return self._body

    def setBody(self, body: Body):
        """Sets the body the satellite orbits. This really only be used when satellites switch bodies after escaping
        one body."""
        self._body = body


class Satellite(_SGP4_Propagator):

    def __init__(self, tle: TwoLineElement):
        super().__init__(tle)
        self._name = tle.getName()
        self._recent_time = None
        self._recent_state = (None, None)

    def getState(self, jd: JulianDate):
        self._recent_time = jd
        state = super().getState(jd)
        self._recent_state = state
        return state

    def name(self) -> str:
        return self._name

    def tle(self) -> TwoLineElement:
        return self._tle

    def epoch(self) -> JulianDate:
        return self._tleEpoch

    def updateTle(self, tle: TwoLineElement) -> None:
        self._name = tle.getName()
        self._tle = tle
        self._tleEpoch = tle.epoch()
        self.initialize(tle)

    def recentTime(self) -> JulianDate:
        return self._recent_time

    def recentState(self):
        return self._recent_state
