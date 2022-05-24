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
        """Sets the body the satellite orbits. This really only be used when satellites switch
        bodies after escaping one body."""
        self._body = body


class Satellite(_SGP4_Propagator):
    """This object describes a satellite that can be modeled with an SGP type propagation model,
    namely the SGP4 or SDP4, described by the NORAD two-line element set.
    The class also contains an internal buffer holding the most recent state and time associated
    with it."""

    def __init__(self, tle: TwoLineElement):
        super().__init__(tle)
        self._name = tle.getName()
        self._recent_time = None
        self._recent_state = (None, None)

    def getState(self, jd: JulianDate):
        """Computes the state vectors for a satellite at a given time. The Satellite object should
        be loaded with an up-to-date TLE as the accuracy becomes less the farther from the TLE
        epoch the desired time is.
        Parameters:
        jd: The time to find the state vectors of the satellite.
        Returns a tuple with the position and velocity vectors as the first and second elements."""
        self._recent_time = jd
        state = super().getState(jd)
        self._recent_state = state
        return state

    def name(self) -> str:
        """Returns the name of the satellite, which is determined by the name in the TLE."""
        return self._name

    def tle(self) -> TwoLineElement:
        """Returns the most recent TLE loaded into the satellite. If the TLE hasn't been updated
        by using updateTle(), this is the TLE used when instantiating the Satellite object."""
        return self._tle

    def epoch(self) -> JulianDate:
        """The epoch found in the most recently loaded TLE."""
        return self._tleEpoch

    def updateTle(self, tle: TwoLineElement) -> None:
        """Updates the TLE to a more up-to-date TLE than the one loaded when instantiating."""
        self._name = tle.getName()
        self._tle = tle
        self._tleEpoch = tle.epoch()
        self.initialize(tle)

    def recentTime(self) -> JulianDate:
        """The time the most recent state was calculated."""
        return self._recent_time

    def recentState(self):
        """The most recent state calculated."""
        return self._recent_state
