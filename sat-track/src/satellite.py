from math import asin, degrees

from _sgp4 import _SGP4_Propagator

from anomalies import atan2
from coordinates import GeoPosition
from elements import OrbitalElements
from body import Body, EARTH_BODY
from spacetime import JulianDate
from tle import TwoLineElement
from pyevspace import EVector

from topos import ijkToSEZ


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


class Satellite:
    """This object describes a satellite that can be modeled with an SGP type propagation model,
    namely the SGP4 or SDP4, described by the NORAD two-line element set.
    The class also contains an internal buffer holding the most recent state and time associated
    with it."""

    def __init__(self, tle: TwoLineElement, body: Body = EARTH_BODY):
        self._propagator = _SGP4_Propagator(tle)
        self._tle = tle
        self._tleEpoch = tle.epoch()
        self._name = tle.getName()
        self._recent_time = None
        self._recent_state = (None, None)
        self._body = body
        self.getState(tle.epoch())

    def __str__(self):
        return f'Name: {self._name}\nTLE:\n{self._tle.getLine1()}\n{self._tle.getLine2()}\nRecent Time: ' \
               f'{self._recent_time}\nRecent State:\nPosition: {self._recent_state[0]}\nVelocity: ' \
               f'{self._recent_state[1]}'

    def getState(self, jd: JulianDate):
        """Computes the state vectors for a satellite at a given time. The Satellite object should
        be loaded with an up-to-date TLE as the accuracy becomes less the farther from the TLE
        epoch the desired time is.
        Parameters:
        jd: The time to find the state vectors of the satellite.
        Returns a tuple with the position and velocity vectors as the first and second elements."""
        self._recent_time = jd
        state = self._propagator.getState(self._tle, jd)
        pos = EVector(state[0][0] * 1000.0, state[0][1] * 1000.0, state[0][2] * 1000.0)
        vel = EVector(state[1][0] * 1000.0, state[1][1] * 1000.0, state[1][2] * 1000.0)
        self._recent_state = (pos, vel)
        return self._recent_state

    def getToposPosition(self, jd: JulianDate, geo: GeoPosition) -> EVector:
        if self._recent_time == jd:
            return ijkToSEZ(self._recent_state[0], jd, geo)
        else:
            self.getState(jd)
            return ijkToSEZ(self._recent_state[0], jd, geo)

    def getAltitude(self, jd: JulianDate, geo: GeoPosition) -> float:
        if self._recent_time == jd:
            sez = ijkToSEZ(self._recent_state[0], jd, geo)
        else:
            self.getState(jd)
            sez = ijkToSEZ(self._recent_state[0], jd, geo)
        return degrees(asin(sez[2] / sez.mag()))

    def getAzimuth(self, jd: JulianDate, geo: GeoPosition) -> float:
        if self._recent_time == jd:
            sez = ijkToSEZ(self._recent_state[0], jd, geo)
        else:
            self.getState(jd)
            sez = ijkToSEZ(self._recent_state[0], jd, geo)
        return degrees(atan2(sez[1], -sez[0]))

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

    def updateTLE(self, tle: TwoLineElement) -> None:
        """Updates the TLE to a more up-to-date TLE than the one loaded when instantiating."""
        self._name = tle.getName()
        self._tle = tle
        self._tleEpoch = tle.epoch()
        self._propagator = _SGP4_Propagator(tle)

    def recentTime(self) -> JulianDate:
        """The time the most recent state was calculated."""
        return self._recent_time

    def recentState(self):
        """The most recent state calculated."""
        return self._recent_state

    def getBody(self) -> Body:
        return self._body

    def setBody(self, body: Body):
        self._body = body
