from elements import OrbitalElements
from body import Body, EARTH_BODY
from juliandate import JulianDate


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
