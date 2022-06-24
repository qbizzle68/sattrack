from elements import OrbitalElements
from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS


class Body:

    def __init__(self, name: str, orbit: OrbitalElements, mu: float, re: float, rev_period: float, offset: float = 0,
                 epoch: JulianDate = None, *, rp: float = 0, parent=None):
        """Initial values needed do describe an orbitable body.
        Parameters:
        name:   Name of the body.
        orbit:  Classical orbital parameters estimating the bodies orbit.
        mu:     The gravitational parameter of the body (G * M) UNITS.
        re:     Equitorial radius in meters.
        rev_period: Revolution period, a.k.a. sidereal day in seconds.
        offset: Offset from prime celestial zero longitude in degrees.
        epoch:  Epoch of the offset, in JulianDays.
        rp:     Polar radius in meters, if different than equitorial radius.
        parent: Parent body, None by default."""
        self._name = name
        self._orbit = orbit
        self._MU = mu
        self._radius_eq = re
        self._radius_pl = re if rp == 0 else rp
        self._flattening = (re - rp) / re
        self._parent = parent
        self._rev_period = rev_period
        self._offset = offset
        self._epoch = epoch

    def __str__(self):
        pass

    def getOffset(self, time: float) -> float:
        pass

    def name(self) -> str:
        return self._name

    def orbit(self) -> OrbitalElements:
        return self._orbit

    def mu(self) -> float:
        return self._MU

    def equitorialRadius(self) -> float:
        return self._radius_eq

    def polarRadius(self) -> float:
        return self._radius_pl

    def flattening(self) -> float:
        return self._flattening

    def parent(self):
        return self._parent

    def siderealPeriod(self) -> float:
        return self._rev_period

    def offsetAngle(self) -> float:
        return self._offset

    def offsetEpoch(self) -> JulianDate:
        return self._epoch

    # other methods here


_earth_elements = OrbitalElements(sma=149.598e9,ecc=0.0167,inc=0.0,raan=348.73936,aop=102.94719,meanAnomaly=100.46435,epoch=J2000)
#_earth_elements = OrbitalElements(149.598e9, 0.0167, 0.0, 348.73936, 102.94719, 100.46435, J2000)
EARTH_BODY = Body("Earth", _earth_elements, EARTH_MU, EARTH_EQUITORIAL_RADIUS, 86164.090531, 280.46061837, J2000,
             rp=EARTH_POLAR_RADIUS, parent=None)  # set parent to sun when created
