from pyevspace import EVector

from sattrack.spacetime.sidereal import earthOffsetAngle
from sattrack.spacetime.juliandate import JulianDate
from sattrack.util.constants import EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS

'''
What does a Body need?
mass (MU)
equitorial radius, polar radius, flattening
elements
parent Body
rotation period (sidereal vs solar?)
name
reference frame (rotating)

what do we need for interfacing?
way to compute position (relative to parent)
    - simply propagating with mean motion or using a more accurate model 
compute reference frame offset (i.e. sidereal time)
radius at a latitude

'''


class Body:
    """
    A class for representing a celestial body. The class holds the most significant values of a body needed for various
    calculations.
    """

    def __init__(self, name: str, mu: float, re: float, revPeriod: float, *, rp: float = 0, parent=None):
        """
        Initializes attributes to the arguments.

        Args:
            name: Name of the body.
            mu: Gravitational parameter of the body in km^3s^-2
            re: Equitorial radius in kilometers.
            revPeriod: Revolution period (sidereal day length) in seconds.
            rp: Polar radius in kilometers (Default = 0). The flattening ratios is set based on this value.
            parent: Parent body of this body.
        """

        self._name = name
        self._MU = mu
        self._radiusEq = re
        self._radiusPl = re if rp == 0 else rp
        self._flattening = (re - rp) / re
        self._parent = parent
        self._revPeriod = revPeriod

    def __str__(self):
        """Returns a string representation of the body."""
        rtn = f'{self._name}\nMU: {"%.4f" % self._MU} km^3s^-2\nEquitorial radius: {"%.3f" % self._radiusEq} km\n' \
              f'Polar radius: {"%.3f" % self._radiusPl} km\nSidereal day: {self._revPeriod} seconds'
        if self._parent is not None:
            rtn += f'\nParent: {self._parent.getName()}'
        return rtn

    def getName(self) -> str:
        """Returns the name of the body."""
        return self._name

    def getMu(self) -> float:
        """Returns the gravitational parameter of the body."""
        return self._MU

    def getEquitorialRadius(self) -> float:
        """Returns the equitorial radius in kilometers."""
        return self._radiusEq

    def getPolarRadius(self) -> float:
        """Returns the polar radius in kilometers."""
        return self._radiusPl

    def getFlattening(self) -> float:
        """Returns the flattening rate of the body if it is not a perfect sphere."""
        return self._flattening

    def getParent(self):
        """Returns the parent body this body orbits."""
        return self._parent

    def getSiderealPeriod(self) -> float:
        """Returns the time of a sidereal day in seconds."""
        return self._revPeriod

    def getOffsetAngle(self, time: JulianDate) -> float:
        """Method for computing the offset angle of the body."""
        pass

    def getPosition(self, time: JulianDate) -> EVector:
        """Method for computing body position, should be set for each individual body."""
        pass

    # other methods here


'''_earth_elements = OrbitalElements(sma=149.598e9, ecc=0.0167, inc=0.0, raan=348.73936, aop=102.94719,
                                  meanAnomaly=100.46435, epoch=J2000)'''
# todo: make a sun body and set as earth's parent
EARTH_BODY = Body("Earth", EARTH_MU, EARTH_EQUITORIAL_RADIUS, 86164.090531, rp=EARTH_POLAR_RADIUS, parent=None)
EARTH_BODY.getOffsetAngle = earthOffsetAngle
