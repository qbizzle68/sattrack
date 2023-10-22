from inspect import signature, Parameter
from typing import Callable, TYPE_CHECKING

from pyevspace import Vector

from sattrack.util.constants import TWOPI, SUN_MU, SUN_RADIUS, EARTH_MU, EARTH_EQUITORIAL_RADIUS,\
    EARTH_SIDEREAL_PERIOD, EARTH_POLAR_RADIUS
from sattrack.bodies.sun import getSunPosition
from sattrack.core.sidereal import siderealTime

if TYPE_CHECKING:
    from sattrack.orbit.satellite import Orbit
    from sattrack.core.juliandate import JulianDate


def _checkCallable(func: Callable, errorName: str):
    """Makes sure an object is a callable with 1 parameter."""

    if isinstance(func, Callable):
        params = signature(func).parameters
        if len(params) != 1:
            ls = [p for p in params.values() if p.default is Parameter.empty and p.kind in
                  (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)]
            if len(ls) != 1:
                raise TypeError(f'{errorName} must be a Callable with 1 non-default positional parameter')
    else:
        raise TypeError(f'{errorName} must be a Callable type', type(func))


class Body:
    __slots__ = '_name', '_MU', '_Re', '_Rp', '_flattening', '_revPeriod', '_orbit', '_parent', \
                '_offsetFunc', '_positionFunc'

    def __init__(self, name: str, MU: float, Re: float, revPeriod: float, offsetFunc: Callable, positionFunc: Callable,
                 *, Rp: float = 0, orbit: 'Orbit' = None, parent: 'Body' = None):
        """Creates an object representing a celestial body like a plant or moon. Any unknown or non-applicable
        parameters should be set to None."""

        _checkCallable(offsetFunc, 'offsetFunc')
        _checkCallable(positionFunc, 'positionFunc')
        if orbit is not None:
            self._orbit = orbit
        if parent is not None:
            self._parent = parent

        self._name = name
        self._MU = MU
        self._Re = Re
        if not Rp:
            self._flattening = (Re - Rp) / Re
            self._Rp = Rp
        else:
            self._flattening = 0
            self._Rp = Re
        self._revPeriod = revPeriod
        self._offsetFunc = offsetFunc
        self._positionFunc = positionFunc

    @property
    def name(self):
        """Returns the name of the celestial body."""
        return self._name

    @property
    def mu(self):
        """Returns the gravitational parameter of the celestial body."""
        return self._MU

    @property
    def equitorialRadius(self):
        """Returns the equitorial radius of the celestial body."""
        return self._Re

    @property
    def polarRadius(self):
        """Returns the polar radius of the celestial body."""
        return self._Rp

    @property
    def flattening(self):
        """Returns the flattening ratio of the celestial body (Re - Rp) / Re."""
        return self._flattening

    @property
    def period(self):
        """Returns the sidereal revolution period of the celestial body."""
        return self._revPeriod

    @property
    def orbit(self):
        """Returns the Orbit of the celestial object, or None if one wasn't set."""
        return self._orbit

    @property
    def parent(self):
        """Returns the parent body if it was set, otherwise None."""
        return self._parent

    def getOffsetAngle(self, time: 'JulianDate') -> float:
        """Returns the offset angle of the prime-meridian of the body from the celestial x-axis at a given time."""

        return self._offsetFunc(time)

    def getPosition(self, time: 'JulianDate') -> Vector:
        """Returns the position of the center of the celestial body relative to the Earth's center."""

        return self._positionFunc(time)

    def getAngularMomentum(self):
        """Returns the angular momentum of the rotation of the body."""

        return TWOPI / self._revPeriod


"""Global Body object representing the Sun."""
SUN_BODY = Body('Sun', SUN_MU, SUN_RADIUS, 0, lambda time: 0, getSunPosition)
"""Global Body object representing the Earth."""
EARTH_BODY = Body('Earth', EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_SIDEREAL_PERIOD, siderealTime,
                  lambda time: Vector(0, 0, 0), Rp=EARTH_POLAR_RADIUS, parent=SUN_BODY)
