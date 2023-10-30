from abc import ABC, abstractmethod
from functools import partial
from inspect import getattr_static

from pyevspace import Vector

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate
    from sattrack.core.coordinates import CelestialCoordinates, GeoPosition, AltAz


class BodyOrbitController(ABC):

    @abstractmethod
    def computePosition(self, time: 'JulianDate') -> Vector:
        pass

    @abstractmethod
    def computeCelestialCoordinates(self, time: 'JulianDate') -> 'CelestialCoordinates':
        pass

    @abstractmethod
    def computeAltAz(self, geo: 'GeoPosition', time: 'JulianDate') -> 'AltAz':
        pass


class Body:
    # __slots__ = '_name', '_MU', '_Re', '_Rp', '_flattening', '_revPeriod', '_controller', '_parent'

    def __init__(self, name: str, MU: float, Re: float, revPeriod: float, controller, *,
                 Rp: float = 0, parent: 'Body' = None):
        self._name = name
        self._MU = MU
        self._Re = Re
        self._controller = controller()
        self._parent = parent

        if not Rp:
            Rp = Re
        self._flattening = (Re - Rp) / Re
        self._Rp = Rp
        self._revPeriod = revPeriod

        # Add all attributes that don't start with an underscore from the controller class.
        attributes = [attr for attr in dir(controller) if attr[0] != '_']
        for attrName in attributes:
            # Not sure why but this only works with the _static version of getattr.
            attr = getattr_static(controller, attrName)
            # Need to bind self or cls to the method depending on the function type.
            if isinstance(attr, staticmethod):
                func = attr
            elif isinstance(attr, classmethod):
                # Need to access the __func__ attribute of the classmethod since it's not callable directly.
                func = partial(attr.__func__, controller)
            else:
                func = partial(attr, self._controller)

            setattr(self, attrName, func)

    @property
    def name(self):
        return self._name

    @property
    def mu(self):
        return self._MU

    @property
    def Re(self):
        return self._Re

    @property
    def Rp(self):
        return self._Rp

    @property
    def flattening(self):
        return self._flattening

    @property
    def period(self):
        return self._revPeriod

    @property
    def parent(self):
        return self._parent

    @property
    def controller(self):
        return self._controller
