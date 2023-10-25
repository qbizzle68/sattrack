import json

from sattrack.bodies.topocentric import AltAz

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate


class Visibility:
    __slots__ = '_illuminated', '_unobscured', '_visible'

    # visible argument is redundant
    def __init__(self, illuminated: bool, unobscured: bool):
        self._illuminated = illuminated
        self._unobscured = unobscured
        self._visible = illuminated and unobscured

    @property
    def illuminated(self):
        return self._illuminated

    @property
    def unobscured(self):
        return self._unobscured

    @property
    def visible(self):
        return self._visible


class PositionInfo:
    __slots__ = '_altAz', '_time', '_visibility'

    def __init__(self, altitude: float, azimuth: float, time: 'JulianDate', illuminated: bool, unobscured: bool):
        self._altAz = AltAz(altitude, azimuth)
        self._visibility = Visibility(illuminated, unobscured)
        self._time = time

    def __str__(self):
        return ' {{:^17}} | {:^12} | {:^{altW}.{p}f} | {:^{azw}.{p}f} {:^5} | {!s:^11} | {!s:^10} | {!s:^7} ' \
            .format(self._time.time(), self._altAz.altitude, self._altAz.azimuth, '({})'.format(self._altAz.direction),
                    bool(self._visibility.illuminated), bool(self._visibility.unobscured),
                    bool(self._visibility.visible), altW=8, azw=6, p=2)

    def __repr(self):
        args = (self._altAz.altitude, self._altAz.azimuth, self._time, self._visibility.illuminated,
                self._visibility.unobscured)
        return f'PositionInfo{str(args)}'

    def toJson(self):
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        return {"altitude": self._altAz.altitude, "azimuth": self._altAz.azimuth, "direction": self._altAz.direction,
                "time": self._time.toDict(), "illuminated": self._visibility.illuminated,
                "unobscured": self._visibility.unobscured, "visible": self._visibility.visible}

    @property
    def altAz(self):
        return self._altAz

    @property
    def altitude(self):
        return self._altAz.altitude

    @property
    def azimuth(self):
        return self._altAz.azimuth

    @property
    def direction(self):
        return self._altAz.direction

    @property
    def time(self):
        return self._time

    @property
    def visibility(self):
        return self._visibility

    @property
    def illuminated(self):
        return self._visibility.illuminated

    @property
    def unobscured(self):
        return self._visibility.unobscured

    @property
    def visible(self):
        return self._visibility.visible
