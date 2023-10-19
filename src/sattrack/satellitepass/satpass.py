import json
from math import asin, degrees, atan2

from _pyevspace import dot

from sattrack.eclipse import getShadowTimes, Shadow
from sattrack.exceptions import NoSatelliteEclipseException, OrbitPathAlwaysUp, SatelliteAlwaysAbove, NoPassException
from sattrack.orbit.orbitpath import OrbitPath, TIME_DIFFERENCE, computeEllipseVectors
from sattrack.satellitepass.info import PositionInfo, Visibility
from sattrack.sun import getSunTimes
from sattrack.util.conversions import atan3

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.orbit.satellite import Orbitable
    from sattrack.coordinates import GeoPosition
    from sattrack.spacetime.juliandate import JulianDate


class SatellitePass:

    __slots__ = '_name', '_infos', '_visibility'

    def __init__(self, infos: list[PositionInfo], name: str = ''):
        self._name = name
        self._infos = {
            'rise': min(infos, key=lambda o: o.time),
            'set': max(infos, key=lambda o: o.time),
            'max': max(infos, key=lambda o: o.altitude),
        }

        illuminated = [info for info in infos if info.illuminated]
        unobscured = [info for info in infos if info.unobscured]
        visible = [info for info in infos if info.visible]
        self._visibility = Visibility(bool(illuminated), bool(unobscured))

        # Set first, last, and max visible infos.
        self._infos['maxVisible'] = max(visible, key=lambda o: o.altitude, default=None)
        self._infos['firstVisible'] = min(visible, key=lambda o: o.time, default=None)
        self._infos['lastVisible'] = max(visible, key=lambda o: o.time, default=None)

        # Set first, last, and max illuminated infos.
        self._infos['maxIlluminated'] = max(illuminated, key=lambda o: o.altitude, default=None)
        self._infos['firstIlluminated'] = min(illuminated, key=lambda o: o.time, default=None)
        self._infos['lastIlluminated'] = max(illuminated, key=lambda o: o.time, default=None)

        # Set first, last, and max unobscured infos.
        self._infos['maxUnobscured'] = max(unobscured, key=lambda o: o.altitude, default=None)
        self._infos['firstUnobscured'] = min(unobscured, key=lambda o: o.time, default=None)
        self._infos['lastUnobscured'] = max(unobscured, key=lambda o: o.time, default=None)

    def _getInfos(self) -> list[(str, PositionInfo)]:
        rtn = [
            ('rise', self._infos['rise']),
            ('set', self._infos['set']),
            ('max', self._infos['max']),
        ]

        basicInfos = [j for _, j in rtn]
        if self._infos['firstUnobscured'] not in basicInfos:
            rtn.append(('first unobscured', self._infos['firstUnobscured']))
        if self._infos['lastUnobscured'] not in basicInfos:
            rtn.append(('last unobscured', self._infos['lastUnobscured']))
        if self._infos['firstIlluminated'] not in basicInfos:
            rtn.append(('first illuminated', self._infos['firstIlluminated']))
        if self._infos['lastIlluminated'] not in basicInfos:
            rtn.append(('last illuminated', self._infos['lastIlluminated']))

        return rtn

    def __str__(self):
        width = 97
        if self._name == '':
            title = '{:^{w}}'.format('Pass details at {}'.format(self._infos['max'].time.date()), w=width)
        else:
            title = '{:^{w}}'.format('Pass details for {}, at {}'.format(self._name, self._infos['max'].time.date()),
                                     w=width)
        heading = ' {:^17} | {:^12} | altitude | {:^12} | illuminated | unobscured | visible ' \
            .format('instance', 'time', 'azimuth')
        lineBreak = '-' * width
        string = title + '\n' + heading + '\n'

        infos = [i for i in self._getInfos() if i[1] is not None]
        infos.sort(key=lambda o: o[1].time)

        for name, info in infos:
            string += '{}\n{}\n'.format(lineBreak, str(info).format(name))
        return string

    def __repr__(self):
        return self.__str__()

    def toJson(self) -> str:
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self) -> dict:
        return {"name": self._name, "riseInfo": self._infos['rise'], "setInfo": self._infos['set'],
                "maxInfo": self._infos['max'], "firstUnobscured": self._infos['firstUnobscured'],
                "lastUnobscured": self._infos['lastUnobscured'], "firstIlluminated": self._infos['firstIlluminated'],
                "lastIlluminated": self._infos['lastIlluminated'], "illuminated": self._visibility.illuminated,
                "unobscured": self._visibility.unobscured, "visible": self._visibility.visible}

    @property
    def name(self) -> str:
        return self._name

    @property
    def visibility(self) -> Visibility:
        return self._visibility

    @property
    def illuminated(self) -> bool:
        return self._visibility.illuminated

    @property
    def unobscured(self) -> bool:
        return self._visibility.unobscured

    @property
    def visible(self) -> bool:
        return self._visibility.visible

    @property
    def riseInfo(self) -> PositionInfo:
        return self._infos['rise']

    @property
    def setInfo(self) -> PositionInfo:
        return self._infos['set']

    @property
    def maxInfo(self) -> PositionInfo:
        return self._infos['max']

    @property
    def firstVisibleInfo(self) -> PositionInfo:
        return self._infos['firstVisible']

    @property
    def lastVisibleInfo(self) -> PositionInfo:
        return self._infos['lastVisible']

    @property
    def maxVisibleInfo(self) -> PositionInfo:
        return self._infos['maxVisible']

    @property
    def firstUnobscuredInfo(self) -> PositionInfo:
        return self._infos['firstUnobscured']

    @property
    def lastUnobscuredInfo(self) -> PositionInfo:
        return self._infos['lastUnobscured']

    @property
    def maxUnobscuredInfo(self) -> PositionInfo:
        return self._infos['maxUnobscured']

    @property
    def firstIlluminatedInfo(self) -> PositionInfo:
        return self._infos['firstIlluminated']

    @property
    def lastIlluminatedInfo(self) -> PositionInfo:
        return self._infos['lastIlluminated']

    @property
    def maxIlluminatedInfo(self) -> PositionInfo:
        return self._infos['maxIlluminated']


class PassController:
    __slots__ = '_sat', '_geo', '_path'

    def __init__(self, sat: 'Orbitable', geo: 'GeoPosition'):
        self._sat = sat
        self._geo = geo
        self._path = OrbitPath(sat, geo)

    def _getMaxTime(self, riseTime: 'JulianDate', setTime: 'JulianDate') -> 'JulianDate':
        # todo: this isn't 100% accurate for every sat, but very accurate for circular LEO ones.
        time = riseTime.future((setTime - riseTime) / 2)
        prevTime = time.future(-1)

        while abs(time - prevTime) > TIME_DIFFERENCE:
            a, b, c = computeEllipseVectors(self._sat, time)
            zeta = self._geo.getZenithVector(time)
            t = atan2(dot(zeta, b), dot(zeta, a))
            prevTime = time
            time = self._sat.timeToNearestAnomaly(t, time, 'true')

        return time

    def _getNextPass(self, time: 'JulianDate', nextOccurrence: bool) -> SatellitePass:

        satTimes = self._path.computeSatellitePassTimes(time, nextOccurrence)
        # # computeSatellitePassTimes will return an empty tuple if there is no pass.
        # if not satTimes:
        #     return ()
        riseTime, setTime = satTimes

        # Some satellites never enter the
        try:
            enterTime, exitTime = getShadowTimes(self._sat, riseTime, Shadow.PENUMBRA)
        except NoSatelliteEclipseException:
            enterTime = setTime.future(0.0001)
            exitTime = enterTime
        sunRiseTime, sunSetTime = getSunTimes(riseTime, self._geo)

        topoState = self._sat.getTopocentricState(self._geo, riseTime)
        riseInfo = PositionInfo(0.0, degrees(atan3(topoState[0][1], -topoState[0][0])), riseTime,
                                riseTime < enterTime, riseTime < sunRiseTime or riseTime > sunSetTime)
        topoState = self._sat.getTopocentricState(self._geo, setTime)
        setInfo = PositionInfo(0.0, degrees(atan3(topoState[0][1], -topoState[0][0])), setTime,
                               setTime < enterTime or setTime > exitTime,
                               setTime < sunRiseTime or setTime > sunSetTime)

        maxTime = self._getMaxTime(riseTime, setTime)
        topoState = self._sat.getTopocentricState(self._geo, maxTime)
        maxInfo = PositionInfo(degrees(asin(topoState[0][2] / topoState[0].mag())),
                               degrees(atan3(topoState[0][1], -topoState[0][0])),
                               riseTime.future((setTime - riseTime) / 2),
                               maxTime < enterTime or maxTime > exitTime,
                               maxTime < sunRiseTime or maxTime > sunSetTime)

        firstIlluminatedInfo = lastIlluminatedInfo = firstUnobscuredInfo = lastUnobscuredInfo = None
        if riseTime < exitTime < setTime:
            topoState = self._sat.getTopocentricState(self._geo, exitTime)
            firstIlluminatedInfo = PositionInfo(degrees(asin(topoState[0][2] / topoState[0].mag())),
                                                degrees(atan3(topoState[0][1], -topoState[0][0])), exitTime, True,
                                                exitTime < sunRiseTime or exitTime >= sunSetTime)
        if riseTime < enterTime < setTime:
            topoState = self._sat.getTopocentricState(self._geo, enterTime)
            lastIlluminatedInfo = PositionInfo(degrees(asin(topoState[0][2] / topoState[0].mag())),
                                               degrees(atan3(topoState[0][1], -topoState[0][0])), enterTime, True,
                                               enterTime < sunRiseTime or enterTime >= sunSetTime)
        if riseTime < sunSetTime < setTime:
            topoState = self._sat.getTopocentricState(self._geo, sunSetTime)
            firstUnobscuredInfo = PositionInfo(degrees(asin(topoState[0][2] / topoState[0].mag())),
                                               degrees(atan3(topoState[0][1], -topoState[0][0])), sunSetTime,
                                               sunSetTime < enterTime or sunSetTime >= exitTime, True)
        if riseTime < sunRiseTime < setTime:
            topoState = self._sat.getTopocentricState(self._geo, sunRiseTime)
            lastUnobscuredInfo = PositionInfo(degrees(asin(topoState[0][2] / topoState[0].mag())),
                                              degrees(atan3(topoState[0][1], -topoState[0][0])), sunRiseTime,
                                              sunRiseTime < enterTime or sunRiseTime >= exitTime, True)

        infos = [riseInfo, setInfo, maxInfo] + \
                [info for info in [firstIlluminatedInfo, lastIlluminatedInfo,
                                   firstUnobscuredInfo, lastUnobscuredInfo] if info is not None]

        return SatellitePass(infos)

    @staticmethod
    def _nextTime(o: ('JulianDate', 'JulianDate')) -> 'JulianDate':
        return o[1].future(0.0001)

    @staticmethod
    def _prevTime(o: ('JulianDate', 'JulianDate')) -> 'JulianDate':
        return o[0].future(-0.0001)

    def getNextPass(self, time: 'JulianDate', number: int = 1) -> SatellitePass:
        if number > 0:
            # If time is a SatellitePass start after set time.
            if isinstance(time, SatellitePass):
                time = time.setInfo.time.future(0.0001)

            # Get times directly.
            if number == 1:
                return self._getNextPass(time, True)

            nextOccurrence = True
            key = self._nextTime

        elif number < 0:
            # If time is a SatellitePass start before rise time.
            if isinstance(time, SatellitePass):
                time = time.riseInfo.time.future(-0.0001)

            # Get times directly.
            if number == -1:
                return self._getNextPass(time, False)

            nextOccurrence = False
            key = self._prevTime

        else:
            raise ValueError(f'number must be a non-zero integer, not {number}')

        tmp = (time, time)
        for i in range(abs(number) - 1):
            tmp = self._path.computeSatellitePassTimes(key(tmp), nextOccurrence)

        return self._getNextPass(key(tmp), nextOccurrence)

    def getPassList(self, time: 'JulianDate', length: int) -> list[SatellitePass]:
        lastTime = time.future(length)
        passList = []
        nextPass = self._getNextPass(time, True)
        time = nextPass.riseInfo.time

        while time < lastTime:
            passList.append(nextPass)
            try:
                nextPass = self._getNextPass(time.future(0.0001), True)
            # If we get an exception in the middle of getting passes return the list we have up until that point.
            except (OrbitPathAlwaysUp, SatelliteAlwaysAbove, NoPassException):
                return passList
            time = nextPass.riseInfo.time

        return passList

    @property
    def orbitable(self) -> 'Orbitable':
        return self._sat

    @property
    def geo(self) -> 'GeoPosition':
        return self._geo

    @property
    def path(self) -> 'OrbitPath':
        return self._path
