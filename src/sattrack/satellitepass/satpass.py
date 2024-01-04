import json
from math import asin, degrees, sqrt, acos

from pyevspace import dot
from pyevspace.core import cross

from sattrack.bodies.sun import Sun
from sattrack.bodies.topocentric import toTopocentric
from sattrack.satellitepass.eclipse import getShadowTimes, Shadow
from sattrack.orbit.exceptions import SatelliteAlwaysAbove, NoPassException, PassedOrbitPathRange
from sattrack.satellitepass.exceptions import NoSatelliteEclipseException
from sattrack.orbit.orbitpath import OrbitPath, computeEllipseVectors
from sattrack.config import TIME_DIFFERENCE
from sattrack.satellitepass.info import PositionInfo, Visibility
from sattrack.util.constants import SECONDS_PER_DAY, MINUTES_PER_DAY, TWOPI, EARTH_ANGULAR_MOMENTUM
from sattrack.util.helpers import atan3, signOf

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.orbit.satellite import Orbitable
    from sattrack.core.coordinates import GeoPosition
    from sattrack.core.juliandate import JulianDate


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
    __slots__ = '_sat', '_geo', '_path', '_ACCURACY_BUFFER'

    def __init__(self, sat: 'Orbitable', geo: 'GeoPosition'):
        self._sat = sat
        self._geo = geo
        self._path = OrbitPath(sat, geo)
        # this value can have issues for a pass that is much less than a
        # minute, but we have to cut it off somewhere
        self._ACCURACY_BUFFER = 10 / SECONDS_PER_DAY

    @property
    def orbitable(self) -> 'Orbitable':
        return self._sat

    @property
    def geo(self) -> 'GeoPosition':
        return self._geo

    @property
    def path(self) -> 'OrbitPath':
        return self._path

    def _computeMaximumAnomaly(self, time):
        a, b, _ = computeEllipseVectors(self._sat, time)
        zeta = self._geo.getZenithVector(time)
        return atan3(dot(zeta, b), dot(zeta, a))

    def _computeCurrentMaximumTime(self, time):

        previousTime = time.future(-1)
        approxTime = time

        # epsilon = 0.1 / SECONDS_PER_DAY
        while abs(approxTime - previousTime) > TIME_DIFFERENCE:
            maxParameter = self._computeMaximumAnomaly(approxTime)
            previousTime = approxTime
            approxTime = self._sat.timeToNearestAnomaly(maxParameter, approxTime, 'true')

        return approxTime

    def _computeNextMaximumTime(self, time: 'JulianDate', nextOccurrence: bool) -> 'JulianDate':
        """time should be a maximum altitude time, otherwise this might not get the next
        occurrence, but the current one again."""

        maxParameter = self._computeMaximumAnomaly(time)
        if nextOccurrence:
            adjustedTime = time.future(self._ACCURACY_BUFFER)
            nextTime = self._sat.timeToNextAnomaly(maxParameter, adjustedTime, 'true')
        else:
            adjustedTime = time.future(-self._ACCURACY_BUFFER)
            nextTime = self._sat.timeToPreviousAnomaly(maxParameter, adjustedTime, 'true')

        nextMaxTime = self._computeCurrentMaximumTime(nextTime)
        if abs(nextMaxTime - time) < (10 / MINUTES_PER_DAY):
            return self._computeNextMaximumTime(nextTime, nextOccurrence)

        return nextMaxTime

    def _computeMaximumTimeExec(self, time: 'JulianDate', number: int, nextOccurrence: bool,
                                maximumSearchPeriod: float) -> 'JulianDate':

        lastPositiveAltitude = time
        if self._sat.getAltitude(self._geo, time) > 0:
            nextTime = self._computeCurrentMaximumTime(time)
        elif number > 0:
            nextTime = self._computeNextMaximumTime(time, True)
        else:
            nextTime = self._computeNextMaximumTime(time, False)

        count = 0
        if self._sat.getAltitude(self._geo, nextTime) > 0:
            lastPositiveAltitude = nextTime
            count = 1

        while count < abs(number):
            if abs(lastPositiveAltitude - nextTime) > maximumSearchPeriod:
                raise TimeoutError(f'Timeout: No satellite pass found within {maximumSearchPeriod} days.')
            nextTime = self._computeNextMaximumTime(nextTime, nextOccurrence)
            if self._sat.getAltitude(self._geo, nextTime) > 0:
                lastPositiveAltitude = nextTime
                count += 1

        return nextTime

    def computeMaximumTime(self, time: 'JulianDate', number: int, maximumSearchPeriod: float = 30) -> 'JulianDate':
        if number == 0:
            raise ValueError(f'number must be a non-zero integer, not {number}')
        directionSign = signOf(number)
        nextOccurrence = bool(directionSign + 1)

        return self._computeMaximumTimeExec(time, number, nextOccurrence, maximumSearchPeriod)

    def _computeIntersectionAnomalyTerms(self, time: 'JulianDate') -> (float, float):
        """Returns the offset term, then maximum term."""

        a, b, c = computeEllipseVectors(self._sat, time)
        zeta = self._geo.getZenithVector(time)
        gamma = self._geo.getPositionVector(time)

        zDotA = dot(zeta, a)
        zDotB = dot(zeta, b)
        numerator = dot(zeta, gamma - c)
        denominator = sqrt(zDotA * zDotA + zDotB * zDotB)

        try:
            lhs = acos(numerator / denominator)
        except ValueError as e:
            if e.args[0] == 'math domain error':
                raise PassedOrbitPathRange()
            raise e

        rhs = atan3(zDotB, zDotA)
        return lhs, rhs

    def _computeIntersectionTimeExec(self, maxTime: 'JulianDate', isRising: bool) -> 'JulianDate':

        time = maxTime
        prevTime = time.future(-1)

        # epsilon = 10 / SECONDS_PER_DAY
        while abs(prevTime - time) > TIME_DIFFERENCE:
            offsetTerm, maxTerm = self._computeIntersectionAnomalyTerms(time)
            if isRising:
                targetAnomaly = (maxTerm - offsetTerm) % TWOPI
            else:
                targetAnomaly = (maxTerm + offsetTerm) % TWOPI

            prevTime = time
            time = self._sat.timeToNearestAnomaly(targetAnomaly, time, 'true')

        return time

    def _refineIntersectionTime(self, time: 'JulianDate') -> 'JulianDate':
        updatedTime = time
        alt = self._sat.getAltitude(self._geo, updatedTime)
        _, topocentricVelocity = self._sat.getTopocentricState(self._geo, updatedTime)
        direction = -1 if topocentricVelocity[2] > 0 else 1

        while abs(alt) > 4.848e-6:  # 1 arc-second
            # todo: fix how we compute topocentric coordinates (publicly and privately)
            position, velocity = self._sat.getState(updatedTime)

            gamma = self._geo.getPositionVector(updatedTime)
            relativePosition = position - gamma
            topocentricPosition = toTopocentric(relativePosition, self._geo, updatedTime)

            coriolis = cross(EARTH_ANGULAR_MOMENTUM, position)
            relativeVelocity = velocity - coriolis
            topocentricVelocity = toTopocentric(relativeVelocity, self._geo, updatedTime)

            angularVelocityVector = cross(topocentricPosition, topocentricVelocity) / topocentricPosition.mag2()
            dt = (alt / angularVelocityVector.mag()) / SECONDS_PER_DAY

            updatedTime = updatedTime.future(dt * direction)

            alt = self._sat.getAltitude(self._geo, updatedTime)

        return updatedTime

    def computeIntersectionTimes(self, time: 'JulianDate', number: int, maximumSearchPeriod: float = 30)\
            -> ('JulianDate', 'JulianDate'):
        if number > 0:
            maxTime = self._computeMaximumTimeExec(time, number, True, maximumSearchPeriod)
        elif number < 0:
            maxTime = self._computeMaximumTimeExec(time, number, False, maximumSearchPeriod)
        else:
            raise ValueError(f'number must be a non-zero integer, not {number}')

        riseTime = self._computeIntersectionTimeExec(maxTime, True)
        setTime = self._computeIntersectionTimeExec(maxTime, False)

        refinedRiseTime = self._refineIntersectionTime(riseTime)
        refinedSetTime = self._refineIntersectionTime(setTime)

        return refinedRiseTime, refinedSetTime

    def _getNextPassExec(self, maxTime: 'JulianDate') -> SatellitePass:
        riseTime = self._computeIntersectionTimeExec(maxTime, True)
        setTime = self._computeIntersectionTimeExec(maxTime, False)
        riseTime = self._refineIntersectionTime(riseTime)
        setTime = self._refineIntersectionTime(setTime)

        try:
            enterTime, exitTime = getShadowTimes(self._sat, riseTime, Shadow.PENUMBRA)
        except NoSatelliteEclipseException:
            enterTime = setTime.future(0.0001)
            exitTime = enterTime
        # fixme: need to consider latitudes where the sun doesn't rise or set
        sunRiseTime, sunSetTime = Sun.computeRiseSetTimes(self._geo, riseTime)

        topoState = self._sat.getTopocentricState(self._geo, riseTime)
        riseInfo = PositionInfo(0.0, degrees(atan3(topoState[0][1], -topoState[0][0])), riseTime,
                                riseTime < enterTime, riseTime < sunRiseTime or riseTime > sunSetTime)
        topoState = self._sat.getTopocentricState(self._geo, setTime)
        setInfo = PositionInfo(0.0, degrees(atan3(topoState[0][1], -topoState[0][0])), setTime,
                               setTime < enterTime or setTime > exitTime,
                               setTime < sunRiseTime or setTime > sunSetTime)

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

        return SatellitePass(infos, self._sat.name)

    def getNextPass(self, time: 'JulianDate', number: int = 1, maximumSearchPeriod: float = 30) -> SatellitePass:
        if number > 0:
            if isinstance(time, SatellitePass):
                time = time.setInfo.time.future(0.0001)

            findNext = True
        elif number < 0:
            if isinstance(time, SatellitePass):
                time = time.riseInfo.time.future(-0.0001)

            findNext = False
        else:
            raise ValueError(f'number must be a non-zero integer, not {number}')

        approxMaximumTime = self._computeMaximumTimeExec(time, number, findNext, maximumSearchPeriod)
        return self._getNextPassExec(approxMaximumTime)

    def getPassList(self, time: 'JulianDate', duration: float, maximumSearchPeriod: float = 30) -> list[SatellitePass]:
        if duration > 0:
            findNext = True
            firstNumber = 1
            nextNumber = 2
        elif duration < 0:
            findNext = False
            firstNumber = -1
            nextNumber = -2
        else:
            raise ValueError('duration must be a non-zero number of days')

        nextMaximumTime = self._computeMaximumTimeExec(time, firstNumber, findNext, maximumSearchPeriod)

        maximumTimeList = []
        while abs(nextMaximumTime - time) < abs(duration):
            maximumTimeList.append(nextMaximumTime)
            try:
                nextMaximumTime = self._computeMaximumTimeExec(nextMaximumTime, nextNumber, findNext,
                                                               maximumSearchPeriod)
            except (SatelliteAlwaysAbove, NoPassException, PassedOrbitPathRange):
                break

        passList = []
        for passTime in maximumTimeList:
            try:
                nextPass = self._getNextPassExec(passTime)
                passList.append(nextPass)
            except (SatelliteAlwaysAbove, NoPassException, PassedOrbitPathRange):
                break

        return passList


"""
tle = TwoLineElement('''ISS (ZARYA)
1 25544U 98067A   23355.71058934  .00017521  00000+0  30729-3 0  9994
2 25544  51.6406 119.6587 0002362  41.7601  99.6210 15.50640663430878''')
geo = GeoPosition(38.0608, -97.9298, 0.0)
jd = JulianDate(2023, 12, 23, 18, 53, 44.938, -6.0)
"""