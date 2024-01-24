import json
from math import asin, degrees, sqrt, acos, cos, sin

from pyevspace import dot, cross, getMatrixEuler, ZXZ, Angles, rotateMatrixFrom, Vector

from sattrack.bodies.position import computeEarthOffsetAngle
from sattrack.bodies.sun import Sun
from sattrack.bodies.topocentric import toTopocentric
from sattrack.core._analysis import AnalyticFunction, Boundary
from sattrack.orbit.elements import smaToMeanMotion, eccentricToMeanAnomaly
from sattrack.satellitepass.eclipse import getShadowTimes, Shadow
from sattrack.orbit.exceptions import SatelliteAlwaysAbove, NoPassException, PassedOrbitPathRange
from sattrack.satellitepass.exceptions import NoSatelliteEclipseException
from sattrack.orbit.orbitpath import OrbitPath, computeEllipseVectors
from sattrack.config import TIME_DIFFERENCE, MAXIMUM_ANOMALY_DISCONTINUITY_GAP, MAXIMUM_ANOMALY_DIFFERENCE, \
    MAXIMUM_ANOMALY_STEP_SIZE, MAXIMUM_ANOMALY_CHUNK_DURATION, MAXIMUM_ANOMALY_REFINE_GAP
from sattrack.satellitepass.info import PositionInfo, Visibility
from sattrack.util.constants import SECONDS_PER_DAY, MINUTES_PER_DAY, TWOPI, EARTH_ANGULAR_MOMENTUM, \
    EARTH_ANGULAR_VELOCITY
from sattrack.util.helpers import atan3, signOf, computeAngleDifference

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.orbit.satellite import Orbitable
    from sattrack.core.coordinates import GeoPosition
    from sattrack.core.juliandate import JulianDate
    from sattrack.orbit.elements import Elements


class ZenithVectorComputer:
    __slots__ = '_referenceTime', '_sinLatitude', '_cosLatitude', '_localSiderealTime',\
                '_EARTH_ANG_VEL'

    def __init__(self, geo: 'GeoPosition', time: 'JulianDate'):
        self._referenceTime = time
        self._sinLatitude = sin(geo.latitudeRadians)
        self._cosLatitude = cos(geo.latitudeRadians)

        earthOffsetAngle = computeEarthOffsetAngle(time)
        self._localSiderealTime = (geo.longitudeRadians + earthOffsetAngle) % TWOPI
        self._EARTH_ANG_VEL = EARTH_ANGULAR_VELOCITY * TWOPI

    def computeZenithVector(self, time: 'JulianDate') -> Vector:
        dt = time - self._referenceTime
        longitudeParameter = self._localSiderealTime + self._EARTH_ANG_VEL * dt

        xTerm = self._cosLatitude * cos(longitudeParameter)
        yTerm = self._cosLatitude * sin(longitudeParameter)

        return Vector(xTerm, yTerm, self._sinLatitude)


class MaximumAnomalyFunction(AnalyticFunction):
    __slots__ = ('_sat', '_geo', '_a', '_b', '_c', '_eccentricity', '_meanAnomaly',
                 '_meanMotion', '_zenithComputer')

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition', time: 'JulianDate'):
        self._sat = satellite
        self._geo = geo

        self._zenithComputer = ZenithVectorComputer(self._geo, time)
        self._updateConstants(time)

    def _updateConstants(self, time: 'JulianDate'):
        elements = self._sat.getElements(time)

        self._a, self._b, self._c = self._computeEllipseVectorsFast(elements)
        self._eccentricity = elements.ecc
        self._meanAnomaly = elements.meanAnomaly

        meanMotion = smaToMeanMotion(elements.sma, self._sat.body.mu)
        self._meanMotion = meanMotion * SECONDS_PER_DAY

    # todo: since this is used elsewhere, it probably should be make a global method
    @staticmethod
    def _computeEllipseVectorsFast(elements: 'Elements') -> (Vector, Vector, Vector):
        matrix = getMatrixEuler(ZXZ, Angles(elements.raan, elements.inc, elements.aop))
        aMag = elements.sma
        bMag = aMag * sqrt(1 - elements.ecc * elements.ecc)
        cMag = aMag * elements.ecc

        a = rotateMatrixFrom(matrix, Vector(aMag, 0, 0))
        b = rotateMatrixFrom(matrix, Vector(0, bMag, 0))
        c = rotateMatrixFrom(matrix, Vector(-cMag, 0, 0))

        return a, b, c

    # def _computeZenithVectorFast(self, time: 'JulianDate') -> Vector:
    #     dt = time - self._referenceTime
    #     longitudeParameter = self._localSiderealTime + self._w * dt
    #
    #     xTerm = self._cosLatitude * cos(longitudeParameter)
    #     yTerm = self._cosLatitude * sin(longitudeParameter)
    #
    #     return Vector(xTerm, yTerm, self._sinLatitude)

    def _computeMeanAnomalyAtMaximum(self, time: 'JulianDate') -> float:
        # zenithVector = self._computeZenithVectorFast(time)
        zenithVector = self._zenithComputer.computeZenithVector(time)
        zDotA = dot(zenithVector, self._a)
        zDotB = dot(zenithVector, self._b)

        eccentricAnomaly = atan3(zDotB, zDotA)
        return eccentricToMeanAnomaly(eccentricAnomaly, self._eccentricity)

    def compute(self, time: 'JulianDate') -> float:
        self._updateConstants(time)

        maximumMeanAnomaly = self._computeMeanAnomalyAtMaximum(time)
        rawAngleDifference = maximumMeanAnomaly - self._meanAnomaly

        return computeAngleDifference(rawAngleDifference)


class MaximumMeanAnomalyFinder:
    __slots__ = '_sat', '_geo', '_referenceTime', '_function', '_nextOccurrence'

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition'):
        self._sat = satellite
        self._geo = geo
        self._referenceTime: 'JulianDate' = None
        self._function: AnalyticFunction = None
        self._nextOccurrence: bool = None

    @staticmethod
    def _refineApproximateTime(boundary: Boundary, epsilon: float) -> Boundary:
        while (boundaryRange := boundary.range()) > epsilon:
            middleTime = boundary.lower.x.future(boundaryRange / 2)
            left, right = boundary.bifurcate(middleTime)
            boundary = left if right.hasSameSign() else right

        return boundary

    def _checkPositiveAltitude(self, time: 'JulianDate') -> bool:
        """Returns True if the Orbitable is above the GeoPosition at time."""

        return self._sat.getAltitude(self._geo, time) > 0

    def _toContinueSearch(self, currentTime: 'JulianDate', endTime: 'JulianDate',
                          duration: float) -> bool:
        """Determines if the search algorithm should continue searching or if the
        current chunk has completed its search."""

        if self._nextOccurrence:
            return endTime - currentTime < duration
        else:
            return currentTime - endTime < duration

    def _searchChunk(self, time: 'JulianDate', step: float, duration: float) -> Boundary | None:
        """Searches a chunk of time incrementally, looking for roots to the instances
        AnalyticFunction. The step parameter defines how much the x-value changes in
        computing the next Point of the AnalyticFunction, positive for the next occurrence,
        negative when looking for the previous. The duration parameter defines how long the
        chunk for the current search period is in solar days from time."""

        basePoint = self._function.computePoint(time)
        nextPoint = basePoint
        while self._toContinueSearch(nextPoint.x, time, duration):
            basePoint = nextPoint
            nextTime = basePoint.x.future(step)
            nextPoint = self._function.computePoint(nextTime)

            if basePoint.y * nextPoint.y <= 0:
                initialBoundary = Boundary(basePoint, nextPoint, self._function)
                boundary = self._refineApproximateTime(initialBoundary, MAXIMUM_ANOMALY_DISCONTINUITY_GAP)

                if boundary.difference() < MAXIMUM_ANOMALY_DIFFERENCE:
                    return boundary

    def _execFindMaximumTimes(self, time: 'JulianDate', nextOccurrence: bool, duration: float,
                              onlyOnce: bool) -> list['JulianDate']:
        """Returns the list of times corresponding to each visible maximum anomaly achieved.
        The time space searched is defined by duration in solar days if onlyOnce is False. If
        onlyOnce is True, duration is the timeout period when the function will return if no
        pass is found within the duration."""

        self._function = MaximumAnomalyFunction(self._sat, self._geo, time)
        self._nextOccurrence = nextOccurrence

        chunkStartTime = time
        if nextOccurrence:
            step = MAXIMUM_ANOMALY_STEP_SIZE
            deltaTime = 0.001
        else:
            step = -MAXIMUM_ANOMALY_STEP_SIZE
            deltaTime = -0.001

        maximumAnomalyList = []
        while abs(chunkStartTime - time) < duration:
            chunkDuration = min((time.future(duration) - chunkStartTime, MAXIMUM_ANOMALY_CHUNK_DURATION))

            result = self._searchChunk(chunkStartTime, step, chunkDuration)
            if result is not None:
                boundary = self._refineApproximateTime(result, MAXIMUM_ANOMALY_REFINE_GAP)
                exactRoot = boundary.lower.x

                if abs(exactRoot - time) > duration:
                    break

                if self._checkPositiveAltitude(exactRoot):
                    if onlyOnce:
                        return [exactRoot]
                    else:
                        maximumAnomalyList.append(exactRoot)

                chunkStartTime = exactRoot.future(deltaTime)
            else:
                chunkStartTime = chunkStartTime.future(step)

        return maximumAnomalyList

    def computeMaximumTime(self, time: 'JulianDate', nextOccurrence: bool = True,
                           timeout: float = 14) -> list['JulianDate']:
        """Compute the next or previous time the satellite achieves the orbit's maximum
        mean anomaly and the satellite's altitude is positive at this time. If nextOccurrence
        is True, the next maximum is found, otherwise the previous is found. The timeout
        parameter is used to determine when to stop looking for a maximum time if none
        are found with the time period.

        This method returns a list with either 0 or 1 JulianDate objects. This is to ensure
        consistency with the class's computeMaximumTimeList method."""

        return self._execFindMaximumTimes(time, nextOccurrence, timeout, True)

    def computeMaximumTimeList(self, time: 'JulianDate', nextOccurrence: bool = True,
                               duration: float = 7) -> list['JulianDate']:
        """Compute all following or preceding times the satellite achieves the orbit's
        maximum mean anomaly and the corresponding altitude is positive. If nextOccurrence
        is True the following maximums are found, otherwise the preceding are found. The
        duration parameter expresses the time period to find all positive altitude maximums."""

        return self._execFindMaximumTimes(time, nextOccurrence, duration, False)


class MeanAnomalyIntersectionFinder:
    __slots__ = '_sat', '_geo', '_zenithComputer', '_a', '_b', '_c', '_eccentricity'

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition'):
        self._sat = satellite
        self._geo = geo
        self._zenithComputer: ZenithVectorComputer = None

    def _computeIntersectionAnomalyTerms(self, time: 'JulianDate') -> (float, float):
        elements = self._sat.getElements(time)
        a, b, c = MaximumAnomalyFunction._computeEllipseVectorsFast(elements)
        zeta = self._zenithComputer.computeZenithVector(time)
        # todo: make a fast way to compute this
        gamma = self._geo.getPositionVector(time)

        zDotA = dot(zeta, a)
        zDotB = dot(zeta, b)
        numerator = dot(zeta, gamma - c)
        denominator = sqrt(zDotA * zDotA + zDotB * zDotB)

        # fixme: if this happens, we could catch this and instead of calling this with time,
        #   call with a time that's a fraction of the way from maximum time to time.
        try:
            offsetAdjustment = cos(numerator / denominator)
        except ValueError as e:
            if e.args[0] == 'math domain error':
                raise PassedOrbitPathRange()
            raise e

        middleMeanAnomaly = atan3(zDotB, zDotA)
        return middleMeanAnomaly, offsetAdjustment

    def _execComputeIntersectionTime(self, passTime: 'JulianDate', isRising: bool) -> 'JulianDate':
        time = passTime
        previousTime = time.future(-1)

        while abs(previousTime - time) > TIME_DIFFERENCE:
            middleMeanAnomaly, offsetAdjustment = self._computeIntersectionAnomalyTerms(time)
            if isRising:
                targetMeanAnomaly = (middleMeanAnomaly - offsetAdjustment) % TWOPI
            else:
                targetMeanAnomaly = (middleMeanAnomaly + offsetAdjustment) % TWOPI

            previousTime = time
            time = self._sat.timeToNearestAnomaly(targetMeanAnomaly, time, 'mean')

        return time

    def _refineIntersectionTime(self, time: 'JulianDate') -> 'JulianDate':
        # todo: use the altitude derivative to compute this more accurately (i hope)
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

    def computeIntersectionTimes(self, time: 'JulianDate') -> ('JulianDate', 'JulianDate'):
        if self._zenithComputer is None:
            self._zenithComputer = ZenithVectorComputer(self._geo, time)

        riseTime = self._execComputeIntersectionTime(time, True)
        setTime = self._execComputeIntersectionTime(time, False)

        refinedRiseTime = self._refineIntersectionTime(riseTime)
        refinedSetTime = self._refineIntersectionTime(setTime)

        return refinedRiseTime, refinedSetTime


class PassTimeController:
    __slots__ = '_sat', '_geo', '_maximumAnomalyFinder', '_anomalyIntersectionFinder'

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition'):
        self._sat = satellite
        self._geo = geo
        self._maximumAnomalyFinder = MaximumMeanAnomalyFinder(satellite, geo)
        self._anomalyIntersectionFinder = MeanAnomalyIntersectionFinder(self._sat, self._geo)

    def _isOrbitAbove(self, time: 'JulianDate') -> bool:
        a, b, c = computeEllipseVectors(self._sat, time)
        gamma = self._geo.getPositionVector(time)
        zeta = self._geo.getZenithVector(time)

        numerator = dot(zeta, gamma - c)
        zDotA = dot(zeta, a)
        zDotB = dot(zeta, b)
        denominator = sqrt(zDotA*zDotA + zDotB*zDotB)

        try:
            acos(numerator / denominator)
        except ValueError as e:
            if e.args[0] == 'math domain error':
                return False
            raise e
        else:
            return True

    def _determineException(self, time: 'JulianDate'):
        satelliteName = self._sat.name or 'Satellite'

        # This will be more accurate descriptively if we know whether or not the
        # orbit rises and sets or not.
        if self._sat.isGeosynchronous():
            if self._isOrbitAbove(time):
                if self._sat.getAltitude(self._geo, time) > 0:
                    raise SatelliteAlwaysAbove(f'{satelliteName} is geo-synchronous and always above the geo-position')
                raise NoPassException(f'{satelliteName} is geo-synchronous and not above the geo-position')
            else:
                raise NoPassException(f'{satelliteName} is geo-synchronous whose orbit is not above the geo-position')

        raise NoPassException(f'{satelliteName} does not pass over the geo-position')

    def _refineToTrueMaximum(self, time: 'JulianDate') -> 'JulianDate':
        # todo: use the derivative of altitude to compute the true maximum
        return time

    def computePassTime(self, time: 'JulianDate', nextOccurrence: bool = True,
                        timeout: float = 14) -> 'JulianDate':
        result = self._maximumAnomalyFinder.computeMaximumTime(time, nextOccurrence, timeout)

        if not result:
            self._determineException(time)

        return self._refineToTrueMaximum(result[0])

    def computePassList(self, time: 'JulianDate', nextOccurrence: bool = True,
                        duration: float = 7) -> 'list[JulianDate]':
        results = self._maximumAnomalyFinder.computeMaximumTimeList(time, nextOccurrence, duration)

        if not results:
            self._determineException(time)

        preciseResults = [self._refineToTrueMaximum(time) for time in results]
        return preciseResults

    def computeIntersectionTimes(self, time: 'JulianDate') -> 'JulianDate, JulianDate':
        return self._anomalyIntersectionFinder.computeIntersectionTimes(time)


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

    def __getattr__(self, name: str):
        if name.endswith('Info'):
            nameKey = name[:-4]
            infos = self._infos
            if nameKey in infos:
                return infos.get(nameKey)

        raise AttributeError()


class PassFinder:
    __slots__ = '_sat', '_geo', '_passController'

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition'):
        self._sat = satellite
        self._geo = geo
        self._passController = PassTimeController(self._sat, self._geo)

    @property
    def orbitable(self) -> 'Orbitable':
        return self._sat

    @property
    def geo(self) -> 'GeoPosition':
        return self._geo

    @staticmethod
    def _computeAzimuth(topocentricState: list[Vector]) -> float:
        azimuthRadians = atan3(topocentricState[0][1], -topocentricState[0][0])

        return degrees(azimuthRadians)

    @staticmethod
    def _computeAltitude(topocentricState: list[Vector]) -> float:
        altitudeRadians = asin(topocentricState[0][2] / topocentricState[0].mag())

        return degrees(altitudeRadians)

    def _generatePositionInfo(self, time: 'JulianDate', illuminated: bool, unobscured: bool) -> PositionInfo:
        topocentricState = self._sat.getTopocentricState(self._geo, time)
        altitude = self._computeAltitude(topocentricState)
        azimuth = self._computeAzimuth(topocentricState)

        return PositionInfo(altitude, azimuth, time, illuminated, unobscured)

    def _computeSatellitePass(self, maximumTime: 'JulianDate') -> SatellitePass:
        riseTime, setTime = self._passController.computeIntersectionTimes(maximumTime)

        try:
            shadowEnterTime, shadowExitTime = getShadowTimes(self._sat, riseTime, Shadow.PENUMBRA)
        except NoSatelliteEclipseException:
            shadowEnterTime = setTime.future(0.0001)
            shadowExitTime = shadowEnterTime
        # fixme: need to consider latitudes where the sun doesn't rise or set
        sunRiseTime, sunSetTime = Sun.computeRiseSetTimes(self._geo, riseTime)

        topocentricState = self._sat.getTopocentricState(self._geo, riseTime)
        azimuth = self._computeAzimuth(topocentricState)
        riseInfo = PositionInfo(0.0, azimuth, riseTime, riseTime < shadowEnterTime,
                                riseTime < sunRiseTime or riseTime > sunSetTime)
        topocentricState = self._sat.getTopocentricState(self._geo, setTime)
        azimuth = self._computeAzimuth(topocentricState)
        setInfo = PositionInfo(0.0, azimuth, setTime, setTime < shadowEnterTime or setTime > shadowExitTime,
                               setTime < sunRiseTime or setTime > sunSetTime)
        topocentricState = self._sat.getTopocentricState(self._geo, maximumTime)
        altitude = self._computeAltitude(topocentricState)
        azimuth = self._computeAzimuth(topocentricState)
        maxInfo = PositionInfo(altitude, azimuth, maximumTime,
                               maximumTime < shadowEnterTime or maximumTime > shadowExitTime,
                               maximumTime < sunRiseTime or maximumTime > sunSetTime)

        firstIlluminatedInfo = lastIlluminatedInfo = firstUnobscuredInfo = lastUnobscuredInfo = None
        if riseTime < shadowExitTime < setTime:
            firstIlluminatedInfo = self._generatePositionInfo(shadowExitTime, True,
                                                              shadowExitTime < sunRiseTime or
                                                              shadowExitTime >= sunSetTime)
        if riseTime < shadowEnterTime < setTime:
            lastIlluminatedInfo = self._generatePositionInfo(shadowEnterTime, True,
                                                             shadowEnterTime < sunRiseTime or
                                                             shadowEnterTime >= sunSetTime)
        if riseTime < sunSetTime < setTime:
            firstUnobscuredInfo = self._generatePositionInfo(sunSetTime, sunSetTime < shadowEnterTime or
                                                             sunSetTime >= shadowExitTime, True)
        if riseTime < sunRiseTime < setTime:
            lastUnobscuredInfo = self._generatePositionInfo(sunRiseTime, sunRiseTime < shadowEnterTime or
                                                            sunRiseTime >= shadowExitTime, True)

        alternateInfos = (firstIlluminatedInfo, lastIlluminatedInfo, firstUnobscuredInfo, lastUnobscuredInfo)
        infos = [riseInfo, setInfo, maxInfo] + [info for info in alternateInfos if info is not None]

        return SatellitePass(infos, self._sat.name)

    def computeNextPass(self, time: 'JulianDate', nextOccurrence: bool = True, timeout: float = 7) -> SatellitePass:
        if isinstance(time, SatellitePass):
            if nextOccurrence:
                time = time.setInfo.time.future(0.001)
            else:
                time = time.riseInfo.time.future(-0.001)

        maximumTime = self._passController.computePassTime(time, nextOccurrence, timeout)

        return self._computeSatellitePass(maximumTime)

    def computePassList(self, time: 'JulianDate', duration: float) -> list[SatellitePass]:
        if duration > 0:
            nextOccurrence = True
        elif duration < 0:
            nextOccurrence = False
        else:
            raise ValueError('duration must be a non-zero number of days')

        maximumTimeList = self._passController.computePassList(time, nextOccurrence, duration)
        satellitePassList = [self._computeSatellitePass(time) for time in maximumTimeList]

        return satellitePassList


class PassController:
    __slots__ = '_sat', '_geo', '_path', '_finder'

    def __init__(self, satellite: 'Orbitable', geo: 'GeoPosition'):
        self._sat = satellite
        self._geo = geo
        self._path = OrbitPath(satellite, geo)
        self._finder = PassFinder(satellite, geo)

    @property
    def orbitable(self):
        return self._sat

    @property
    def geo(self):
        return self._geo

    @property
    def path(self):
        return self._path

    def getNextPass(self, time: 'JulianDate', number: int = 1, maximumSearchPeriod: float = 14) -> SatellitePass:
        if number > 0:
            nextOccurrence = True
        elif number < 0:
            nextOccurrence = False
        else:
            raise ValueError(f'number must be a non-zero integer, not {number}')

        return self._finder.computeNextPass(time, nextOccurrence, maximumSearchPeriod)

    def getPassList(self, time: 'JulianDate', duration: float, UNUSED: float = 0) -> list[SatellitePass]:
        if duration == 0:
            raise ValueError(f'duration must be a non-zero number of days')

        return self._finder.computePassList(time, duration)


class PassControllerOld:
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
        # This range is way too big. This may be called right before a pass occurs, so the nextMaxTime happens
        # within this range, and 'looks' like an error, but is valid.
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

"""
tle = TwoLineElement('''STARLINK-G7-9 SINGLE
1 72001C 24002B   24003.19929884  .00290917  00000 0  17947-2 0    07
2 72001  53.1623 248.1741 0006574 256.7333 109.8990 15.76335933    12''')
geo = GeoPosition(38.0826774, -97.9036822)
jd = JulianDate(2024, 1, 3, 21, 12, 13.60554800000682, -6.0)
"""