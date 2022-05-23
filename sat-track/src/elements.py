from juliandate import JulianDate
from constants import EARTH_MU
from anomalies import trueToMean, meanToTrue
from math import sqrt, radians, degrees


class OrbitalElements:

    def __init__(self, sma: float, ecc: float, inc: float, raan: float, aop: float, trueAnomaly: float,
                 epoch: JulianDate = 0):
        self._sma = sma
        self._ecc = ecc
        self._inc = inc
        self._raan = raan
        self._aop = aop
        self._trueAnomaly = trueAnomaly
        self._epoch = epoch

    def __str__(self) -> str:
        return (f'Semi-major axis: {self._sma}, Eccentricity: {self._ecc}, Inclination: {self._inc}\n'
                + f'RAAN: {self._raan}, Argument of Perigee: {self._aop}, True Anomaly: {self._trueAnomaly}\n'
                + f'Epoch: {str(self._epoch)}')

    def moveTo(self, jd: JulianDate):
        """Moves a satellites angular position to a specified time using the mean motion of the object. Time is
        determined by a future or past Julian Date.s
        Parameters:
        jd: The Julian Date you wish to move the object to."""
        self.moveBy(jd.difference(self._epoch))

    def moveBy(self, dt: float):
        """Moves a satellites angular position by a number of solar days, using the mean motion of the object.
        A negative value indicates a past position.
        Parameters:
        dt: The number of solar days to move the object to."""
        dtp = 86400 * dt
        M1rad = radians(trueToMean(self._trueAnomaly, self._ecc)) + dtp * sqrt(EARTH_MU / (self._sma ** 3))
        self._trueAnomaly = meanToTrue(degrees(M1rad), self._ecc)
        self._epoch = self._epoch.future(dt)

    def timeUntil(self, trueAnomaly: float) -> float:
        """Returns the time until the object achieves a specified true anomaly.
            If the trueAnomaly parameter is equal to the current true anomaly, the returned value will be 0,
            NOT the orbital period. I.e. the angle needed to travel will be 0, not 360 degrees.
        Parameters:
        trueAnomaly:    The time computed will be the time it takes to travel to this true anomaly."""
        M0 = trueToMean(self._trueAnomaly, self._ecc)
        M1 = trueToMean(trueAnomaly, self._ecc)
        da = M1 - M0 if M1 >= M0 else M1 + 360.0 - M0
        return radians(da) / sqrt(EARTH_MU / (self._sma ** 3))

    def setSma(self, sma: float):
        self._sma = sma

    def getSma(self) -> float:
        return self._sma

    def setEcc(self, ecc: float):
        self._ecc = ecc

    def getEcc(self):
        return self._ecc

    def setInc(self, inc: float):
        self._inc = inc

    def getInc(self):
        return self._inc

    def setRaan(self, raan: float):
        self._raan = raan

    def getRaan(self):
        return self._raan

    def setAop(self, aop: float):
        self._aop = aop

    def getAop(self):
        return self._aop

    def setTrueAnomaly(self, trueAnom: float):
        self._trueAnomaly = trueAnom

    def getTrueAnomaly(self):
        return self._trueAnomaly

    def setEpoch(self, epoch: JulianDate):
        self._epoch = epoch

    def getEpoch(self):
        return self._epoch
