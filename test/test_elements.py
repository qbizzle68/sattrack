import unittest
from math import pi, degrees

from sattrack.orbit.sgp4 import TwoLineElement, elementsFromState

from sattrack.core.juliandate import JulianDate
from sattrack.orbit.elements import elementsFromTle, Elements, radiusAtPeriapsis, radiusAtApoapsis, trueToMeanAnomaly, \
    trueToEccentricAnomaly, meanToTrueAnomaly, meanToEccentricAnomaly, eccentricToTrueAnomaly, eccentricToMeanAnomaly, \
    radiusAtAnomaly, flightAngleAtAnomaly, velocityAtAnomaly, nextMeanAnomaly, previousMeanAnomaly, nextTrueAnomaly, \
    previousTrueAnomaly, nearestTrueAnomaly, nearestMeanAnomaly, meanAnomalyAtTime, trueAnomalyAtTime, computeAnomaly, \
    smaToMeanMotion, meanMotionToSma
from sattrack.orbit.satellite import Satellite
from sattrack.util.constants import EARTH_EQUITORIAL_RADIUS, EARTH_MU


class TestElements(unittest.TestCase):

    tle = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23301.52764947  .00019202  00000+0  34666-3 0  9994
2 25544  51.6435  28.0714 0001006 114.0425 254.5243 15.49747991422472''')
        cls.sat = Satellite(cls.tle)
        cls.jd = JulianDate(2023, 10, 28, 11, 0, 0, -5)

    def testFromTLE(self):
        raan, inc, aop, ecc, sma, meanAnomaly = elementsFromTle(self.tle, self.jd)

        self.assertAlmostEqual(0.4778769422144868, raan)
        self.assertAlmostEqual(0.9013491122536916, inc)
        self.assertAlmostEqual(1.999410681682025, aop)
        self.assertAlmostEqual(9.830359251717191e-05, ecc)
        self.assertAlmostEqual(6788.48667486485, sma)
        self.assertAlmostEqual(5.41253602664068, meanAnomaly)

        elements = Elements.fromTle(self.tle, self.jd)
        self.assertAlmostEqual(0.4778769422144868, elements.raan)
        self.assertAlmostEqual(0.9013491122536916, elements.inc)
        self.assertAlmostEqual(1.999410681682025, elements.aop)
        self.assertAlmostEqual(9.830359251717191e-05, elements.ecc)
        self.assertAlmostEqual(6788.48667486485, elements.sma)
        self.assertAlmostEqual(5.41253602664068, elements.meanAnomaly)

    def testElements(self):
        elements = Elements(pi / 2, pi / 4, pi, 0, EARTH_EQUITORIAL_RADIUS + 400, pi, self.jd)
        elementsDegrees = Elements.fromDegrees(90, 45, 180, 0, EARTH_EQUITORIAL_RADIUS + 400, 180, self.jd)

        self.assertAlmostEqual(elements.raan, elementsDegrees.raan)
        self.assertAlmostEqual(elements.inc, elementsDegrees.inc)
        self.assertAlmostEqual(elements.aop, elementsDegrees.aop)
        self.assertAlmostEqual(elements.ecc, elementsDegrees.ecc)
        self.assertAlmostEqual(elements.sma, elementsDegrees.sma)
        self.assertAlmostEqual(elements.meanAnomaly, elementsDegrees.meanAnomaly)
        self.assertAlmostEqual(elements.trueAnomaly, elementsDegrees.trueAnomaly)

        elementsDegrees = Elements.fromDegrees(90, 45, 180, 0, EARTH_EQUITORIAL_RADIUS + 400, 180, self.jd,
                                               degrees(elements.trueAnomaly))
        self.assertAlmostEqual(elementsDegrees.trueAnomaly, elements.trueAnomaly)

    def testFromState(self):
        state = self.sat.getState(self.jd)
        elements = Elements.fromState(*state, self.jd)

        self.assertAlmostEqual(elements.raan, 0.47826202632911824)
        self.assertAlmostEqual(elements.inc, 0.9011289790793156)
        self.assertAlmostEqual(elements.aop, 2.067480647498184)
        self.assertAlmostEqual(elements.ecc, 0.000834494512920963)
        self.assertAlmostEqual(elements.sma, 6792.313051673292)
        self.assertAlmostEqual(elements.meanAnomaly, 5.344712934053373)
        self.assertAlmostEqual(elements.trueAnomaly, 5.34336580365274)

        raan, inc, aop, ecc, sma, meanAnomaly, trueAnomaly = elementsFromState(*state, EARTH_MU)

        self.assertAlmostEqual(raan, 0.47826202632911824)
        self.assertAlmostEqual(inc, 0.9011289790793156)
        self.assertAlmostEqual(aop, 2.067480647498184)
        self.assertAlmostEqual(ecc, 0.000834494512920963)
        self.assertAlmostEqual(sma, 6792.313051673292)
        self.assertAlmostEqual(meanAnomaly, 5.344712934053373)
        self.assertAlmostEqual(elements.trueAnomaly, 5.34336580365274)

    def testUpdateAnomaly(self):
        elements = Elements(pi / 2, pi / 4, pi, 0, EARTH_EQUITORIAL_RADIUS + 400, pi, self.jd)
        anomUpdate = 0.0
        jd = self.jd.future(0.5)

        # True anomaly update.
        elements.updateAnomaly('true', anomUpdate, jd)
        self.assertAlmostEqual(elements.trueAnomaly, anomUpdate)
        self.assertAlmostEqual(elements.meanAnomaly, anomUpdate)
        self.assertAlmostEqual(elements.epoch, jd)

        # Mean anomaly update.
        anomUpdate = pi
        jd = self.jd.future(1.0)
        elements.updateAnomaly('mean', anomUpdate, jd)
        self.assertAlmostEqual(elements.trueAnomaly, anomUpdate)
        self.assertAlmostEqual(elements.meanAnomaly, anomUpdate)
        self.assertAlmostEqual(elements.epoch, jd)

        with self.assertRaises(ValueError):
            elements.updateAnomaly('foo', 0.0, self.jd)

    def testString(self):
        elements = Elements.fromTle(self.tle, self.jd)
        ans = ' elements |  raan   |   inc   |   aop   |   ecc    |   sma    | mean anom | true anom |' \
              '        epoch         \n  values  | 27.3803 | 51.6435 | 114.558 | 0.000098 | 6788.487 |' \
              '  310.115  |  310.107  | 2023/10/28 11:00:00.0 -5 UTC'
        self.assertEqual(str(elements), ans)

    def testRadiusMethods(self):
        periapsis = radiusAtPeriapsis(6700, 0.1)
        self.assertAlmostEqual(periapsis, 6030)

        apoapsis = radiusAtApoapsis(6700, 0.1)
        self.assertAlmostEqual(apoapsis, 7370)

    def testAnomalies(self):
        # All anomalies should be equal.
        self.assertAlmostEqual(trueToMeanAnomaly(0.0, 0.1), 0.0)
        self.assertAlmostEqual(trueToEccentricAnomaly(0.0, 0.0), 0.0)
        self.assertAlmostEqual(meanToTrueAnomaly(0.0, 0.99), 0.0)
        self.assertAlmostEqual(meanToEccentricAnomaly(0.0, 0.5), 0.0)
        self.assertAlmostEqual(eccentricToTrueAnomaly(0.0, 0.25), 0.0)
        self.assertAlmostEqual(eccentricToMeanAnomaly(0.0, 0.0), 0.0)

        self.assertAlmostEqual(trueToMeanAnomaly(pi / 2, 0.1), 1.3711301619226748)
        self.assertAlmostEqual(trueToEccentricAnomaly(pi / 2, 0.1), 1.4706289056333368)
        self.assertAlmostEqual(meanToTrueAnomaly(pi / 2, 0.1), 1.7707963267948965)
        # This tests when eccentricity > MAXIMUM_ECCENTRICITY, which forces the Newton method.
        self.assertAlmostEqual(meanToTrueAnomaly(pi / 2, 0.2), 1.9606920626749207)
        self.assertAlmostEqual(meanToEccentricAnomaly(pi / 2, 0.1), 1.6703016694822845)
        self.assertAlmostEqual(eccentricToTrueAnomaly(pi / 2, 0.1), 1.6709637479564563)
        self.assertAlmostEqual(eccentricToMeanAnomaly(pi / 2, 0.1), 1.4707963267948965)

    def testValuesAtAnomalies(self):
        self.assertAlmostEqual(radiusAtAnomaly(6700, 0.0, 1.2345), 6700)
        self.assertAlmostEqual(radiusAtAnomaly(6700, 0.01, pi / 4), 6652.291197835827)

        # Compare with other methods.
        self.assertAlmostEqual(radiusAtPeriapsis(6700, 0.1), radiusAtAnomaly(6700, 0.1, 0))
        self.assertAlmostEqual(radiusAtApoapsis(6700, 0.1), radiusAtAnomaly(6700, 0.1, pi))

        self.assertAlmostEqual(flightAngleAtAnomaly(0.1, 0), 0)
        self.assertAlmostEqual(flightAngleAtAnomaly(0.1, pi), 0)
        self.assertAlmostEqual(flightAngleAtAnomaly(0.1, pi / 4), 0.0659451227998346)
        self.assertAlmostEqual(flightAngleAtAnomaly(0.5, pi / 2), -flightAngleAtAnomaly(0.5, 3 * pi / 2))

        self.assertAlmostEqual(velocityAtAnomaly(6700, 6652.291197835827, EARTH_MU), 7.768264900982522)

    def testNextPreviousAnomalies(self):
        ans = JulianDate(2023, 10, 28, 11, 58, 3.8709309697151184, -5)
        self.assertAlmostEqual(nextMeanAnomaly(15.5, pi / 4, self.jd, 3 * pi / 2, self.jd), ans)
        ans = JulianDate(2023, 10, 28, 11, 58, 3.870971202850342, -5)
        self.assertAlmostEqual(nextMeanAnomaly(15.5, pi, self.jd, pi / 4, self.jd), ans)
        ans = JulianDate(2023, 10, 28, 10, 25, 9.677387773990631, -5)
        self.assertAlmostEqual(previousMeanAnomaly(15.5, pi / 4, self.jd, 3 * pi / 2, self.jd), ans)

        ans = JulianDate(2023, 10, 28, 12, 3, 0.01492917537689209, -5)
        self.assertAlmostEqual(nextTrueAnomaly(15.5, 0.1, pi / 4, self.jd, 3 * pi / 2, self.jd), ans)
        ans = JulianDate(2023, 10, 28, 10, 30, 5.821385979652405, -5)
        self.assertAlmostEqual(previousTrueAnomaly(15.5, 0.1, pi / 4, self.jd, 3 * pi / 2, self.jd), ans)

    def testNearestAnomaly(self):
        ans = JulianDate(2023, 10, 28, 10, 30, 5.821385979652405, -5)
        self.assertAlmostEqual(nearestTrueAnomaly(15.5, 0.1, pi / 4, self.jd, 3 * pi / 2), ans)
        ans = JulianDate(2023, 10, 28, 10, 25, 9.677387773990631, -5)
        self.assertAlmostEqual(nearestMeanAnomaly(15.5, pi / 4, self.jd, 3 * pi / 2), ans)
        ans = JulianDate(2023, 10, 28, 11, 23, 13.548372387886047, -5)
        self.assertAlmostEqual(nearestMeanAnomaly(15.5, pi / 4, self.jd, 3 * pi / 4), ans)
        ans = JulianDate(2023, 10, 28, 10, 25, 9.677387773990631, -5)
        self.assertAlmostEqual(nearestMeanAnomaly(15.5, 3 * pi / 2, self.jd, 3 * pi / 4), ans)
        ans = JulianDate(2023, 10, 28, 11, 34, 50.32258540391922, -5)
        self.assertAlmostEqual(nearestMeanAnomaly(15.5, 3 * pi / 2, self.jd, pi / 4), ans)

    def testComputeAnomaly(self):
        ans = meanAnomalyAtTime(15.5, pi / 4, self.jd, self.jd.future(0.25))
        self.assertAlmostEqual(ans, 1.2563050991728986)

        ans = trueAnomalyAtTime(15.5, 0.1, pi / 4, self.jd, self.jd.future(0.25))
        self.assertAlmostEqual(ans, 1.3121384540782834)

        state = self.sat.getState(self.jd)
        ans = computeAnomaly('true', *state)
        self.assertAlmostEqual(ans, 5.343365803652744)
        ans = computeAnomaly('mean', *state)
        self.assertAlmostEqual(ans, 5.344712934053373)

        with self.assertRaises(ValueError):
            computeAnomaly('foo', *state)

    def testSma(self):
        self.assertAlmostEqual(smaToMeanMotion(6700, EARTH_MU), 0.001151215647092755)
        self.assertAlmostEqual(meanMotionToSma(0.001151215647092755, EARTH_MU), 6700)
