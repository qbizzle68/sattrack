import unittest
from math import degrees

from _pyevspace import Vector

from sattrack.bodies.exceptions import SunRiseSetException
from sattrack.bodies.position import JulianTimes, computeNutationDeltas, computeTrueObliquity
from sattrack.bodies.sun import computeEarthHeliocentricLongitude, computeEarthHeliocentricLatitude, \
    computeSunDistance, computeSunGeocentricLongitude, computeSunGeocentricLatitude, computeSunApparentLongitude, \
    computeSunRightAscension, computeSunDeclination, computeSunCoordinates, Sun, computeSunPosition, computeSunData, \
    computeSunRiseSetTimes, computeSunTransitInfo, computeSunTwilightTimes, Twilight, computeTwilightType
from sattrack.core.coordinates import GeoPosition
from sattrack.core.juliandate import JulianDate


class TestSun(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        jd = JulianDate(2003, 10, 17, 12, 30, 30, -7)
        cls.jd = jd
        cls.jt = JulianTimes(jd, 67)
        cls.geo = GeoPosition(39.742476, -105.1786, 1.83014)

    def testHeliocentricAngles(self):
        lng = computeEarthHeliocentricLongitude(self.jt)
        lat = computeEarthHeliocentricLatitude(self.jt)
        self.assertAlmostEqual(degrees(lng), 24.0182616917, 5)
        self.assertAlmostEqual(degrees(lat), -0.0001011219, 5)

    def testRadius(self):
        self.assertAlmostEqual(computeSunDistance(self.jt), 0.9965422974, 5)

    def testGeocentricAngles(self):
        lng = computeSunGeocentricLongitude(self.jt)
        lat = computeSunGeocentricLatitude(self.jt)
        self.assertAlmostEqual(degrees(lng), 204.0182616917, 5)
        self.assertAlmostEqual(degrees(lat), 0.0001011219)

    def testNutation(self):
        nutationObliquity, nutationLongitude = computeNutationDeltas(self.jt)
        self.assertAlmostEqual(degrees(nutationObliquity), -0.00399840, 5)
        self.assertAlmostEqual(degrees(nutationLongitude), 0.00166657, 5)

    def testObliquity(self):
        trueObliquity = computeTrueObliquity(self.jt)
        self.assertAlmostEqual(degrees(trueObliquity), 23.440465, 5)

    def testApparentLongitude(self):
        apparentLongitude = computeSunApparentLongitude(self.jt)
        self.assertAlmostEqual(degrees(apparentLongitude), 204.0085519281, 5)

    def testCoordinates(self):
        # From individual methods.
        rightAscension = computeSunRightAscension(self.jt)
        declination = computeSunDeclination(self.jt)
        self.assertAlmostEqual(degrees(rightAscension), 202.22741, 5)
        self.assertAlmostEqual(degrees(declination), -9.31434, 5)

        # From single method.
        coords = computeSunCoordinates(self.jt)
        self.assertAlmostEqual(degrees(coords.rightAscensionRadians), 202.22741, 5)
        self.assertAlmostEqual(coords.declination, -9.31434, 5)

        # From Sun object.
        # This computation doesn't set deltat precisely, so we use less precision.
        coords = Sun.computeCelestialCoordinates(self.jd)
        self.assertAlmostEqual(degrees(coords.rightAscensionRadians), 202.22741, 3)
        self.assertAlmostEqual(coords.declination, -9.31434, 3)

    def testPosition(self):
        ans = Vector(-136182794.60902783, -55651344.97781926, -24128878.125889584)

        pos = computeSunPosition(self.jd)
        self.assertAlmostEqual(pos[0], ans[0])
        self.assertAlmostEqual(pos[1], ans[1])
        self.assertAlmostEqual(pos[2], ans[2])

        pos = Sun.computePosition(self.jd)
        self.assertAlmostEqual(pos[0], ans[0])
        self.assertAlmostEqual(pos[1], ans[1])
        self.assertAlmostEqual(pos[2], ans[2])

    @unittest.skip('this isn\'t implemented correctly yet')
    def testAltAz(self):
        pass

    def testAngleData(self):
        data = computeSunData(self.geo, self.jd)
        riseTime = JulianDate(2003, 10, 17, 7, 12, 0, -6)
        setTime = JulianDate(2003, 10, 17, 18, 18, 0, -6)
        transit = JulianDate(2003, 10, 17, 12, 46, 0, -6)
        altitude = 41

        self.assertAlmostEqual(data[0].value, riseTime.value, delta=1 / 1400)
        self.assertAlmostEqual(data[1].value, setTime.value, delta=1 / 1400)
        self.assertAlmostEqual(data[2].value, transit.value, delta=1 / 1400)
        self.assertAlmostEqual(degrees(data[3]), altitude, 0)

        # Functional based computation.
        riseComp, setComp = computeSunRiseSetTimes(self.geo, self.jd)
        self.assertAlmostEqual(riseComp.value, riseTime.value, delta=1 / 1400)
        self.assertAlmostEqual(setComp.value, setTime.value, delta=1 / 1400)

        transitComp, altComp = computeSunTransitInfo(self.geo, self.jd)
        self.assertAlmostEqual(transitComp.value, transit.value, delta=1 / 1400)
        self.assertAlmostEqual(degrees(altComp), altitude, 0)

        # Class based computation.
        riseComp, setComp = Sun.computeRiseSetTimes(self.geo, self.jd)
        self.assertAlmostEqual(riseComp.value, riseTime.value, delta=1 / 1400)
        self.assertAlmostEqual(setComp.value, setTime.value, delta=1 / 1400)

        transitComp, altComp = Sun.computeTransitInfo(self.geo, self.jd)
        self.assertAlmostEqual(transitComp.value, transit.value, delta=1 / 1400)
        self.assertAlmostEqual(degrees(altComp), altitude, 0)

        # Test the else branch in computeSunTwilightTimes
        riseComp, setComp = computeSunTwilightTimes(self.geo, self.jd, Twilight.Day)
        self.assertAlmostEqual(riseComp.value, riseTime.value, delta=1 / 1400)
        self.assertAlmostEqual(setComp.value, setTime.value, delta=1 / 1400)
        riseComp, setComp = computeSunTwilightTimes(self.geo, self.jd, Twilight.Night)
        self.assertAlmostEqual(riseComp.value, riseTime.value, delta=1 / 1400)
        self.assertAlmostEqual(setComp.value, setTime.value, delta=1 / 1400)

    def testSunAlways(self):
        north = GeoPosition(80, 0)
        south = GeoPosition(-80, 0)
        northSummer = JulianDate(2023, 6, 22, 0, 0, 0)
        northWinter = JulianDate(2023, 12, 23, 0, 0, 0)

        # Sun always up in north in summer.
        with self.assertRaises(SunRiseSetException):
            computeSunRiseSetTimes(north, northSummer)

        # Sun always down in south in winter.
        with self.assertRaises(SunRiseSetException):
            computeSunRiseSetTimes(south, northSummer)

        # Sun always down in north in winter.
        with self.assertRaises(SunRiseSetException):
            computeSunRiseSetTimes(north, northWinter)

        # Sun always up in south in summer.
        with self.assertRaises(SunRiseSetException):
            computeSunRiseSetTimes(south, northWinter)

    def testTwilightTimes(self):
        civilStartComp, civilEndComp = computeSunTwilightTimes(self.geo, self.jd, Twilight.Civil)
        civilStart = JulianDate(2003, 10, 17, 6, 45, 0, -6)
        civilEnd = JulianDate(2003, 10, 17, 18, 46, 0, -6)
        self.assertAlmostEqual(civilStartComp.value, civilStart.value, delta=1 / 1440)
        self.assertAlmostEqual(civilEndComp.value, civilEnd.value, delta=1 / 1440)

        # From Sun object.
        civilStartComp, civilEndComp = Sun.computeTwilightTimes(self.geo, self.jd, Twilight.Civil)
        self.assertAlmostEqual(civilStartComp.value, civilStart.value, delta=1 / 1440)
        self.assertAlmostEqual(civilEndComp.value, civilEnd.value, delta=1 / 1440)

        nauticalStartComp, nauticalEndComp = computeSunTwilightTimes(self.geo, self.jd, Twilight.Nautical)
        nauticalStart = JulianDate(2003, 10, 17, 6, 14, 0, -6)
        nauticalEnd = JulianDate(2003, 10, 17, 19, 17, 0, -6)
        self.assertAlmostEqual(nauticalStartComp.value, nauticalStart.value, delta=1 / 1440)
        self.assertAlmostEqual(nauticalEndComp.value, nauticalEnd.value, delta=1 / 1400)

        # From Sun object.
        nauticalStartComp, nauticalEndComp = Sun.computeTwilightTimes(self.geo, self.jd, Twilight.Nautical)
        self.assertAlmostEqual(nauticalStartComp.value, nauticalStart.value, delta=1 / 1440)
        self.assertAlmostEqual(nauticalEndComp.value, nauticalEnd.value, delta=1 / 1400)

        astroStartComp, astroEndComp = computeSunTwilightTimes(self.geo, self.jd, Twilight.Astronomical)
        astroStart = JulianDate(2003, 10, 17, 5, 42, 0, -6)
        astroEnd = JulianDate(2003, 10, 17, 19, 48, 0, -6)
        self.assertAlmostEqual(astroStartComp.value, astroStart.value, delta=1 / 1440)
        self.assertAlmostEqual(astroEndComp.value, astroEnd.value, delta=1 / 1440)

        # From Sun object.
        astroStartComp, astroEndComp = Sun.computeTwilightTimes(self.geo, self.jd, Twilight.Astronomical)
        self.assertAlmostEqual(astroStartComp.value, astroStart.value, delta=1 / 1440)
        self.assertAlmostEqual(astroEndComp.value, astroEnd.value, delta=1 / 1440)

    def testTwilightType(self):
        dayTime = JulianDate(2003, 10, 17, 12, 0, 0, -6)
        nightTime = JulianDate(2003, 10, 17, 0, 0, 0, -6)
        civilRise = JulianDate(2003, 10, 17, 7, 0, 0, -6)
        civilSet = JulianDate(2003, 10, 17, 18, 30, 0, -6)
        nauticalRise = JulianDate(2003, 10, 17, 6, 30, 0, -6)
        nauticalSet = JulianDate(2003, 10, 17, 19, 0, 0, -6)
        astroRise = JulianDate(2003, 10, 17, 6, 0, 0, -6)
        astroSet = JulianDate(2003, 10, 17, 19, 30, 0, -6)

        self.assertEqual(computeTwilightType(self.geo, dayTime), Twilight.Day)
        self.assertEqual(computeTwilightType(self.geo, nightTime), Twilight.Night)
        self.assertEqual(computeTwilightType(self.geo, civilRise), Twilight.Civil)
        self.assertEqual(computeTwilightType(self.geo, civilSet), Twilight.Civil)
        self.assertEqual(computeTwilightType(self.geo, nauticalRise), Twilight.Nautical)
        self.assertEqual(computeTwilightType(self.geo, nauticalSet), Twilight.Nautical)
        self.assertEqual(computeTwilightType(self.geo, astroRise), Twilight.Astronomical)
        self.assertEqual(computeTwilightType(self.geo, astroSet), Twilight.Astronomical)

        self.assertEqual(Sun.computeTwilightType(self.geo, dayTime), Twilight.Day)
        self.assertEqual(Sun.computeTwilightType(self.geo, nightTime), Twilight.Night)
        self.assertEqual(Sun.computeTwilightType(self.geo, civilRise), Twilight.Civil)
        self.assertEqual(Sun.computeTwilightType(self.geo, civilSet), Twilight.Civil)
        self.assertEqual(Sun.computeTwilightType(self.geo, nauticalRise), Twilight.Nautical)
        self.assertEqual(Sun.computeTwilightType(self.geo, nauticalSet), Twilight.Nautical)
        self.assertEqual(Sun.computeTwilightType(self.geo, astroRise), Twilight.Astronomical)
        self.assertEqual(Sun.computeTwilightType(self.geo, astroSet), Twilight.Astronomical)
