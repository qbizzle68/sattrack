import unittest
from unittest import mock

from sattrack.core.coordinates import GeoPosition, AltAz
from sattrack.core.juliandate import JulianDate
from sattrack.orbit.exceptions import SatelliteAlwaysAbove, NoPassException
from sattrack.orbit.orbitpath import OrbitPath
from sattrack.orbit.satellite import Satellite
# noinspection PyUnresolvedReferences
from sattrack.orbit.tle import TwoLineElement
from sattrack.satellitepass.info import Visibility
from sattrack.satellitepass.satpass import PassController


class TestPassController(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23302.71093069  .00024763  00000+0  44330-3 0  9993
2 25544  51.6441  22.2092 0000603  68.5079  66.2522 15.49821000422656''')
        cls.tle = tle
        cls.iss = Satellite(tle)
        cls.jd = tle.epoch
        cls.geo = GeoPosition(38, -97)

    def testBasicPass(self):
        pc = PassController(self.iss, self.geo)
        # Use pass list to test that code as well.
        plist = pc.getPassList(self.jd, 1)
        np = plist[0]

        riseTime = JulianDate(2023, 10, 29, 17, 4, 47.24222928285599, 0.0)
        riseAltAz = AltAz(0.0, 257.3565371228019)
        riseVisibility = Visibility(True, False)
        self.assertEqual(np.riseInfo.time, riseTime)
        self.assertAlmostEqual(np.riseInfo.altAz.altitude, riseAltAz.altitude)
        self.assertAlmostEqual(np.riseInfo.altAz.azimuth, riseAltAz.azimuth)
        self.assertEqual(np.riseInfo.visibility.visible, riseVisibility.visible)
        self.assertEqual(np.riseInfo.visibility.illuminated, riseVisibility.illuminated)
        self.assertEqual(np.riseInfo.visibility.unobscured, riseVisibility.unobscured)

        setTime = JulianDate(2023, 10, 29, 17, 8, 10.866391360759735, 0.0)
        setAltAz = AltAz(0.0, 220.1725996642073)
        setVisibility = Visibility(True, False)
        self.assertEqual(np.setInfo.time, setTime)
        self.assertEqual(np.setInfo.altAz.altitude, setAltAz.altitude)
        self.assertEqual(np.setInfo.altAz.azimuth, setAltAz.azimuth)
        self.assertEqual(np.setInfo.visibility.visible, setVisibility.visible)
        self.assertEqual(np.setInfo.visibility.illuminated, setVisibility.illuminated)
        self.assertEqual(np.setInfo.visibility.unobscured, setVisibility.unobscured)

        maxTime = JulianDate(2023, 10, 29, 17, 6, 29.05431032180786, 0.0)
        maxAltAz = AltAz(1.000219925243141, 236.1121919959671)
        maxVisibility = Visibility(True, False)
        self.assertEqual(np.maxInfo.time, maxTime)
        self.assertEqual(np.maxInfo.altAz.altitude, maxAltAz.altitude)
        self.assertEqual(np.maxInfo.altAz.azimuth, maxAltAz.azimuth)
        self.assertEqual(np.maxInfo.visibility.visible, maxVisibility.visible)
        self.assertEqual(np.maxInfo.visibility.illuminated, maxVisibility.illuminated)
        self.assertEqual(np.maxInfo.visibility.unobscured, maxVisibility.unobscured)

        self.assertIsNone(np.firstVisibleInfo)
        self.assertIs(np.firstIlluminatedInfo, np.riseInfo)
        self.assertIsNone(np.firstUnobscuredInfo)
        self.assertIsNone(np.maxVisibleInfo)
        self.assertIs(np.maxIlluminatedInfo, np.maxInfo)
        self.assertIsNone(np.maxUnobscuredInfo)
        self.assertIsNone(np.lastVisibleInfo)
        self.assertIs(np.lastIlluminatedInfo, np.setInfo)
        self.assertIsNone(np.lastUnobscuredInfo)

        self.assertEqual(np.name, 'ISS (ZARYA)')
        self.assertEqual(np.visibility.visible, False)
        self.assertEqual(np.visibility.illuminated, True)
        self.assertEqual(np.visibility.unobscured, False)
        self.assertEqual(np.visible, False)
        self.assertEqual(np.illuminated, True)
        self.assertEqual(np.unobscured, False)

    @unittest.skip('this will get complicated quickly, handle it later')
    @mock.patch('sattrack.satellitepass.satpass.OrbitPath.computeSatellitePassTimes')
    @mock.patch('sattrack.satellitepass.satpass.getShadowTimes')
    @mock.patch('sattrack.satellitepass.satpass.Sun')
    def testComplexPass(self, sunMock, shadowMock, pathMock):
        pathMock.return_value = (JulianDate(2023, 10, 17, 16, 58, 44.75004941225052, 0.0),
                                 JulianDate(2023, 10, 17, 17, 7, 17.914878129959106, 0.0))
        shadowMock.return_value = (JulianDate(2023, 10, 17, 17, 6, 0, 0),
                                   JulianDate(2023, 10, 17, 17, 0, 0, 0))
        sunMock.computeRiseSetTimes.return_value = (JulianDate(2023, 10, 17, 17, 4, 0, 0),
                                                    JulianDate(2023, 10, 17, 17, 2, 0, 0))

        tle = TwoLineElement('''ISS (ZARYA)
1 25544U 98067A   23290.70730211  .00020227  00000+0  35718-3 0  9996
2 25544  51.6416  81.6490 0004772 114.5122 328.6540 15.50333434420790''')
        sat = Satellite(tle)
        jd = tle.epoch
        pc = PassController(sat, self.geo)

        np = pc.getNextPass(jd)
        self.assertTrue(True)

    @mock.patch('sattrack.satellitepass.satpass.OrbitPath.computeSatellitePassTimes')
    @mock.patch('sattrack.satellitepass.satpass.PassController._getMaxTime')
    def testPassString(self, maxMock, pathMock):
        pathMock.return_value = (JulianDate(2023, 10, 17, 16, 58, 44.75004941225052, 0.0),
                                 JulianDate(2023, 10, 17, 17, 7, 17.914878129959106, 0.0))
        maxMock.return_value = JulianDate(2023, 10, 17, 17, 3, 0, 0)

        tle = TwoLineElement('''ISS (ZARYA)
        1 25544U 98067A   23290.70730211  .00020227  00000+0  35718-3 0  9996
        2 25544  51.6416  81.6490 0004772 114.5122 328.6540 15.50333434420790''')
        sat = Satellite(tle)
        jd = tle.epoch
        pc = PassController(sat, self.geo)

        np = pc.getNextPass(jd)
        answer = '                Pass details for ISS (ZARYA), at 2023/10/17 17:03:01.332 +0.0 UTC                '\
        '\n     instance      |     time     | altitude |   azimuth    | illuminated | unobscured | visible '\
        '\n-------------------------------------------------------------------------------------------------'\
        '\n       rise        | 16:58:44.75  |   0.00   | 318.38 (NW)  |    True     |   False    |  False  '\
        '\n-------------------------------------------------------------------------------------------------'\
        '\n        max        | 17:03:01.332 |   9.34   | 10.54   (N)  |    True     |   False    |  False  '\
        '\n-------------------------------------------------------------------------------------------------'\
        '\n        set        | 17:07:17.915 |   0.00   | 63.59  (ENE) |    True     |   False    |  False  \n'

        self.assertEqual(str(np), answer)

    @mock.patch('sattrack.satellitepass.satpass.OrbitPath.computeSatellitePassTimes')
    @mock.patch('sattrack.satellitepass.satpass.PassController._getMaxTime')
    def testJson(self, maxMock, pathMock):
        pathMock.return_value = (JulianDate(2023, 10, 17, 16, 58, 44.75004941225052, 0.0),
                                 JulianDate(2023, 10, 17, 17, 7, 17.914878129959106, 0.0))
        maxMock.return_value = JulianDate(2023, 10, 17, 17, 3, 0, 0)

        tle = TwoLineElement('''ISS (ZARYA)
            1 25544U 98067A   23290.70730211  .00020227  00000+0  35718-3 0  9996
            2 25544  51.6416  81.6490 0004772 114.5122 328.6540 15.50333434420790''')
        sat = Satellite(tle)
        jd = tle.epoch
        pc = PassController(sat, self.geo)

        np = pc.getNextPass(jd)
        answer = '{"name": "ISS (ZARYA)", "riseInfo": {"altitude": 0.0, "azimuth": 318.38147029934737, ' \
                 '"direction": "NW", "time": {"dayNumber": 2460235, "dayFraction": 0.2074623848311603}, ' \
                 '"illuminated": true, "unobscured": false, "visible": false}, "setInfo": {"altitude": 0.0, ' \
                 '"azimuth": 63.58782096285168, "direction": "ENE", "time": {"dayNumber": 2460235, ' \
                 '"dayFraction": 0.2134017925709486}, "illuminated": true, "unobscured": false, "visible": false}, ' \
                 '"maxInfo": {"altitude": 9.338585638749889, "azimuth": 10.53764999084694, "direction": "N", "time": ' \
                 '{"dayNumber": 2460235, "dayFraction": 0.2104320889338851}, "illuminated": true, "unobscured": ' \
                 'false, "visible": false}, "firstUnobscured": null, "lastUnobscured": null, "firstIlluminated": ' \
                 '{"altitude": 0.0, "azimuth": 318.38147029934737, "direction": "NW", "time": {"dayNumber": 2460235, ' \
                 '"dayFraction": 0.2074623848311603}, "illuminated": true, "unobscured": false, "visible": false}, ' \
                 '"lastIlluminated": {"altitude": 0.0, "azimuth": 63.58782096285168, "direction": "ENE", "time": ' \
                 '{"dayNumber": 2460235, "dayFraction": 0.2134017925709486}, "illuminated": true, "unobscured": ' \
                 'false, "visible": false}, "illuminated": true, "unobscured": false, "visible": false}'

        self.assertEqual(np.toJson(), answer)

    def testNthPass(self):
        jd = JulianDate(2023, 10, 30, 0, 0, 0)
        pc = PassController(self.iss, self.geo)
        plist = pc.getPassList(jd, 3)

        # This time is between plist index 9 and 10 (pass 10 and 11).
        middleTime = JulianDate(2023, 10, 31, 10, 0, 0)
        p = pc.getNextPass(middleTime, 1)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[10].maxInfo.time.value)
        p = pc.getNextPass(plist[10], 1)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[11].maxInfo.time.value)
        # Move to the next range of valid Orbit pass times.
        p = pc.getNextPass(middleTime, 8)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[17].maxInfo.time.value)
        p = pc.getNextPass(plist[10], 8)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[18].maxInfo.time.value)

        # Previous
        p = pc.getNextPass(middleTime, -1)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[9].maxInfo.time.value)
        p = pc.getNextPass(plist[9], -1)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[8].maxInfo.time.value)
        # Move to the previous range of valid Orbit pass times.
        p = pc.getNextPass(middleTime, -8)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[2].maxInfo.time.value)
        p = pc.getNextPass(plist[9], -8)
        self.assertAlmostEqual(p.maxInfo.time.value, plist[1].maxInfo.time.value)

    def testExceptions(self):
        tle = TwoLineElement('''NOAA 12 [-]             
1 21263U 91032A   23302.55791971  .00000231  00000+0  11445-3 0  9996
2 21263  98.5430 294.1456 0012980 323.4152  36.6141 14.26185652687862''')
        sat = Satellite(tle)
        pc = PassController(sat, self.geo)

        with self.assertRaises(ValueError):
            pc.getNextPass(self.jd, 0)

        # This will internally raise a NoSatelliteEclipseException, but will return a valid pass.
        np = pc.getNextPass(tle.epoch)
        self.assertTrue(np)

        # This will raise NoPassException internally, and return us an empty list.
        tle = TwoLineElement('''ISS (ZARYA)
1 25544U 98067A   23290.70730211  .00020227  00000+0  35718-3 0  9996
2 25544  51.6416  81.6490 0004772 114.5122 328.6540 15.50333434420790''')
        sat = Satellite(tle)
        geo = GeoPosition(85, 0)
        pc = PassController(sat, geo)
        plist = pc.getPassList(tle.epoch, 1)
        self.assertFalse(plist)

    def testAllPassBranches(self):
        jd = JulianDate(2023, 10, 17, 19, 0, 0, -5)
        spacewayTle = TwoLineElement('''SPACEWAY 2
1 28903U 05046B   23289.61779785  .00000119  00000+0  00000+0 0  9990
2 28903   3.0252  83.0468 0000855 140.0482 245.3432  1.00271510 41579''')
        spaceway = Satellite(spacewayTle)
        geo = GeoPosition(38, -97, 0)
        pc = PassController(spaceway, geo)

        # No orbit pass times, orbit path is up, sat is geo, sat is up
        with self.assertRaises(SatelliteAlwaysAbove):
            pc.getNextPass(jd)

        # No orbit pass times, orbit path is up, sat is geo, sat is not up
        geo = GeoPosition(38, 83, 0)
        pc = PassController(spaceway, geo)
        with self.assertRaises(NoPassException):
            pc.getNextPass(jd)

        # No orbit pass times, orbit path is up, sat is not geo
        customTle = TwoLineElement('''SPACEWAY 2
1 28903U 05046B   23289.61779785  .00000119  00000+0  00000+0 0  9990
2 28903   3.0252  83.0468 0000855 140.0482 245.3432  1.10271510 41570''')
        customSat = Satellite(customTle)
        geo = GeoPosition(38, -97, 0)
        pc = PassController(customSat, geo)
        np = pc.getNextPass(jd)
        self.assertTrue(np)

        # No orbit pass times, path is down
        ixpeTle = TwoLineElement('''IXPE
1 49954U 21121A   23289.01703978  .00006881  00000+0  53942-3 0  9995
2 49954   0.2310  48.2272 0010151 248.4626  63.2067 14.95117562101124''')
        ixpe = Satellite(ixpeTle)
        geo = GeoPosition(38.0, -97.0, 0)
        pc = PassController(ixpe, geo)
        with self.assertRaises(NoPassException):
            pc.getNextPass(jd)

        issTle = TwoLineElement('''ISS (ZARYA)
        1 25544U 98067A   23290.70730211  .00020227  00000+0  35718-3 0  9996
        2 25544  51.6416  81.6490 0004772 114.5122 328.6540 15.50333434420790''')
        iss = Satellite(issTle)
        geo = GeoPosition(75.0, -97.0, 0)
        pc = PassController(iss, geo)
        with self.assertRaises(NoPassException):
            pc.getNextPass(jd)

        geo = GeoPosition(-75.0, -97.0, 0)
        pc = PassController(iss, geo)
        with self.assertRaises(NoPassException):
            pc.getNextPass(jd)

        # Orbit pass times exist, sat is geo, sat is up.
        beidouTle = TwoLineElement('''BEIDOU-2 IGSO-1
1 36828U 10036A   23290.60952736 -.00000076  00000+0  00000+0 0  9993
2 36828  54.2324 171.7763 0037688 191.9314 346.6367  1.00298034 48443''')
        beidou = Satellite(beidouTle)
        geo = GeoPosition(38.0, 83.0, 0)
        pc = PassController(beidou, geo)
        np = pc.getNextPass(jd)
        self.assertTrue(np)

        # Orbit pass times exist, sat is geo, sat is not up.
        geo = GeoPosition(38.0, -97.0, 0)
        pc = PassController(beidou, geo)
        with self.assertRaises(NoPassException):
            pc.getNextPass(jd)

        # Orbit pass times exist, sat is not geo
        pc = PassController(iss, geo)
        np = pc.getNextPass(jd)
        self.assertTrue(np)

    def testProperties(self):
        pc = PassController(self.iss, self.geo)

        self.assertIs(pc.orbitable, self.iss)
        self.assertIs(pc.geo, self.geo)
        self.assertIsInstance(pc.path, OrbitPath)
