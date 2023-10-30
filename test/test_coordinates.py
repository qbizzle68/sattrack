import unittest
from math import radians, pi

from pyevspace import Vector, Matrix

from sattrack.core.coordinates import GeoPosition, CelestialCoordinates, geocentricToGeodetic, geodeticToGeocentric, \
    computeSubPoint, AltAz
from sattrack.core.juliandate import J2000


class TestGeoPosition(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.geo = GeoPosition(38.1, -97.9, 0.5)

    def testInitialization(self):
        # Test key functions.
        self.assertAlmostEqual(self.geo.latitudeRadians, radians(38.1))
        self.assertAlmostEqual(self.geo.longitudeRadians, radians(-97.9))
        # Test _other values are set correctly
        self.assertAlmostEqual(self.geo.latitude, 38.1)
        self.assertAlmostEqual(self.geo.longitude, -97.9)

        self.assertEqual(self.geo.elevation, 0.5)

        # Test latitude out of range.
        with self.assertRaises(ValueError):
            GeoPosition(-100, 25)

    def testString(self):
        geoString = str(self.geo)
        self.assertEqual(geoString, 'latitude: 38.1, longitude: -97.9, elevation: 0.5')

    def testJson(self):
        geoJson = self.geo.toJson()
        json = '{"latitude": 38.1, "longitude": -97.9, "elevation": 0.5}'
        self.assertEqual(geoJson, json)

    def testGeocentric(self):
        self.assertAlmostEqual(self.geo.latitudeGeocentric, 0.6617116098246585)
        self.assertAlmostEqual(self.geo.radius, 6370.0355049089185)

    def testVectors(self):
        time = J2000
        self.assertSequenceEqual(self.geo.computeVelocityVector(),
                                 Vector(0.0, 0.46451036904122134, 0.0))
        self.assertSequenceEqual(self.geo.getZenithVector(time),
                                 Vector(-0.7861514538796758, -0.03510869339868142, 0.6170358751407488))
        self.assertSequenceEqual(self.geo.getPositionVector(time),
                                 Vector(-5020.582288599724, -224.21390100268943, 3914.1836880499277))

    def testReferenceFrame(self):
        time = J2000
        matrix = Matrix((-0.6164214792840293, 0.044614475679551506, -0.7861514538796757),
                        (-0.027528732044877716, -0.9990042785493157, -0.03510869339868141),
                        (-0.7869350219613372, 0.0, 0.6170358751407488))
        ref = self.geo.getReferenceFrame(time)
        for i in range(3):
            for j in range(3):
                with self.subTest(row=i, col=j):
                    self.assertAlmostEqual(ref.matrix[i][j], matrix[i][j])

    def testParts(self):
        latParts = (12, 34, 56.789)
        lat = latParts[0] + latParts[1] / 60 + latParts[2] / 3600
        lngParts = (-12, 34, 56.789)
        lng = lngParts[0] - lngParts[1] / 60 - lngParts[2] / 3600
        geo = GeoPosition(lat, lng)
        parts = geo.parts

        for angleExpected, angleValue in zip((latParts, lngParts), parts):
            for expected, actual in zip(angleExpected, angleValue):
                with self.subTest():
                    self.assertAlmostEqual(actual, expected)


class TestCelestialCoordinates(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.coords = CelestialCoordinates(12, -45)

    def testInit(self):
        self.assertAlmostEqual(self.coords.latitudeRadians, -pi / 4)
        self.assertAlmostEqual(self.coords.longitudeRadians, pi)

        # Test declination out of range.
        with self.assertRaises(ValueError):
            CelestialCoordinates(25, -100)

    def testString(self):
        coordString = 'right-ascension: 12, declination: -45'
        self.assertEqual(str(self.coords), coordString)

    def testJson(self):
        json = '{"right-ascension": 12, "declination": -45}'
        self.assertEqual(self.coords.toJson(), json)

    def testNewProperties(self):
        self.assertAlmostEqual(self.coords.rightAscension, self.coords.longitude)
        self.assertAlmostEqual(self.coords.declination, self.coords.latitude)

        self.assertAlmostEqual(self.coords.rightAscensionRadians, self.coords.longitudeRadians)
        self.assertAlmostEqual(self.coords.declinationRadians, self.coords.latitudeRadians)

    def testParts(self):
        latParts = (12, 34, 56.789)
        lat = latParts[0] + latParts[1] / 60 + latParts[2] / 3600
        lngParts = (-12, 34, 56.789)
        lng = lngParts[0] - lngParts[1] / 60 - lngParts[2] / 3600
        geo = CelestialCoordinates(lng, lat)
        parts = geo.parts

        for angleExpected, angleValue in zip((lngParts, latParts), parts):
            for expected, actual in zip(angleExpected, angleValue):
                with self.subTest():
                    self.assertAlmostEqual(actual, expected)


class TestCoordinatesModule(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.geo = GeoPosition(45, 180)
        cls.equator = GeoPosition(0, 0)
        cls.pole = GeoPosition(90, -180)

    def testToGeodetic(self):
        self.assertEqual(geocentricToGeodetic(0), 0)
        self.assertEqual(geocentricToGeodetic(90), 90)
        self.assertAlmostEqual(geocentricToGeodetic(45), 45.19242142177217)

        with self.assertRaises(ValueError):
            geocentricToGeodetic(100)
            geocentricToGeodetic(-100)

    def testToGeocentric(self):
        self.assertEqual(geodeticToGeocentric(0), 0)
        self.assertEqual(geodeticToGeocentric(90), 90)
        self.assertAlmostEqual(geodeticToGeocentric(45.19242142177217), 45)

        with self.assertRaises(ValueError):
            geodeticToGeocentric(100)
            geodeticToGeocentric(-100)

    def testSubPoint(self):
        position = Vector(-5020.582288599724, -224.21390100268943, 3914.1836880499277)
        geo = computeSubPoint(position, J2000)
        self.assertAlmostEqual(geo.latitude, 38.1)
        self.assertAlmostEqual(geo.longitude, -97.9)


class TestAltAz(unittest.TestCase):

    def testInit(self):
        # Test modulo of azimuth.
        altAz = AltAz(45, 405)
        self.assertAlmostEqual(altAz.azimuth, 45)
        self.assertAlmostEqual(altAz.altitude, 45)
        self.assertEqual(altAz.direction, 'NE')

        # Test altitude out of range.
        with self.assertRaises(ValueError):
            AltAz(100, 45)

    def testAzimuthString(self):
        self.assertEqual(AltAz.azimuthAngleString(0), 'N')
        self.assertEqual(AltAz.azimuthAngleString(15), 'NNE')
        self.assertEqual(AltAz.azimuthAngleString(35), 'NE')
        self.assertEqual(AltAz.azimuthAngleString(60), 'ENE')
        self.assertEqual(AltAz.azimuthAngleString(80), 'E')
        self.assertEqual(AltAz.azimuthAngleString(105), 'ESE')
        self.assertEqual(AltAz.azimuthAngleString(125), 'SE')
        self.assertEqual(AltAz.azimuthAngleString(150), 'SSE')
        self.assertEqual(AltAz.azimuthAngleString(175), 'S')
        self.assertEqual(AltAz.azimuthAngleString(200), 'SSW')
        self.assertEqual(AltAz.azimuthAngleString(215), 'SW')
        self.assertEqual(AltAz.azimuthAngleString(240), 'WSW')
        self.assertEqual(AltAz.azimuthAngleString(260), 'W')
        self.assertEqual(AltAz.azimuthAngleString(285), 'WNW')
        self.assertEqual(AltAz.azimuthAngleString(305), 'NW')
        self.assertEqual(AltAz.azimuthAngleString(330), 'NNW')
        self.assertEqual(AltAz.azimuthAngleString(355), 'N')

    def testString(self):
        altAz = AltAz(45, 45)
        ans = 'altitude: 45, azimuth: 45, direction: NE'
        self.assertEqual(str(altAz), ans)

    def testJson(self):
        altAz = AltAz(45, 45)
        ans = '{"altitude": 45, "azimuth": 45, "direction": "NE"}'
        self.assertEqual(altAz.toJson(), ans)
