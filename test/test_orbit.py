import unittest
from math import pi

from pyevspace import Vector, Matrix

from sattrack.core.coordinates import GeoPosition
from sattrack.core.juliandate import JulianDate
from sattrack.orbit.elements import Elements
from sattrack.orbit.satellite import Orbit, Satellite
from sattrack.orbit.tle import TwoLineElement
from test.close import vectorIsClose


class TestOrbit(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        jd = JulianDate(2023, 10, 29, 12, 0, 0)
        cls.jd = jd
        elements = Elements(pi / 2, pi / 4, pi, 0.01, 6700, 0, jd)
        cls.elements = elements
        cls.orbit = Orbit(elements, 'orbit')

    def testProperties(self):
        self.assertEqual(self.orbit.name, 'orbit')

    def testOrbitAnomalies(self):
        anomaly = self.orbit.anomalyAtTime(self.jd.future(0.5), 'true')
        self.assertAlmostEqual(anomaly, 5.739947586418642)
        anomaly = self.orbit.anomalyAtTime(self.jd.future(0.5), 'mean')
        self.assertAlmostEqual(anomaly, 5.75021880414991)
        with self.assertRaises(ValueError):
            self.orbit.anomalyAtTime(self.jd, 'blah')

        time = self.orbit.timeToNextAnomaly(pi / 2, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 12, 22, 27.09484577178955, 0))
        time = self.orbit.timeToNextAnomaly(pi / 2, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 12, 22, 44.467473328113556, 0))
        with self.assertRaises(ValueError):
            self.orbit.timeToNextAnomaly(0, self.jd, 'blah')

        time = self.orbit.timeToPreviousAnomaly(pi / 2, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 10, 51, 29.22487199306488, 0))
        time = self.orbit.timeToPreviousAnomaly(pi / 2, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 10, 51, 46.59753978252411, 0))
        with self.assertRaises(ValueError):
            self.orbit.timeToPreviousAnomaly(0, self.jd, 'blah')

        time = self.orbit.timeToNearestAnomaly(pi, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 12, 45, 28.934986889362335, 0))
        time = self.orbit.timeToNearestAnomaly(pi, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 12, 45, 28.934986889362335, 0))
        with self.assertRaises(ValueError):
            self.orbit.timeToNearestAnomaly(0, self.jd, 'blah')

    def testState(self):
        state = self.orbit.getState(self.jd)
        self.assertTrue(vectorIsClose(state[0], Vector(0, -6633, 0)))
        self.assertTrue(vectorIsClose(state[1], Vector(5.508832636164127, 0, -5.508832636164127)))

        topoState = self.orbit.getTopocentricState(GeoPosition(38, -97), self.jd)
        self.assertTrue(vectorIsClose(topoState[0],
                                      Vector(-3537.1212706294814, 3372.669755228793, -10870.799633434872)))
        self.assertTrue(vectorIsClose(topoState[1], Vector(2.7679258634811816, -4.327057005371725, -5.40504388193948)))

    def testElements(self):
        elements = self.orbit.getElements(self.jd)

        self.assertAlmostEqual(elements.raan, self.elements.raan)
        self.assertAlmostEqual(elements.inc, self.elements.inc)
        self.assertAlmostEqual(elements.aop, self.elements.aop)
        self.assertAlmostEqual(elements.ecc, self.elements.ecc)
        self.assertAlmostEqual(elements.sma, self.elements.sma)
        self.assertAlmostEqual(elements.meanAnomaly, self.elements.meanAnomaly)
        self.assertAlmostEqual(elements.trueAnomaly, self.elements.trueAnomaly)

        elements = self.orbit.getElements(self.jd.future(0.5))
        self.assertAlmostEqual(elements.meanAnomaly, 5.75021880414991)
        self.assertAlmostEqual(elements.trueAnomaly, 5.739947586418642)
        self.assertAlmostEqual(elements.epoch.value, self.jd.future(0.5).value)

        # Test elements property.
        self.assertIs(self.orbit.elements, self.elements)

    def testReferenceFrame(self):
        frame = self.orbit.getReferenceFrame(self.jd)
        tmp = 0.7071067811865476
        matrix = Matrix((0, tmp, tmp), (-1, 0, 0), (0, -tmp, tmp))
        self.assertEqual(frame.matrix, matrix)

    def testApsides(self):
        self.assertAlmostEqual(self.orbit.getPeriapsis(self.jd), 6633.0)
        self.assertAlmostEqual(self.orbit.getApoapsis(self.jd), 6767.0)

    def testIsGeo(self):
        self.assertFalse(self.orbit.isGeosynchronous())

        elements = Elements(pi, 0, 0, 0, 42164.16963419928, 0, self.elements.epoch)
        orbit = Orbit(elements)
        self.assertTrue(orbit.isGeosynchronous())


class TestSatellite(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        issTle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23302.71093069  .00024763  00000+0  44330-3 0  9993
2 25544  51.6441  22.2092 0000603  68.5079  66.2522 15.49821000422656''')
        cls.issTle = issTle
        cls.iss = Satellite(issTle)
        cls.jd = issTle.epoch
        geoTle = TwoLineElement('''TDRS 3                  
1 19548U 88091B   23301.50102956 -.00000305  00000+0  00000+0 0  9997
2 19548  13.3438 347.5683 0037774 329.4957 211.0578  1.00269464115737''')
        cls.geoSat = Satellite(geoTle)
        cls.geo = GeoPosition(38, -97)

    def testProperties(self):
        self.assertEqual(self.iss.name, 'ISS (ZARYA)')
        self.assertEqual(self.iss.tle, self.issTle)

    def testAnomalies(self):

        anomaly = self.iss.anomalyAtTime(self.jd.future(0.5), 'true')
        self.assertAlmostEqual(anomaly, 5.863433066840235)
        anomaly = self.iss.anomalyAtTime(self.jd.future(0.5), 'mean')
        self.assertAlmostEqual(anomaly, 5.8634735325189595)
        with self.assertRaises(ValueError):
            self.iss.anomalyAtTime(self.jd, 'blah')


        time = self.iss.timeToNextAnomaly(pi / 2, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 17, 9, 51.479162871837616, 0.0))
        time = self.iss.timeToNextAnomaly(pi / 2, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 17, 9, 51.58598184585571, 0.0))
        with self.assertRaises(ValueError):
            self.iss.timeToNextAnomaly(0, self.jd, 'blah')

        time = self.iss.timeToPreviousAnomaly(pi / 2, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 15, 37, 5.373244285583496, 0.0))
        time = self.iss.timeToPreviousAnomaly(pi / 2, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 15, 37, 5.480063259601593, 0.0))
        with self.assertRaises(ValueError):
            self.iss.timeToPreviousAnomaly(0, self.jd, 'blah')

        time = self.iss.timeToNearestAnomaly(pi, self.jd, 'true')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 17, 33, 3.112451434135437, 0.0))
        time = self.iss.timeToNearestAnomaly(pi, self.jd, 'mean')
        self.assertAlmostEqual(time, JulianDate(2023, 10, 29, 17, 33, 3.112451434135437, 0.0))
        with self.assertRaises(ValueError):
            self.iss.timeToNearestAnomaly(0, self.jd, 'blah')

    def testState(self):
        state = self.iss.getState(self.jd)
        self.assertTrue(vectorIsClose(state[0], Vector(-5562.624299946008, 958.2688858116214, 3775.4430675013127)))
        self.assertTrue(vectorIsClose(state[1], Vector(-3.7711857438936827, -5.1564500912660085, -4.234379662212585)))

        topoState = self.iss.getTopocentricState(self.geo, self.jd)
        self.assertTrue(vectorIsClose(topoState[0], Vector(114.75410015890384, -2516.753840201706, -64.2762481688684)))
        self.assertTrue(vectorIsClose(topoState[1], Vector(6.360056538268955, 3.4864179492270213, 1.262728487402717)))

    def testReferenceFrame(self):
        frame = self.iss.getReferenceFrame(self.jd)
        matrix = Matrix((0.12094180935084459, -0.9473727298885425, 0.29640848404604575),
                        (0.6730447425414484, -0.14122361818090612, -0.7259935703607913),
                        (0.7296463892107323, 0.28729914780838567, 0.6205443951727515))
        # So vectorIsClose just iterates over the vector, so it will work for any iterable.
        self.assertTrue(vectorIsClose(frame.matrix[0], matrix[0]))
        self.assertTrue(vectorIsClose(frame.matrix[1], matrix[1]))
        self.assertTrue(vectorIsClose(frame.matrix[2], matrix[2]))

    def testApsides(self):
        self.assertAlmostEqual(self.iss.getPeriapsis(self.jd), 6787.879667413829)
        self.assertAlmostEqual(self.iss.getApoapsis(self.jd), 6788.698335067378)

    def testIsGeo(self):
        self.assertFalse(self.iss.isGeosynchronous())
        self.assertTrue(self.geoSat.isGeosynchronous())
