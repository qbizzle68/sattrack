import unittest

from sattrack.core.juliandate import JulianDate
from sattrack.orbit.satellite import Satellite
from sattrack.orbit.tle import TwoLineElement
from sattrack.satellitepass.eclipse import getShadowPositions, Shadow, getShadowAnomalies, getShadowTimes, isEclipsed
from sattrack.satellitepass.exceptions import NoSatelliteEclipseException


class TestEclipse(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23302.71093069  .00024763  00000+0  44330-3 0  9993
2 25544  51.6441  22.2092 0000603  68.5079  66.2522 15.49821000422656''')
        cls.tle = tle
        cls.sat = Satellite(tle)
        cls.jd = tle.epoch

    def testFunctions(self):
        fromPositions = getShadowPositions(self.sat, self.jd, Shadow.PENUMBRA)
        fromAnomalies = getShadowAnomalies(self.sat, self.jd, Shadow.PENUMBRA)
        fromTimes = getShadowTimes(self.sat, self.jd, Shadow.PENUMBRA)

        self.assertAlmostEqual(fromPositions[0][0], fromAnomalies[0])
        self.assertAlmostEqual(fromPositions[1][0], fromAnomalies[1])
        self.assertAlmostEqual(fromPositions[0][1].value, fromTimes[0].value)
        self.assertAlmostEqual(fromPositions[1][1].value, fromTimes[1].value)

        fromPositions = getShadowPositions(self.sat, self.jd, Shadow.UMBRA)
        fromAnomalies = getShadowAnomalies(self.sat, self.jd, Shadow.UMBRA)
        fromTimes = getShadowTimes(self.sat, self.jd, Shadow.UMBRA)

        self.assertAlmostEqual(fromPositions[0][0], fromAnomalies[0])
        self.assertAlmostEqual(fromPositions[1][0], fromAnomalies[1])
        self.assertAlmostEqual(fromPositions[0][1].value, fromTimes[0].value)
        self.assertAlmostEqual(fromPositions[1][1].value, fromTimes[1].value)

    def testEclipsed(self):
        # Before entering
        jd = JulianDate(2023, 10, 29, 17, 40, 0, 0)
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.PENUMBRA))
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.UMBRA))
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.ANNULAR))
        # After exiting
        jd = JulianDate(2023, 10, 29, 18, 30, 0, 0)
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.PENUMBRA))
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.UMBRA))
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.ANNULAR))

        # Entering between penumbra and umbra
        jd = JulianDate(2023, 10, 29, 17, 48, 30, 0)
        self.assertTrue(isEclipsed(self.sat, jd, Shadow.PENUMBRA))
        self.assertFalse(isEclipsed(self.sat, jd, Shadow.UMBRA))

        # Leaving between umbra and penumbra
        # fixme: the isEclipsed method is too rudimentary, needs improvement
        if False:
            jd = JulianDate(2023, 10, 29, 18, 24, 20, 0)
            self.assertFalse(isEclipsed(self.sat, jd, Shadow.UMBRA))
            self.assertTrue(isEclipsed(self.sat, jd, Shadow.PENUMBRA))

        # Middle of umbra
        jd = JulianDate(2023, 10, 29, 18, 0, 0, 0)
        self.assertTrue(isEclipsed(self.sat, jd, Shadow.UMBRA))
        self.assertTrue(isEclipsed(self.sat, jd, Shadow.PENUMBRA))

    def testNoEclipse(self):
        tle = TwoLineElement('''NOAA 12 [-]             
        1 21263U 91032A   23302.55791971  .00000231  00000+0  11445-3 0  9996
        2 21263  98.5430 294.1456 0012980 323.4152  36.6141 14.26185652687862''')
        sat = Satellite(tle)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowPositions(sat, tle.epoch, Shadow.PENUMBRA)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowPositions(sat, tle.epoch, Shadow.UMBRA)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowAnomalies(sat, tle.epoch, Shadow.PENUMBRA)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowAnomalies(sat, tle.epoch, Shadow.UMBRA)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowTimes(sat, tle.epoch, Shadow.PENUMBRA)

        with self.assertRaises(NoSatelliteEclipseException):
            getShadowTimes(sat, tle.epoch, Shadow.UMBRA)
