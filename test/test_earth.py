import unittest

from pyevspace import Vector

from sattrack.bodies.earth import Earth
from sattrack.bodies.position import JulianTimes
from sattrack.core.juliandate import J2000


class TestEarth(unittest.TestCase):

    def test_earthPosition(self):
        self.assertEqual(Earth.computePosition(J2000), Vector(0, 0, 0))

    def test_earthBaseMethods(self):
        self.assertIsNone(Earth.computeCelestialCoordinates(J2000))
        self.assertIsNone(Earth.computeAltAz(None, J2000))

    def test_siderealTime(self):
        self.assertAlmostEqual(Earth.computeMeanSiderealTime(J2000), 4.894961212735793)
        self.assertAlmostEqual(Earth.computeApparentSiderealTime(J2000), 4.894899280735596)
        self.assertAlmostEqual(Earth.computeLocalSiderealTime(J2000, 15), 5.156698668534745)
        self.assertAlmostEqual(Earth.computeLocalSiderealTime(J2000, -135), 2.5387047905432514)

        self.assertAlmostEqual(Earth.computeApparentSiderealTime(J2000), Earth.computeLocalSiderealTime(J2000, 0.0))

    def test_withJulianTimes(self):
        jt = JulianTimes(J2000)

        self.assertAlmostEqual(Earth.computeMeanSiderealTime(jt), 4.894961212735793)
        self.assertAlmostEqual(Earth.computeApparentSiderealTime(jt), 4.894899280735596)
        self.assertAlmostEqual(Earth.computeLocalSiderealTime(jt, 15), 5.156698668534745)
        self.assertAlmostEqual(Earth.computeLocalSiderealTime(jt, -135), 2.5387047905432514)
