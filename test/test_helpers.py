import unittest
from math import pi

from sattrack.util.helpers import atan3, signOf, computeAngleDifference, hasSameSign


# It's possible these might be used in the public API, so testing them is probably wise.
class TestHelpers(unittest.TestCase):

    def test_atan3(self):
        # Just need to focus on testing the signs of each quadrant.
        self.assertAlmostEqual(atan3(1, 1), pi / 4)
        self.assertAlmostEqual(atan3(1, -1), 3 * pi / 4)
        self.assertAlmostEqual(atan3(-1, -1), 5 * pi / 4)
        self.assertAlmostEqual(atan3(-1, 1), 7 * pi / 4)

    def test_signOf(self):
        self.assertEqual(signOf(-5), -1)
        self.assertEqual(signOf(1.2345), 1)
        self.assertEqual(signOf(0), 1)

    def test_angleDifference(self):
        self.assertAlmostEqual(computeAngleDifference(-pi / 4), -pi / 4)
        self.assertAlmostEqual(computeAngleDifference(pi / 4), pi / 4)
        self.assertAlmostEqual(computeAngleDifference(3 * pi / 4), 3 * pi / 4)
        self.assertAlmostEqual(computeAngleDifference(-3 * pi / 4), -3 * pi / 4)
        self.assertAlmostEqual(computeAngleDifference(5 * pi / 4), - 3 * pi / 4)
        self.assertAlmostEqual(computeAngleDifference(-5 * pi / 4), 3 * pi / 4)
        self.assertAlmostEqual(computeAngleDifference(9 * pi / 4), pi / 4)
        self.assertAlmostEqual(computeAngleDifference(-13 * pi / 4), 3 * pi / 4)

    def test_sameSign(self):
        self.assertTrue(hasSameSign(5, 15))
        self.assertTrue(hasSameSign(2, 4.5))
        self.assertTrue(hasSameSign(-pi, -1))
        self.assertFalse(hasSameSign(1.4, -88))
        self.assertFalse(hasSameSign(-1.23, 4.56))
