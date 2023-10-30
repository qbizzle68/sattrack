import unittest

from sattrack.bodies.body import BodyOrbitController, Body
from sattrack.bodies.earth import Earth, EarthController
from sattrack.bodies.sun import Sun
from sattrack.util.constants import EARTH_MU, EARTH_EQUITORIAL_RADIUS, EARTH_POLAR_RADIUS, EARTH_FLATTENING, \
    EARTH_SIDEREAL_PERIOD


# This is to test aspects of the abstract Body class that are not reached with the current instantiated Body types.
class TestBody(unittest.TestCase):

    def test_classMethod(self):
        class DummyController(BodyOrbitController):
            def __init__(self):
                pass

            def computePosition(self, time):
                pass

            def computeCelestialCoordinates(self, time):
                pass

            def computeAltAz(self, geo, time):
                pass

            @classmethod
            def dummyClassMethod(cls, arg):
                return arg

        dummy = Body('dummy', 0, 1, 0, DummyController)
        self.assertEqual(dummy.dummyClassMethod('success'), 'success')

    def test_properties(self):
        self.assertEqual(Earth.name, 'Earth')
        self.assertEqual(Earth.mu, EARTH_MU)
        self.assertEqual(Earth.Re, EARTH_EQUITORIAL_RADIUS)
        self.assertEqual(Earth.Rp, EARTH_POLAR_RADIUS)
        self.assertAlmostEqual(Earth.flattening, EARTH_FLATTENING)
        self.assertEqual(Earth.period, EARTH_SIDEREAL_PERIOD)
        self.assertEqual(Earth.parent, Sun)
        self.assertIsInstance(Earth.controller, EarthController)
