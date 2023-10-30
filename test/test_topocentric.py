import unittest
from unittest import mock

from pyevspace import Vector
from sattrack.orbit.sgp4 import TwoLineElement

from sattrack.bodies.topocentric import toTopocentric, toTopocentricOffset, toTopocentricState, fromTopocentric, \
    getAltitude, getAzimuth
from sattrack.core.coordinates import GeoPosition
from sattrack.core.juliandate import JulianDate
from sattrack.orbit.satellite import Satellite
from test.close import vectorIsClose


class TestTopocentric(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23301.52764947  .00019202  00000+0  34666-3 0  9994
2 25544  51.6435  28.0714 0001006 114.0425 254.5243 15.49747991422472''')
        cls.sat = Satellite(tle)
        cls.geo = GeoPosition(38, -100)
        cls.jd = JulianDate(2023, 10, 28, 11, 0, 0, -5)

    def testToTopocentric(self):
        vec = Vector(1, 0, 0)
        topo = toTopocentric(vec, self.geo, self.jd)
        ans = Vector(-0.6146654841499265, -0.05685858361419406, -0.7867359430356059)
        self.assertTrue(vectorIsClose(topo, ans))

        topo = toTopocentricOffset(Vector(1000, 1000, 1000), self.geo, self.jd)
        ans = Vector(-1388.4114146273155, -1055.2408257786808, -6496.307187751361)
        self.assertTrue(vectorIsClose(topo, ans))

        state = self.sat.getState(self.jd)
        topoState = toTopocentricState(*state, self.geo, self.jd)
        ansPos = Vector(-4156.282885074589, -4759.2849048480775, -3852.7853838789615)
        ansVel = Vector(2.066486874293925, 1.8210211148306046, 6.824010869605797)
        self.assertTrue(vectorIsClose(topoState[0], ansPos))
        self.assertTrue(vectorIsClose(topoState[1], ansVel))

    def testFromTopocentric(self):
        topo = Vector(1000, 1000, 1000)
        vec = fromTopocentric(topo, self.geo, self.jd)
        ans = Vector(-6482.546445501111, -632.4347171146414, 3733.093663941972)
        self.assertTrue(vectorIsClose(vec, ans))

    @mock.patch('sattrack.bodies.topocentric.toTopocentricOffset')
    def testAltitude(self, topoMock):
        topoMock.return_value = Vector(1000, 0, 1000)
        altComp = getAltitude(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(altComp, 45.0)

        topoMock.return_value = Vector(1000, 1000, 0)
        altComp = getAltitude(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(altComp, 0.0)

        topoMock.return_value = Vector(0, 0, 1000)
        altComp = getAltitude(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(altComp, 90.0)

    @mock.patch('sattrack.bodies.topocentric.toTopocentricOffset')
    def testAzimuth(self, topoMock):
        topoMock.return_value = Vector(100, 100, 50)
        azComp = getAzimuth(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(azComp, 135)

        # Test below horizon.
        topoMock.return_value = Vector(100, 100, -50)
        azComp = getAzimuth(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(azComp, 135)

        # Test quadrant is still positive value.
        topoMock.return_value = Vector(-100, -100, 100)
        azComp = getAzimuth(self.sat, self.geo, self.jd)
        self.assertAlmostEqual(azComp, 315)


if __name__ == '__main__':
    unittest.main()
