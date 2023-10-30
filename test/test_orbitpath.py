import unittest

from sattrack.orbit.sgp4 import TwoLineElement

from sattrack.core.coordinates import GeoPosition
from sattrack.core.juliandate import JulianDate
from sattrack.orbit.orbitpath import OrbitPath
from sattrack.orbit.satellite import Satellite


class TestOrbitPath(unittest.TestCase):

    tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   23301.73882576  .00023270  00000+0  41790-3 0  9992
2 25544  51.6433  27.0249 0000531  82.8568  24.6803 15.49767794422509''')
    iss = Satellite(tle)

    def testNormalISS(self):

        jd = JulianDate(2023, 10, 28, 12, 0, 0, -5)
        geo = GeoPosition(38, -97)

        path = OrbitPath(self.iss, geo)
        timeInputs = [jd.future(i) for i in range(10)]
        answers = [
            (JulianDate(2023, 10, 28, 1, 13, 48.78761887549081, -5),
             JulianDate(2023, 10, 28, 12, 35, 27.601347863674164, -5)),
            (JulianDate(2023, 10, 29, 0, 50, 22.64815270899453, -5),
             JulianDate(2023, 10, 29, 12, 12, 0.4846993088722229, -5)),
            (JulianDate(2023, 10, 31, 0, 3, 30.681670904146756, -5),
             JulianDate(2023, 10, 31, 11, 25, 6.618046760559082, -5)),
            (JulianDate(2023, 10, 31, 23, 40, 4.871512949466705, -5),
             JulianDate(2023, 11, 1, 11, 1, 39.853357672691345, -5)),
            (JulianDate(2023, 11, 1, 23, 16, 39.18607771396637, -5),
             JulianDate(2023, 11, 2, 10, 38, 13.190740048885345, -5)),
            (JulianDate(2023, 11, 2, 22, 53, 13.632245063781738, -5),
             JulianDate(2023, 11, 3, 10, 14, 46.622147262096405, -5)),
            (JulianDate(2023, 11, 3, 22, 29, 48.21560740470886, -5),
             JulianDate(2023, 11, 4, 9, 51, 20.13840615749359, -5)),
            (JulianDate(2023, 11, 4, 22, 6, 22.99039900302887, -5),
             JulianDate(2023, 11, 5, 9, 27, 53.779186606407166, -5)),
            (JulianDate(2023, 11, 5, 21, 42, 57.859215438365936, -5),
             JulianDate(2023, 11, 6, 9, 4, 27.433284223076043, -5)),
            (JulianDate(2023, 11, 6, 21, 19, 32.873273491859436, -5),
             JulianDate(2023, 11, 7, 8, 41, 1.1384376883470395, -5))]

        for time, answer in zip(timeInputs, answers):
            with self.subTest():
                times = path.computeOrbitPassTimes(time)
                self.assertAlmostEqual(times[0].value, answer[0].value)
                self.assertAlmostEqual(times[1].value, answer[1].value)

        riseTime, setTime = path.computeSatellitePassTimes(jd, True)
        riseAnswer = JulianDate(2023, 10, 29, 2, 18, 22.019396424280785, -5)
        setAnswer = JulianDate(2023, 10, 29, 2, 28, 46.29301697014489, -5)
        self.assertAlmostEqual(riseTime.value, riseAnswer.value)
        self.assertAlmostEqual(setTime.value, setAnswer.value)

        riseTime, setTime = path.computeSatellitePassTimes(jd, False)
        riseAnswer = JulianDate(2023, 10, 28, 11, 13, 37.966099083423615, -5)
        setAnswer = JulianDate(2023, 10, 28, 11, 23, 35.074507892131805, -5)
        self.assertAlmostEqual(riseTime.value, riseAnswer.value)
        self.assertAlmostEqual(setTime.value, setAnswer.value)

    # @unittest.skip
    def testTwoRootDisappear(self):
        tle = TwoLineElement('''EROS C3
1 54880U 22179A   23147.44232990  .00006430  00000+0  20916-3 0  9999
2 54880 139.3581   8.5149 0012925 260.5407  99.3945 15.29442116 22656''')
        sat = Satellite(tle)

        # Case 1-a.
        geo = GeoPosition(62.5818464830717, 119.94391125013482, 0)
        jd = JulianDate(2023, 5, 31, 22, 20, 13.789000000004307, -5.0)
        path = OrbitPath(sat, geo)
        times = path.computeOrbitPassTimes(jd)
        self.assertFalse(times)

        # Case 1-b.
        geo = GeoPosition(62.5818464830717, 3.0543140675080167, 0)
        jd = JulianDate(2023, 5, 31, 22, 34, 44.8859999999986, -5.0)
        path = OrbitPath(sat, geo)
        times = path.computeOrbitPassTimes(jd)
        self.assertFalse(times)

        # Case 1-c.
        geo = GeoPosition(62.5818464830717, 3.0543140675080167, 0)
        jd = JulianDate(2023, 5, 31, 22, 34, 35.8859999999986, -5.0)
        path = OrbitPath(sat, geo)
        times = path.computeOrbitPassTimes(jd)
        self.assertFalse(times)
