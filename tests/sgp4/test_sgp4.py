import os
import re
import unittest

from sattrack._sgp4 import TwoLineElement, WGS72

from sattrack.orbit import Satellite


class TestSGP4(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        p = os.path.dirname(__file__)

        tlePath = os.path.join(p, 'SGP4-VER.TLE')
        tleDict = {}
        with open(tlePath, 'r') as f:
            fileData = f.readlines()

        fileIter = iter(fileData)
        for row in fileIter:
            if row[0] == '1':
                row1 = row
                nxt = next(fileIter)
                row2, data = nxt[:69], nxt[69:].strip()
                s = f'{row1}{row2}'
                start, stop, step = re.split(r'\s+', data)
                tle = TwoLineElement(s, WGS72)
                tleDict[row2[2:7]] = (tle, eval(start), eval(stop), eval(step))
        cls._tleDict = tleDict

        resultsPath = os.path.join(p, 'tcppver.out')
        resultsDict = {}
        with open(resultsPath, 'r') as f:
            fileData = f.readlines()

        fileIter = iter(fileData)
        catNum = 0
        subDict = {}
        for row in fileIter:
            if row[-3:-1] == 'xx':
                if len(subDict) != 0:
                    resultsDict[catNum] = subDict
                subDict = {}
                catNum = row[0:5]
                continue

            data = re.split(r'\s+', row.strip())
            nums = [eval(d) for d in data[:7]]
            dt, rx, ry, rz, vx, vy, vz = nums
            subDict[dt] = ((rx, ry, rz), (vx, vy, vz))
        cls._resultsDict = resultsDict

    def almostEqualVector(self, catNum, timeStamp, posVec, velVec, results):
        msg = 'test state of sat number: {} at timeStamp: {}, {{}}, {{}}-component'\
            .format(catNum, timeStamp)
        self.assertAlmostEqual(results[0][0], posVec[0], delta=1, msg=msg.format('position', 'x'))
        self.assertAlmostEqual(results[0][1], posVec[1], delta=1, msg=msg.format('position', 'y'))
        self.assertAlmostEqual(results[0][2], posVec[2], delta=1, msg=msg.format('position', 'z'))
        self.assertAlmostEqual(results[1][0], velVec[0], delta=1, msg=msg.format('velocity', 'x'))
        self.assertAlmostEqual(results[1][1], velVec[1], delta=1, msg=msg.format('velocity', 'y'))
        self.assertAlmostEqual(results[1][2], velVec[2], delta=1, msg=msg.format('velocity', 'z'))

    def test_tle_output(self):
        for catNum, resultsDict in self._resultsDict.items():
            try:
                tle, start, stop, step = self._tleDict[catNum]
                sat = Satellite(tle)
                jd = tle.epoch
                for epoch, results in resultsDict.items():
                    state = sat.getState(jd.future(epoch / 1440.0))
                    self.almostEqualVector(catNum, epoch, *state, results)
            except KeyError:
                pass
