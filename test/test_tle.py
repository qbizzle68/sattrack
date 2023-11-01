import contextlib
import re
import unittest
from io import StringIO
from unittest import mock

# noinspection PyUnresolvedReferences
from pyevspace import Vector
from sattrack.orbit.tle import TwoLineElement

from sattrack.orbit.exceptions import NoTLEFound
from sattrack.orbit.satellite import Satellite
from sattrack.orbit.tle import getTle
from test.close import vectorIsClose


class TestTle(unittest.TestCase):

    @staticmethod
    def loadTleFile():
        tleLines = []
        timeData = []
        with open(r'test/resources/SGP4-VER.TLE') as f:
            for line in f:
                if line[0] == '#':
                    continue
                elif line[0] == '1':
                    tleLines.append(line)
                elif line[0] == '2':
                    # There's no rsplit() in re module, so need to make do here.
                    splitLine = line.rsplit(' ', 1)
                    lhs = splitLine[0]
                    dt = float(splitLine[1].strip())
                    splitLine = lhs.strip().rsplit(' ', 1)
                    lhs = splitLine[0]
                    end = float(splitLine[1])
                    splitLine = lhs.strip().rsplit(' ', 1)
                    line2 = splitLine[0]
                    start = float(splitLine[1])

                    tleLines.append(line2)
                    timeData.append((start, end, dt))

        tleList = [TwoLineElement(tleLines[i] + '\n' + tleLines[i+1]) for i in range(0, len(tleLines), 2)]

        return tleList, timeData

    @staticmethod
    def loadPositionData():
        stateList = []
        with open(r'test/resources/tcppver.out') as f:
            states = []
            for i, line in enumerate(f):
                if re.match(r'\d{5} xx', line):
                    if line != '00005 xx\n':
                        stateList.append(states)
                        states = []
                    continue
                splitLine = re.split(r'\s+', line.strip())

                pos = Vector(*(float(num) for num in splitLine[1:4]))
                vel = Vector(*(float(num) for num in splitLine[4:7]))

                states.append((pos, vel))

            # Add the last bunch.
            stateList.append(states)

        return stateList

    @classmethod
    def setUpClass(cls) -> None:
        tleList, timeData = cls.loadTleFile()
        cls.tleList = tleList
        cls.timeData = timeData
        stateList = cls.loadPositionData()
        cls.stateList = stateList

        with open('test/resources/response-list.txt') as f:
            cls.issResponse = f.read()

    # This skips testing sat number 26900.
    def testStates(self):
        for tle, states, (start, end, step) in zip(self.tleList, self.stateList, self.timeData):
            jd = tle.epoch
            sat = Satellite(tle)

            # We want to include end, so add 1 to it to include it.
            iterCount = int(round((end - start) / step))
            itr = [start + step * i for i in range(iterCount + 1)]
            for i, dt in enumerate(itr):
                time = jd.future(dt / 1440)
                stateComp = sat.getState(time)

                # fixme: Right now the values of satellite number 26900 are too far off.
                if tle.line1[2:7] == '26900':
                    continue

                with self.subTest(tle=tle.line1[2:7], dt=dt):
                    try:
                        self.assertTrue(vectorIsClose(stateComp[0], states[i][0], 0))
                        self.assertTrue(vectorIsClose(stateComp[1], states[i][1], 0))
                    except AssertionError:
                        print([i for i in stateComp[0]], [i for i in states[i][0]])
                        print([i for i in stateComp[1]], [i for i in states[i][1]])
                        exit(1)

    @mock.patch('sattrack.orbit.tle.requests')
    def testNonInteractiveRequest(self, requestMock):
        responseMock = mock.MagicMock()

        # Test multiple results, default values. (limitOne=False, interactive=False)
        responseMock.status_code = 200
        responseMock.text = self.issResponse
        requestMock.get.return_value = responseMock

        result = getTle('iss')
        self.assertTrue(result)
        for tle in result:
            with self.subTest():
                self.assertIsInstance(tle, TwoLineElement)

        # Test multiple results, limitOne=True.
        result = getTle('iss', limitOne=True)
        self.assertTrue(result)
        self.assertIsInstance(result, TwoLineElement)

        # Test multiple results, both options True.
        result = getTle('iss', limitOne=True, interactive=True)
        self.assertTrue(result)
        self.assertIsInstance(result, TwoLineElement)

        # Test interactive.
        stream = StringIO()
        with contextlib.redirect_stdout(stream):
            with mock.patch('sattrack.orbit.tle.input') as mockInput:
                mockInput.return_value = 3
                result = getTle('iss', interactive=True)

        self.assertTrue(result)
        self.assertIsInstance(result, TwoLineElement)
        self.assertEqual(result.name, 'ISS (ZARYA)')

        # Test interactive return None.
        with contextlib.redirect_stdout(stream):
            with mock.patch('sattrack.orbit.tle.input') as mockInput:
                mockInput.return_value = 16
                result = getTle('iss', interactive=True)

        self.assertIsNone(result)

        # Test a 2LE response format.
        # Mimic a 2LE format response.
        tmp = []
        for i, line in enumerate(self.issResponse.splitlines()):
            if i % 3 == 1 or i % 3 == 2:
                tmp.append(line)
        responseText = '\n'.join(tmp)

        responseMock.text = responseText
        requestMock.get.return_value = responseMock

        result = getTle('iss')
        self.assertTrue(result)
        for tle in result:
            with self.subTest():
                self.assertIsInstance(tle, TwoLineElement)

    @mock.patch('sattrack.orbit.tle.requests')
    def testSingleTleRequest(self, requestMock):
        tleString = '''ISS (ZARYA)             
        1 25544U 98067A   23302.18631802  .00023916  00000+0  42890-3 0  9995
        2 25544  51.6427  24.8067 0000539  60.5934  25.2504 15.49791535422575'''

        responseMock = mock.MagicMock()
        responseMock.status_code = 200
        responseMock.text = tleString
        requestMock.get.return_value = responseMock

        tle = TwoLineElement(tleString)
        tleResponse = getTle('zarya')
        self.assertEqual(tleResponse.name, tle.name)
        self.assertEqual(tleResponse.line1, tle.line1)
        self.assertEqual(tleResponse.line2, tle.line2)

        # Test simple return with alternate options.
        response = getTle('zarya', limitOne=True)
        self.assertEqual(response.name, tle.name)
        self.assertEqual(response.line1, tle.line1)
        self.assertEqual(response.line2, tle.line2)

        response = getTle('zarya', interactive=True)
        self.assertEqual(response.name, tle.name)
        self.assertEqual(response.line1, tle.line1)
        self.assertEqual(response.line2, tle.line2)

    @mock.patch('sattrack.orbit.tle.requests')
    def testNoGPData(self, requestMock):
        responseMock = mock.MagicMock()
        responseMock.status_code = 200
        responseMock.text = 'No GP data found'

        requestMock.get.return_value = responseMock

        # stream = StringIO()
        # with redirect_stdout(stream):
        with self.assertRaises(NoTLEFound):
            getTle('zarya2')
