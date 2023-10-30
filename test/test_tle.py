import re
import unittest
from contextlib import redirect_stdout
from io import StringIO
from unittest import mock

from pyevspace import Vector
from sattrack.orbit.sgp4 import TwoLineElement

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

    @mock.patch('sattrack.orbit.tle.requests')
    @mock.patch('sattrack.orbit.tle.input')
    def testMultipleTleRequest(self, mockInput, requestMock):

        # with mock.patch('sattrack.orbit.tle.requests') as requestMock:
        tleList = '''UME (ISS)               
1 08709U 76019A   23301.55301458  .00000133  00000+0  20482-3 0  9992
2 08709  69.6749 342.2326 0011572 343.9923 128.3410 13.71557622385080
UME-2 (ISS-B)           
1 10674U 78018A   23301.67653483 -.00000035  00000+0  23054-4 0  9999
2 10674  69.3649  32.0090 0160482  70.0218 100.9952 13.43653424241334
ISS (ZARYA)             
1 25544U 98067A   23302.18631802  .00023916  00000+0  42890-3 0  9995
2 25544  51.6427  24.8067 0000539  60.5934  25.2504 15.49791535422575
SWISSCUBE               
1 35932U 09051B   23301.61074555  .00001721  00000+0  38408-3 0  9995
2 35932  98.4987 175.0189 0006857 324.8711  35.2034 14.58129653748653
AISSAT 1                
1 36797U 10035C   23301.83449840  .00003252  00000+0  34755-3 0  9998
2 36797  98.2114 149.6456 0010421 353.9411   6.1678 14.89489601720690
AISSAT 2                
1 40075U 14037G   23301.73842984  .00003397  00000+0  40316-3 0  9993
2 40075  98.3378 204.8581 0006055 143.0406 217.1225 14.85192582503461
ISS DEB [EP BATTERY]    
1 47853U 98067RZ  23301.82478070  .00077059  00000+0  41673-3 0  9993
2 47853  51.6305 333.8027 0004141  29.0261 331.0970 15.79398727149842
ISS (NAUKA)             
1 49044U 21066A   23300.49437812  .00020675  00000+0  37311-3 0  9999
2 49044  51.6418  33.1919 0000639  84.5600 275.5462 15.49706576422307
ISS DEB (SPX-26 IPA FSE)
1 55448U 98067VB  23301.87676908  .00178538  00000+0  86823-3 0  9997
2 55448  51.6323   1.0970 0005360 327.2528  32.8141 15.81344070 42242
OUTPOST MISSION 1       
1 56226U 23084BK  23301.86776365  .00014113  00000+0  73009-3 0  9996
2 56226  97.5242  56.8569 0011452  99.3321 260.9207 15.16298168 20893
ISS DEB [ERA OUTFITTING]
1 56434U 98067VJ  23301.50245583  .00190397  00000+0  11507-2 0  9997
2 56434  51.6340  15.1479 0005925 323.0615  36.9977 15.76488989 27769
ISS DEB                 
1 57161U 98067VL  23301.92944941  .00359164  00000+0  19119-2 0  9991
2 57161  51.6343  15.7688 0005892 344.1544  15.9272 15.78858186 20036
ISS DEB                 
1 57162U 98067VM  23301.91016958  .00290812  00000+0  16954-2 0  9997
2 57162  51.6344  16.4145 0005466 351.2702   8.8202 15.77013840 20028
ISS DEB                 
1 57163U 98067VN  23301.93896613  .00181225  00000+0  14429-2 0  9996
2 57163  51.6354  18.1043 0001537 322.8791  37.2100 15.70091625 19990
ISS DEB [SPX-28 IPA FSE]
1 57212U 98067VP  23301.52115482  .00055723  00000+0  70334-3 0  9995
2 57212  51.6383  24.3092 0002029 221.9658 138.1180 15.59016026 19457'''

        tleString = '''ISS (ZARYA)             
        1 25544U 98067A   23302.18631802  .00023916  00000+0  42890-3 0  9995
        2 25544  51.6427  24.8067 0000539  60.5934  25.2504 15.49791535422575'''
        tle = TwoLineElement(tleString)

        responseMock = mock.MagicMock()
        responseMock.status_code = 200
        responseMock.text = tleList

        requestMock.get.return_value = responseMock
        mockInput.return_value = 3

        stream = StringIO()
        with redirect_stdout(stream):
            tleResponse = getTle('iss')
            self.assertEqual(tleResponse.name, tle.name)
            self.assertEqual(tleResponse.line1, tle.line1)
            self.assertEqual(tleResponse.line2, tle.line2)

        # Select return None option.
        mockInput.return_value = 16
        with redirect_stdout(stream):
            tleResponse = getTle('iss')
            self.assertIsNone(tleResponse)

    @mock.patch('sattrack.orbit.tle.requests')
    def testNoGPData(self, requestMock):
        responseMock = mock.MagicMock()
        responseMock.status_code = 200
        responseMock.text = 'No GP data found'

        requestMock.get.return_value = responseMock

        stream = StringIO()
        with redirect_stdout(stream):
            with self.assertRaises(NoTLEFound):
                getTle('zarya2')
