from juliandate import JulianDate

class TwoLineElement:

    def __init__(self, tle: str):
        tokens = tle.split('\n')
        if len(tokens) != 3:
            raise Exception("Incorrect number of lines")  # make a custom version for this
        for i in range(1, 3):
            if not self._checksum(tokens[i]):
                raise Exception(f"Checksum failed for line {tokens[i]}")  # make a custom version for this
        line1Tokens, line2Tokens = self._checkTokens(tokens[1], tokens[2])
        self._parseLines(line1Tokens, line2Tokens)

    def _checksum(self, line: str) -> bool:
        checksum = 0
        for ch in line[:-1]:
            val = ord(ch)
            if val in range(0x30, 0x3A):
                checksum += (val - 0x30)
            elif val == 0x2D:
                checksum += 1
        return True if (checksum % 10) == (ord(line[-1]) - 0x30) else False

    def _checkTokens(self, line1: str, line2: str) -> tuple[list[str], list[str]]:
        line1Tokens = [i for i in line1.split(' ') if i != '']
        if len(line1Tokens) != 9:
            raise Exception(f"Bad number of tokens in line 1 {len(line1Tokens)}.")  # make custom version for this
        line1lens = [(1, 2), (6, 7), (6, 9), (14, 15), (9, 11), (7, 9), (7, 9), (1, 2), (2, 6)]
        for tok, vals in zip(line1Tokens, line1lens):
            if not len(tok) in range(vals[0], vals[1]):
                raise Exception(f"Token of bad length {tok}, {vals}")  # make a custom version and describe which token

        line2Tokens = [i for i in line2.split(' ') if i != '']
        if not len(line2Tokens) in [8, 9]:
            raise Exception(f"Bad number of tokens in line 2. {len(line2Tokens)}y")  # make custom version for this
        line2lens = [(1, 2), (5, 6), (6, 9), (6, 9), (7, 8), (6, 9), (6, 9)]
        for tok, vals in zip(line2Tokens, line2lens):
            if not len(tok) in range(vals[0], vals[1]):
                raise Exception(f"Token of bad length {tok}, {vals}")  # make a custom version and describe which token
        if len(line2Tokens) == 8:
            if not len(line2Tokens[-1]) in range(16, 18):
                raise Exception(f"bad token {line2Tokens[-1]}, {(16, 17)}")  # make a custom version and describe which token
        else:
            for tok, vals in zip(line2Tokens[-2:], [(10, 12), (2, 6)]):
                if not len(tok) in range(vals[0], vals[1]):
                    raise Exception(f"Token of bad length {tok}, {vals}")  # make a custom version and describe which token
        return line1Tokens, line2Tokens

    def _parseLines(self, line1Tokens: list[str], line2Tokens: list[str]) -> None:
        #   LINE 1
        self._catNum = int(line1Tokens[1][:-1])
        self._class = line1Tokens[1][-1]
        self._cospar = line1Tokens[2]
        epochYear = int(line1Tokens[3][:2])
        epochYear += 2000 if epochYear < 57 else 1900
        epochDay = float(line1Tokens[3][2:])
        self._epoch = JulianDate(12, 31, epochYear - 1, 0, 0, 0).future(epochDay)
        self._nDot = float(line1Tokens[4])
        self._nDDot = float('.' + line1Tokens[5][:-2] + 'e' + line1Tokens[5][-2:])
        self._bStar = float('.' + line1Tokens[6][:-2] + 'e' + line1Tokens[6][-2:])
        self._ephemeris = int(line1Tokens[7])
        self._setNumber = int(line1Tokens[8][:-1])
        #   line 2
        self._inc = float(line2Tokens[2])
        self._raan = float(line2Tokens[3])
        self._ecc = float('0.' + line2Tokens[4])
        self._aop = float(line2Tokens[5])
        self._meanAnom = float(line2Tokens[6])
        if len(line2Tokens) == 8:
            self._meanMotion = float(line2Tokens[7][:11])
            self._revNum = int(line2Tokens[7][11:-1])
        else:
            self._meanMotion = float(line2Tokens[7])
            self._revNum = int(line2Tokens[8][:-1])
