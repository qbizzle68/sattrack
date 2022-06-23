from sattrack.spacetime.juliandate import JulianDate
import requests


class TwoLineElement:
    """An implementation of the NORAD two-line elements sets used by the simplified general perturbation models to
    propagate satellite state vectors. The data in the TLE are 'mean elements' and only make sense when used with the
    appropriate model, namely the SGP4 and SDP4."""

    def __init__(self, tle: str):
        tokens = tle.splitlines()
        if len(tokens) != 3:
            raise Exception("Incorrect number of lines")  # make a custom version for this
        for i in range(1, 3):
            if not self._checksum(tokens[i]):
                raise Exception(f"Checksum failed for line {tokens[i]}")  # make a custom version for this
        line1Tokens, line2Tokens = self._checkTokens(tokens[1], tokens[2])
        self._parseLines(line1Tokens, line2Tokens)
        self._name = tokens[0]
        self._line1 = tokens[1]
        self._line2 = tokens[2]

    def __str__(self) -> str:
        return f'{self._name}\n{self._line1}\n{self._line2}'

    @staticmethod
    def _checksum(line: str) -> bool:
        """Checks the sum of each line with the checksum provided in the TLE."""
        checksum = 0
        for ch in line[:-1]:
            val = ord(ch)
            if val in range(0x30, 0x3A):
                checksum += (val - 0x30)
            elif val == 0x2D:
                checksum += 1
        return True if (checksum % 10) == (ord(line[-1]) - 0x30) else False

    @staticmethod
    def _checkTokens(line1: str, line2: str) -> tuple[list[str], list[str]]:
        """Checks the tokens of the lines and ensures they are of the right format."""
        line1Tokens = [i for i in line1.split(' ') if i != '']
        if len(line1Tokens) != 9:
            raise Exception(f"Bad number of tokens in line 1 {len(line1Tokens)}.")  # make custom version for this
        line1lens = [(1, 2), (6, 7), (6, 9), (14, 15), (9, 11), (7, 9), (7, 9), (1, 2), (2, 6)]
        for tok, vals in zip(line1Tokens, line1lens):
            if not len(tok) in range(vals[0], vals[1]):
                raise Exception(f"Token of bad length {tok}, {vals}")  # make a custom version and describe which token

        line2Tokens = [i for i in line2.split(' ') if i != '']
        if not len(line2Tokens) in [8, 9]:
            raise Exception(
                f"Bad number of tokens in line 2. {len(line2Tokens)}y")  # todo: make custom version for this
        line2lens = [(1, 2), (5, 6), (6, 9), (6, 9), (7, 8), (6, 9), (6, 9)]
        for tok, vals in zip(line2Tokens, line2lens):
            if not len(tok) in range(vals[0], vals[1]):
                raise Exception(
                    f"Token of bad length {tok}, {vals}")  # todo: make a custom version and describe which token
        if len(line2Tokens) == 8:
            if not len(line2Tokens[-1]) in range(16, 18):
                raise Exception(
                    f"bad token {line2Tokens[-1]}, {(16, 17)}")  # todo: make a custom version and describe
                # which token
        else:
            for tok, vals in zip(line2Tokens[-2:], [(10, 12), (2, 6)]):
                if not len(tok) in range(vals[0], vals[1]):
                    raise Exception(
                        f"Token of bad length {tok}, {vals}")  # todo: make a custom version and describe which token
        return line1Tokens, line2Tokens

    def _parseLines(self, line1Tokens: list[str], line2Tokens: list[str]) -> None:
        """Parse tokenized versions of the lines into their constituent elements."""
        #   LINE 1
        self._catNum = int(line1Tokens[1][:-1])
        self._class = line1Tokens[1][-1]
        self._cospar = line1Tokens[2]
        epochYear = int(line1Tokens[3][:2])
        epochYear += 2000 if epochYear < 57 else 1900
        epochDay = float(line1Tokens[3][2:])
        self._epoch = JulianDate(12, 31, epochYear - 1, 0, 0, 0).future(epochDay)
        self._nDot = float(line1Tokens[4])
        self._nDDot = float('0.' + line1Tokens[5][:-2] + 'e' + line1Tokens[5][-2:])
        self._bStar = float('0.' + line1Tokens[6][:-2] + 'e' + line1Tokens[6][-2:])
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

    def setName(self, name: str):
        self._name = name

    def getName(self) -> str:
        return self._name

    def getLine1(self) -> str:
        return self._line1

    def getLine2(self) -> str:
        return self._line2

    def catalogNumber(self) -> int:
        """Return the catalog number of the satellite."""
        return self._catNum

    def classification(self) -> str:
        """Return the classification of the satellite.
        U: Unclassified, C: Classified, S: Secret"""
        return self._class

    def cosparID(self) -> str:
        """Returns the COSPAR ID (international designator) of the satellite of the form:

            XXYYYZZZ, where:

            XX:     Last two digits of launch year.
            YYY:    Launch number of the year.
            ZZZ:    Piece of the launch (staring with single letters)."""
        return self._cospar

    def epoch(self) -> JulianDate:
        """Returns the epoch of the TLE as a Julian date."""
        return self._epoch

    def meanMotionDot(self) -> float:
        """Returns the first derivative of mean motion / 2, also known as the ballistic coefficient,
        in radians / s^2."""
        return self._nDot

    def meanMotionDDot(self) -> float:
        """Returns the second derivative of mean motion / 6 in radians / s^3."""
        return self._nDDot

    def bStar(self) -> float:
        """Returns the B*, the drag term, or radiation pressure coefficient."""
        return self._bStar

    def ephemeris(self) -> int:
        """Returns the ephemeris type, should always be 0 for distributed TLE data."""
        return self._ephemeris

    def setNumber(self) -> int:
        """Returns the element set number, incremented when a new TLE is generated."""
        return self._revNum

    def inclination(self) -> float:
        """Returns the inclination of the orbit in degrees."""
        return self._inc

    def raan(self) -> float:
        """Returns the right-ascension of the ascending node in degrees."""
        return self._raan

    def eccentricity(self) -> float:
        """Returns the eccentricity of the orbit."""
        return self._ecc

    def argumentOfPeriapsis(self) -> float:
        """Returns the argument of periapsis in degrees."""
        return self._aop

    def meanAnomaly(self) -> float:
        """Returns the mean anomaly in degrees."""
        return self._meanAnom

    def meanMotion(self) -> float:
        """Returns the mean motion in revolutions per day."""
        return self._meanMotion

    def revolutionNumber(self) -> int:
        """Returns the revolution number of the satellite."""
        return self._revNum


CELESTRAK_URL = "https://celestrak.com/NORAD/elements/gp.php?{}={}&FORMAT=TLE"


def getTLE(value: str, query: str = 'name') -> TwoLineElement | None:
    """Retrieves a TLE from the Celestrak online repository of continually updating TLEs via HTTP.
    Parameters:
    query:  The querying type, possible values are:
        CATNR: Catalog Number (1 to 9 digits). Allows return of data fora  single catalog number
        INTDES: International Designator (yyyy-nnn). Allows return of data for all objects
                associated with a particular launch.
        GROUP:  Groups of satellites provided on the CelesTrak CurrentDate page.
        NAME:   Satellite Name (Default). Allows searching for satellites by parts of their name.
        SPECIAL:Special data sets for the GEO Protected Zone (GPZ) or GPZ Plus:
    value:  The value to search for, which depends on the querying type."""
    response = requests.get(CELESTRAK_URL.format(query.upper(), value.replace(' ', '%20')))
    lines = response.text.splitlines()

    if response.text == 'No GP data found':
        raise Exception("No GP data found")
    elif len(lines) == 3:
        return TwoLineElement(response.text)
    else:
        lineGroups = [(lines[i], lines[i + 1], lines[i + 2]) for i in range(0, len(lines), 3)]
        print("Multiple TLEs found.")
        for l, i in zip(lineGroups, range(1, len(lineGroups) + 1)):
            print(f'{i}:  {l[0]}')
        print(f'{len(lineGroups) + 1}:  None')
        x = 0
        while x not in range(1, len(lineGroups) + 2):
            x = int(input("Which TLE do you wish to import?: "))
        if x == len(lineGroups) + 1:
            return None
        i, j, k = lineGroups[x - 1]
        tle = i + '\n' + j + '\n' + k
        return TwoLineElement(tle)
