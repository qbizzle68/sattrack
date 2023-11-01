import re

import requests

# noinspection PyUnresolvedReferences
from sattrack.orbit.sgp4 import TwoLineElement
from sattrack.orbit.exceptions import NoTLEFound


CELESTRAK_URL = "https://celestrak.com/NORAD/elements/gp.php"


class TLEResponseIterator:
    __slots__ = '_obj', '_n', '_length'

    def __init__(self, obj: list[str]):
        self._obj = obj
        self._n = 0
        self._length = len(obj)

    def __iter__(self):
        return self

    def __next__(self):
        if self._n >= self._length:
            raise StopIteration

        idx = self._n
        pattern = r'[12] \d{5}[A-Z]?'
        if not re.match(pattern, self._obj[idx]):
            # Assume there is a name associated with the TLE.
            inc = 3
        else:
            # Assume there is not a name associated with the TLE.
            inc = 2

        rtn = '\n'.join(self._obj[idx:idx + inc])
        self._n += inc

        return rtn


def getTle(value: str, query: str = 'NAME', *, limitOne=False, interactive=False) -> TwoLineElement | None:
    """Retrieves a TLE from the Celestrak online repository of continually updating TLEs via HTTP.

    Args:
        value: The value to search for, which depends on the querying type.
        query: The querying type (Default = 'name'), possible values are:
            CATNR: Catalog Number (1 to 9 digits). Allows return of data fora  single catalog number
            INTDES: International Designator (yyyy-nnn). Allows return of data for all objects associated with a
                    particular launch.
            GROUP:  Groups of satellites provided on the CelesTrak CurrentDate page.
            NAME:   Satellite Name (Default). Allows searching for satellites by parts of their name.
            SPECIAL:Special data sets for the GEO Protected Zone (GPZ) or GPZ Plus.
        limitOne: Only returns a single TLE if multiple are returned (returns the first TLE in response).
        interactive: If False and limitOne is not True, return a list of TwoLineElement object if the request returns
            multiple TLEs. If True and limitOne is not True, asks user in put interactively to select which TLE to
            select from the response."""

    params = {query.upper(): value, "FORMAT": "TLE"}
    response = requests.get(CELESTRAK_URL, params=params)

    response.raise_for_status()
    if response.text == 'No GP data found':
        raise NoTLEFound('No GP data found')

    tleIterator = TLEResponseIterator(response.text.splitlines())
    tleStrings = list(tleIterator)

    if len(tleStrings) == 1 or limitOne:
        return TwoLineElement(tleStrings[0])

    if interactive is False:
        # return [TwoLineElement(tle) for tle in tleStrings]
        tmp = []
        for tle in tleStrings:
            print(f"'{tle}'")
            tmp.append(TwoLineElement(tle))
        return tmp
    else:
        return _selectTleInteractive(tleStrings)


def _selectTleInteractive(tleStrings: list[str]) -> TwoLineElement | None:
    idx = 0
    tleCount = len(tleStrings)
    while idx < 1 or idx > tleCount + 1:
        for i, tle in enumerate(tleStrings, 1):
            # This undoes what we did iterating over the original response, but we need the parts.
            tleParts = tle.splitlines()
            print(f'{i}:', end='')
            if len(tleParts) == 3:
                print('', tleParts[0])
            else:
                print('\t', '\n\t'.join(tleParts))
        print(tleCount + 1, ': None')
        idx = int(input('Which TLE do you wish to import?: '))

    if idx == tleCount + 1:
        return None
    return TwoLineElement(tleStrings[idx - 1])
