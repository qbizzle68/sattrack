import requests

from sattrack.orbit.sgp4 import TwoLineElement
from sattrack.orbit.exceptions import NoTLEFound

_CELESTRAK_URL = "https://celestrak.com/NORAD/elements/gp.php?{}={}&FORMAT=TLE"


def getTle(value: str, query: str = 'name') -> TwoLineElement | None:
    """
    Retrieves a TLE from the Celestrak online repository of continually updating TLEs via HTTP.

    Args:
        value: The value to search for, which depends on the querying type.
        query: The querying type (Default = 'name'), possible values are:
            CATNR: Catalog Number (1 to 9 digits). Allows return of data fora  single catalog number
            INTDES: International Designator (yyyy-nnn). Allows return of data for all objects associated with a
                    particular launch.
            GROUP:  Groups of satellites provided on the CelesTrak CurrentDate page.
            NAME:   Satellite Name (Default). Allows searching for satellites by parts of their name.
            SPECIAL:Special data sets for the GEO Protected Zone (GPZ) or GPZ Plus.
    Returns:
        A TwoLineElement object from the search parameters.
    """

    response = requests.get(_CELESTRAK_URL.format(query.upper(), value.replace(' ', '%20')))
    lines = response.text.splitlines()

    if response.text == 'No GP data found':
        raise NoTLEFound("No GP data found")
    elif len(lines) == 3:
        return TwoLineElement(response.text.strip())
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
