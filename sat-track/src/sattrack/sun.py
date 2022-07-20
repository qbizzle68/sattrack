from enum import Enum
from functools import total_ordering
from math import sin, cos, pi, radians, tan, asin, degrees, sqrt

from pyevspace import EVector

from sattrack.spacetime.juliandate import JulianDate, J2000
from sattrack.util.constants import AU, TWOPI
from sattrack.util.conversions import atan3

DELTAT = 157.2  # seconds from https://webspace.science.uu.nl/~gent0113/deltat/deltat.htm


@total_ordering
class TwilightType(Enum):
    """
    Enumerated value to distinguish between the various sun angles relative to a horizon.
    The classifications for sun angle to twilight type is:
    Sun < -18 degrees --> Night time
    Sun < -12 degrees --> Astronomical twilight
    Sun < -6 degrees --> Nautical twilight
    Sun < -50 arc-minutes --> Civil twilight
    Otherwise --> Day time
    Most object can't be seen until their position is in civil twilight.
    """

    Day = 0
    Civil = 1
    Nautical = 2
    Astronomical = 3
    Night = 4

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented


def getSunPosition(time: JulianDate) -> EVector:
    """
    Very simple algorithm for computing the position of the Sun in a geocentric equitorial reference frame.

    Args:
        time: Time to find the Sun's position.

    Returns:
        The Sun's position vector in kilometers.
    """

    #   time since noon TT on Jan 1, 2000
    n = time.difference(J2000)
    #   mean longitude of the Sun
    L = 4.89495042 + 0.0172027924 * n
    #   mean anomaly of the Sun
    g = 6.240040768 + 0.0172019703 * n
    #   ecliptic longitude of the Sun
    lmbda = L + 0.0334230552 * sin(g) + 0.00034906585 * sin(2 * g)
    #   distance of the Sun from the Earth
    R = 1.00014 - 0.01671 * cos(g) - 0.00014 * cos(2 * g)
    #   obliquity of the ecliptic (axial tilt)
    e = 0.4090877234 - 6.981317008e-9 * n

    return EVector(
        R * AU * cos(lmbda),
        R * AU * cos(e) * sin(lmbda),
        R * AU * sin(e) * sin(lmbda)
    )


def getSunPosition2(time: JulianDate):
    """
    A more in depth algorithm for finding the Sun's position vector.

    Args:
        time: Time to find the Sun's position.

    Returns:
        The Sun's position vector in kilometers.
    """

    JDE = time.future(DELTAT / 86400)
    # JC = time.difference(J2000) / 36525.0
    JCE = JDE.difference(J2000) / 36525.0
    JME = JCE / 10.0

    L = expandTableValues(LTable, JME) % TWOPI
    B = expandTableValues(BTable, JME) % TWOPI
    R = expandTableValues(RTable, JME)

    theta = (L + pi) % TWOPI
    beta = -B
    dPsi, dEpsilon = getNutationDeltas(JCE)
    epsilon0 = 0.0
    for i in range(11):
        epsilon0 += UTable[i] * ((JME/10.0) ** i)
    epsilon = radians(epsilon0 / 3600.0) + dEpsilon
    dTau = -radians(20.4898 / (3600 * R))
    lmbda = theta + dPsi + dTau
    alpha = atan3(sin(lmbda) * cos(epsilon) - tan(beta) * sin(epsilon), cos(lmbda))
    delta = asin(sin(beta) * cos(epsilon) + cos(beta) * sin(epsilon) * sin(lmbda))
    print('alpha:', degrees(alpha))
    print('delta:', degrees(delta))

    zComp = sin(delta)
    xComp = sqrt((cos(delta) ** 2) / (1 + (tan(alpha) ** 2)))
    if pi / 2 < alpha < 3*pi/2:
        xComp = -abs(xComp)
    else:
        xComp = abs(xComp)

    yComp = xComp * tan(alpha)
    if alpha < pi:
        yComp = abs(yComp)
    else:
        yComp = -abs(yComp)

    '''
    This code works, I only need the right ascension and declination for my calculations for now.
    phi = radians(geo.getLatitude())
    nu0Raw = 280.46061837 + 360.98564736629 * time.difference(J2000) + 0.000387933 * JC*JC - (JC*JC*JC / 38710000)
    nu0 = radians(nu0Raw) % TWOPI
    nu = nu0 + dPsi * cos(epsilon)
    H = (nu + radians(geo.getLongitude()) - alpha) % TWOPI
    ksi = radians(8.794 / (3600 * R))
    u = atan(0.99664719 * tan(phi))
    x = cos(u) + geo.getElevation() * cos(phi) / 6378140
    y = 0.99664719 * sin(u) + geo.getElevation() * sin(phi) / 6378140
    dAlpha = atan2(-x * sin(ksi) * sin(H), cos(delta) - x * sin(ksi) * cos(H))
    alphaP = alpha + dAlpha
    deltaP = atan2((sin(delta) - y*sin(ksi))*cos(dAlpha), cos(delta) - x*sin(ksi)*cos(H))
    HP = H - dAlpha
    e0 = asin(sin(phi)*sin(deltaP) + cos(phi)*cos(deltaP)*cos(HP))
    # put atmospheric refraction correction here if temp and pressure daa is known
    e = e0
    zenith = pi/2 - e
    azimuth = (atan2(sin(HP), cos(HP)*sin(phi) - tan(deltaP)*cos(phi)) + pi) % TWOPI
    return azimuth, zenith'''
    return EVector(xComp, yComp, zComp)


def expandTableValues(table, JME):
    """returns in radians"""
    arr = [tableSum(ti, JME) for ti in table]
    sums = 0.0
    for i in range(len(table)):
        sums += arr[i] * (JME ** i)
    return sums / 1e8


def tableSum(table, JME):
    """returns in radians"""
    sums = 0.0
    for ti in table:
        sub = ti[0] * cos(ti[1] + ti[2] * JME)
        sums += sub
    return sums


def getNutationDeltas(JCE: float):
    """returns in radians"""
    arrX = getXArray(JCE)
    '''arrDPsi, arrDEpsilon = [], []
    for row in NutationTable:
        sumTerm = 0.0:
        for i in range(5):
            sumTerm += arrX[i] * row[i]
    '''

    dPsi, dEpsilon = 0.0, 0.0
    for ti in NutationTable:
        sumTerm = 0.0
        for i in range(5):
            sumTerm += arrX[i] * ti[i]
        dPsi += (ti[5] + ti[6] * JCE) * sin(sumTerm)
        dEpsilon += (ti[7] + ti[8] * JCE) * cos(sumTerm)
    return radians(dPsi / 36000000.0), radians(dEpsilon / 36000000.0)


def getXArray(JCE: float):
    """returns in radians"""
    X = []
    for i in range(5):
        Xi = 0.0
        for j in range(4):
            Xi += XTable[i][j] * (JCE ** j)
        X.append(radians(Xi))
    return X


LTable = [
    [
        [175347046, 0,          0],
        [3341656,   4.6692568,  6283.07585],
        [34894,     4.6261,     12566.1517],
        [3497,      2.7441,     5753.3849],
        [3418,      2.8289,     3.5231],
        [3136,      3.6277,     77713.7715],
        [2676,      4.4181,     7860.4194],
        [2343,      6.1352,     3930.2097],
        [1324,      0.7425,     11506.7698],
        [1273,      2.0371,     529.691],
        [1199,      1.1096,     1577.3435],
        [990,       5.233,      5884.927],
        [902,       2.045,      26.298],
        [857,       3.508,      398.149],
        [780,       1.179,      5223.694],
        [753,       2.533,      5507.553],
        [505,       4.583,      18849.228],
        [492,       4.205,      775.523],
        [357,       2.92,       0.067],
        [317,       5.849,      11790.629],
        [284,       1.899,      796.298],
        [271,       0.315,      10977.079],
        [243,       0.345,      5486.778],
        [206,       4.806,      2544.314],
        [205,       1.869,      5573.143],
        [202,       2.458,      6069.777],
        [156,       0.833,      213.299],
        [132,       3.411,      2942.463],
        [126,       1.083,      20.775],
        [115,       0.645,      0.98],
        [103,       0.636,      4694.003],
        [102,       0.976,      15720.839],
        [102,       4.267,      7.114],
        [99,        6.21,       2146.17],
        [98,        0.68,       155.42],
        [86,        5.98,       161000.69],
        [85,        1.3,        6275.96],
        [85,        3.67,       71430.7],
        [80,        1.81,       17260.15],
        [79,        3.04,       12036.46],
        [75,        1.76,       5088.63],
        [74,        3.5,        3154.69],
        [74,        4.68,       801.82],
        [70,        0.83,       9437.76],
        [62,        3.98,       8827.39],
        [61,        1.82,       7084.9],
        [57,        2.78,       6286.6],
        [56,        4.39,       14143.5],
        [56,        3.47,       6279.55],
        [52,        0.19,       12139.55],
        [52,        1.33,       1748.02],
        [51,        0.28,       5856.48],
        [49,        0.49,       1194.45],
        [41,        5.37,       8429.24],
        [41,        2.4,        19651.05],
        [39,        6.17,       10447.39],
        [37,        6.04,       10213.29],
        [37,        2.57,       1059.38],
        [36,        1.71,       2352.87],
        [36,        1.78,       6812.77],
        [33,        0.59,       17789.85],
        [30,        0.44,       83996.85],
        [30,        2.74,       1349.87],
        [25,        3.16,       4690.48]
    ],
    [
        [628331966747.0, 0,         0],
        [206059,        2.678235,   6283.07585],
        [4303,          2.6351,     12566.1517],
        [425,           1.59,       3.523],
        [119,           5.769,      26.298],
        [109,           2.966,      1577.344],
        [93,            2.59,       18849.23],
        [72,            1.14,       529.69],
        [68,            1.87,       398.15],
        [67,            4.41,       5507.55],
        [59,            2.89,       5223.69],
        [56,            2.17,       155.42],
        [45,            0.4,        796.3],
        [36,            0.47,       775.52],
        [29,            2.65,       7.11],
        [21,            5.34,       0.98],
        [19,            1.85,       5486.78],
        [19,            4.97,       213.3],
        [17,            2.99,       6275.96],
        [16,            0.03,       2544.31],
        [16,            1.43,       2146.17],
        [15,            1.21,       10977.08],
        [12,            2.83,       1748.02],
        [12,            3.26,       5088.63],
        [12,            5.27,       1194.45],
        [12,            2.08,       4694],
        [11,            0.77,       553.57],
        [10,            1.3,        6286.6],
        [10,            4.24,       1349.87],
        [9,             2.7,        242.73],
        [9,             5.64,       951.72],
        [8,             5.3,        2352.87],
        [6,             2.65,       9437.76],
        [6,             4.67,       4690.48]
    ],
    [
        [52919, 0,      0],
        [8720,  1.0721, 6283.0758],
        [309,   0.867,  12566.152],
        [27,    0.05,   3.52],
        [16,    5.19,   26.3],
        [16,    3.68,   155.42],
        [10,    0.76,   18849.23],
        [9,     2.06,   77713.77],
        [7,     0.83,   775.52],
        [5,     4.66,   1577.34],
        [4,     1.03,   7.11],
        [4,     3.44,   5573.14],
        [3,     5.14,   796.3],
        [3,     6.05,   5507.55],
        [3,     1.19,   242.73],
        [3,     6.12,   529.69],
        [3,     0.31,   398.15],
        [3,     2.28,   553.57],
        [2,     4.38,   5223.69],
        [2,     3.75,   0.98]
    ],
    [
        [289,   5.844,  6283.076],
        [35,    0,      0],
        [17,    5.49,   12566.15],
        [3,     5.2,    155.42],
        [1,     4.72,   3.52],
        [1,     5.3,    18849.23],
        [1,     5.97,   242.73]
    ],
    [
        [114,   3.142,  0],
        [8,     4.13,   6283.08],
        [1,     3.84,   12566.15]
    ],
    [
        [1, 3.14, 0]
    ]
]

BTable = [
    [
        [280,   3.199,  84334.662],
        [102,   5.422,  5507.553],
        [80,    3.88,   5223.69],
        [44,    3.7,    2352.87],
        [32,    4,      1577.34]
    ],
    [
        [9, 3.9,    5507.55],
        [6, 1.73,   5223.69]
    ]
]

RTable = [
    [
        [100013989, 0,          0],
        [1670700,   3.0984635,  6283.07585],
        [13956,     3.05525,    12566.1517],
        [3084,      5.1985,     77713.7715],
        [1628,      1.1739,     5753.3849],
        [1576,      2.8469,     7860.4194],
        [925,       5.453,      11506.77],
        [542,       4.564,      3930.21],
        [472,       3.661,      5884.927],
        [346,       0.964,      5507.553],
        [329,       5.9,        5223.694],
        [307,       0.299,      5573.143],
        [243,       4.273,      11790.629],
        [212,       5.847,      1577.344],
        [186,       5.022,      10977.079],
        [175,       3.012,      18849.228],
        [110,       5.055,      5486.778],
        [98,        0.89,       6069.78],
        [86,        5.69,       15720.84],
        [86,        1.27,       161000.69],
        [65,        0.27,       7260.15],
        [63,        0.92,       529.69],
        [57,        2.01,       83996.85],
        [56,        5.24,       71430.7],
        [49,        3.25,       2544.31],
        [47,        2.58,       775.52],
        [45,        5.54,       9437.76],
        [43,        6.01,       6275.96],
        [39,        5.36,       4694],
        [38,        2.39,       8827.39],
        [37,        0.83,       19651.05],
        [37,        4.9,        12139.55],
        [36,        1.67,       12036.46],
        [35,        1.84,       2942.46],
        [33,        0.24,       7084.9],
        [32,        0.18,       5088.63],
        [32,        1.78,       398.15],
        [28,        1.21,       6286.6],
        [28,        1.9,        6279.55],
        [26,        4.59,       10447.39]
    ],
    [
        [103019, 1.10749, 6283.07585],
        [1721, 1.0644, 12566.1517],
        [702, 3.142, 0],
        [32, 1.02, 18849.23],
        [31, 2.84, 5507.55],
        [25, 1.32, 5223.69],
        [18, 1.42, 1577.34],
        [10, 5.91, 10977.08],
        [9, 1.42, 6275.96],
        [9, 0.27, 5486.78]
    ],
    [
        [4359,  5.7846, 6283.0758],
        [124,   5.579,  12566.152],
        [12,    3.14,   0],
        [9,     3.63,   77713.77],
        [6,     1.87,   5573.14],
        [3,     5.47,   18849.23]
    ],
    [
        [145,   4.273,  6283.076],
        [7,     3.92,   12566.15]
    ],
    [
        [4, 2.56, 6283.08]
    ]
]

XTable = [
    [297.85036, 445267.111480, -0.0019142, 1.0/189474.0],
    [357.527772, 35999.050340, -0.0001603, 1.0/-300000.0],
    [134.96298, 477198.867398, 0.0086972, 1.0/56250.0],
    [93.27191, 483202.017538, -0.0036825, 1.0/327270.0],
    [125.04452, -1934.136261, 0.0020708, 1.0/450000.0]
]

NutationTable = [
    [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
    [-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1],
    [0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5],
    [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
    [0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1],
    [0, 0, 1, 0, 0, 712, 0.1, -7, 0],
    [-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6],
    [0, 0, 0, 2, 1, -386, -0.4, 200, 0],
    [0, 0, 1, 2, 2, -301, 0, 129, -0.1],
    [-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3],
    [-2, 0, 1, 0, 0, -158, 0, 0, 0],
    [-2, 0, 0, 2, 1, 129, 0.1, -70, 0],
    [0, 0, -1, 2, 2, 123, 0, -53, 0],
    [2, 0, 0, 0, 0, 63, 0, 0, 0],
    [0, 0, 1, 0, 1, 63, 0.1, -33, 0],
    [2, 0, -1, 2, 2, -59, 0, 26, 0],
    [0, 0, -1, 0, 1, -58, -0.1, 32, 0],
    [0, 0, 1, 2, 1, -51, 0, 27, 0],
    [-2, 0, 2, 0, 0, 48, 0, 0, 0],
    [0, 0, -2, 2, 1, 46, 0, -24, 0],
    [2, 0, 0, 2, 2, -38, 0, 16, 0],
    [0, 0, 2, 2, 2, -31, 0, 13, 0],
    [0, 0, 2, 0, 0, 29, 0, 0, 0],
    [-2, 0, 1, 2, 2, 29, 0, -12, 0],
    [0, 0, 0, 2, 0, 26, 0, 0, 0],
    [-2, 0, 0, 2, 0, -22, 0, 0, 0],
    [0, 0, -1, 2, 1, 21, 0, -10, 0],
    [0, 2, 0, 0, 0, 17, -0.1, 0, 0],
    [2, 0, -1, 0, 1, 16, 0, -8, 0],
    [-2, 2, 0, 2, 2, -16, 0.1, 7, 0],
    [0, 1, 0, 0, 1, -15, 0, 9, 0],
    [-2, 0, 1, 0, 1, -13, 0, 7, 0],
    [0, -1, 0, 0, 1, -12, 0, 6, 0],
    [0, 0, 2, -2, 0, 11, 0, 0, 0],
    [2, 0, -1, 2, 1, -10, 0, 5, 0],
    [2, 0, 1, 2, 2, -8, 0, 3, 0],
    [0, 1, 0, 2, 2, 7, 0, -3, 0],
    [-2, 1, 1, 0, 0, -7, 0, 0, 0],
    [0, -1, 0, 2, 2, -7, 0, 3, 0],
    [2, 0, 0, 2, 1, -7, 0, 3, 0],
    [2, 0, 1, 0, 0, 6, 0, 0, 0],
    [-2, 0, 2, 2, 2, 6, 0, -3, 0],
    [-2, 0, 1, 2, 1, 6, 0, -3, 0],
    [2, 0, -2, 0, 1, -6, 0, 3, 0],
    [2, 0, 0, 0, 1, -6, 0, 3, 0],
    [0, -1, 1, 0, 0, 5, 0, 0, 0],
    [-2, -1, 0, 2, 1, -5, 0, 3, 0],
    [-2, 0, 0, 0, 1, -5, 0, 3, 0],
    [0, 0, 2, 2, 1, -5, 0, 3, 0],
    [-2, 0, 2, 0, 1, 4, 0, 0, 0],
    [-2, 1, 0, 2, 1, 4, 0, 0, 0],
    [0, 0, 1, -2, 0, 4, 0, 0, 0],
    [-1, 0, 1, 0, 0, -4, 0, 0, 0],
    [-2, 1, 0, 0, 0, -4, 0, 0, 0],
    [1, 0, 0, 0, 0, -4, 0, 0, 0],
    [0, 0, 1, 2, 0, 3, 0, 0, 0],
    [0, 0, -2, 2, 2, -3, 0, 0, 0],
    [-1, -1, 1, 0, 0, -3, 0, 0, 0],
    [0, 1, 1, 0, 0, -3, 0, 0, 0],
    [0, -1, 1, 2, 2, -3, 0, 0, 0],
    [2, -1, -1, 2, 2, -3, 0, 0, 0],
    [0, 0, 3, 2, 2, -3, 0, 0, 0],
    [2, -1, 0, 2, 2, -3, 0, 0, 0]
]

UTable = [
    84381.448, -4680.93, -1.55, 1999.25,
    51.38, -249.67, -39.05, 7.12,
    27.87, 5.79, 2.45
]
