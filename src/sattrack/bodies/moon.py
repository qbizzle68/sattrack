from math import cos, sin, radians, tan, asin, degrees

from pyevspace import Vector

from sattrack.bodies._tables import MOON_LR_TERM_TABLE, MOON_B_TERM_TABLE
from sattrack.bodies.position import computeNutationDeltas, computeTrueObliquity
from sattrack.core.juliandate import JulianDate
from sattrack.core.coordinates import CelestialCoordinates
from sattrack.util.constants import DELTAT, TWOPI, RAD_TO_HOURS
from sattrack.util.helpers import atan3


def _computeMoonMeanLongitude(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    return 3.8103408236230014 + JCE * \
        (8399.709111633996 + JCE *
         (-2.7551767571982487e-05 + JCE *
          (3.239043153721282e-08 + JCE * -2.67713171763403e-10)))


def _computeMeanElongationMoon(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    return 5.198466529842604 + JCE * \
        (7771.377144833719 + JCE *
         (-3.2845351193281285e-05 + JCE *
          (3.197346706519396e-08 + JCE * -1.5436512200896206e-10)))


def _computeSunMeanAnomaly(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    return 6.24006012726235 + JCE * \
        (628.3019551672274 + JCE *
         (-2.68082573106329e-06 + JCE * 7.126701723129153e-10))


def _computeMoonMeanAnomaly(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    return 2.3555556368542616 + JCE * \
        (8328.691424759154 + JCE *
         (0.00015256621123383232 + JCE *
          (2.504095111829911e-07 + JCE * -1.1863303779189298e-09)))


def _computeMoonArgumentLatitude(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    return 1.6279051579829402 + JCE * \
        (8433.46615806092 + JCE *
         (-6.37725855386208e-05 + JCE *
          (-4.949884435605018e-09 + JCE * 2.0216715339731145e-11)))


def _computeLRTerms(time: 'JulianDate') -> (float, float):
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525
    E = 1 + JCE * (-0.002516 + JCE * -0.0000074)

    D = _computeMeanElongationMoon(time)
    M = _computeSunMeanAnomaly(time)
    MPrime = _computeMoonMeanAnomaly(time)
    F = _computeMoonArgumentLatitude(time)

    lTerm = 0
    rTerm = 0
    for d, m, mPrime, f, l, r in MOON_LR_TERM_TABLE:
        arg = d * D + m * M + mPrime * MPrime + f * F
        modifier = 1
        if m == 1 or m == -1:
            modifier = E
        elif m == 2 or m == -2:
            modifier = E * E

        lTerm += (l * modifier * sin(arg))
        rTerm += (r * modifier * cos(arg))

    return radians(lTerm / 1000000), rTerm / 1000


def _computeBTerm(time: 'JulianDate') -> float:
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525
    E = 1 + JCE * (-0.002516 + JCE * -0.0000074)

    D = _computeMeanElongationMoon(time)
    M = _computeSunMeanAnomaly(time)
    MPrime = _computeMoonMeanAnomaly(time)
    F = _computeMoonArgumentLatitude(time)

    bTerm = 0
    for d, m, mPrime, f, b in MOON_B_TERM_TABLE:
        arg = d * D + m * M + mPrime * MPrime + f * F
        modifier = 1
        if m == 1 or m == -1:
            modifier = E
        elif m == 2 or m == -2:
            modifier = E * E

        bTerm += b * modifier * sin(arg)

    return radians(bTerm / 1000000)


def _computeCoordinateDeltas(time: 'JulianDate') -> (float, float, float):
    JDE = time.value + DELTAT / 86400
    JCE = (JDE - 2451545) / 36525

    a1 = 2.0900317792632097 + 2.301199165462003 * JCE
    a2 = 0.9265952998837896 + 8364.739847732933 * JCE
    a3 = 5.470734540376226 + 8399.684725296609 * JCE

    moonLongitude = _computeMoonMeanLongitude(time)
    moonLatitude = _computeMoonArgumentLatitude(time)
    moonMeanAnomaly = _computeMoonMeanAnomaly(time)

    deltaL = 3958 * sin(a1) + 1962 * sin(moonLongitude - moonLatitude) + 318 * sin(a2)
    deltaB = -2235 * sin(moonLongitude) + 382 * sin(a3) + 175 * sin(a1 - moonLatitude) \
        + 175 * sin(a1 + moonLatitude) + 127 * sin(moonLongitude - moonMeanAnomaly) \
        - 115 * sin(moonLongitude + moonMeanAnomaly)

    return radians(deltaL / 1000000), radians(deltaB / 1000000)


def computeMoonLongitude(time: 'JulianDate') -> float:
    longitude = _computeMoonMeanLongitude(time)
    l, _ = _computeLRTerms(time)
    deltaL, _ = _computeCoordinateDeltas(time)

    return (longitude + l + deltaL) % TWOPI


def computeMoonLatitude(time: 'JulianDate') -> float:
    b = _computeBTerm(time)
    _, deltaB = _computeCoordinateDeltas(time)

    return (b + deltaB) % TWOPI


def computeMoonDistance(time: 'JulianDate') -> float:
    _, r = _computeLRTerms(time)

    return 385000.56 + r


def computeMoonApparentLongitude(time: 'JulianDate') -> float:
    longitude = computeMoonLongitude(time)
    deltaPsi, _ = computeNutationDeltas(time)

    return (longitude + deltaPsi) % TWOPI


def computeMoonRightAscension(time: 'JulianDate') -> float:
    apparentLongitude = computeMoonApparentLongitude(time)
    obliquity = computeTrueObliquity(time)
    latitude = computeMoonLatitude(time)

    numerator = sin(apparentLongitude) * cos(obliquity) - tan(latitude) * sin(obliquity)

    return atan3(numerator, cos(apparentLongitude))


def computeMoonDeclination(time: 'JulianDate') -> float:
    apparentLongitude = computeMoonApparentLongitude(time)
    obliquity = computeTrueObliquity(time)
    latitude = computeMoonLatitude(time)

    term1 = sin(latitude) * cos(obliquity)
    term2 = cos(latitude) * sin(obliquity) * sin(apparentLongitude)

    return asin(term1 + term2)


def _computeMoonCoordinatesFast(time: 'JulianDate') -> (float, float):
    apparentLongitude = computeMoonApparentLongitude(time)
    obliquity = computeTrueObliquity(time)
    latitude = computeMoonLatitude(time)

    numerator = sin(apparentLongitude) * cos(obliquity) - tan(latitude) * sin(obliquity)
    rightAscension = atan3(numerator, cos(apparentLongitude))

    term1 = sin(latitude) * cos(obliquity)
    term2 = cos(latitude) * sin(obliquity) * sin(apparentLongitude)
    declination = asin(term1 + term2)

    return rightAscension, declination


def computeMoonPosition(time: 'JulianDate') -> Vector:
    rightAscension, declination = _computeMoonCoordinatesFast(time)
    moonDistance = computeMoonDistance(time)

    zTerm = moonDistance * sin(declination)
    topoProjection = moonDistance * cos(declination)
    xTerm = topoProjection * cos(rightAscension)
    yTerm = topoProjection * sin(rightAscension)

    return Vector(xTerm, yTerm, zTerm)


def computeMoonCoordinates(time: 'JulianDate') -> CelestialCoordinates:
    rightAscension, declination = _computeMoonCoordinatesFast(time)

    raHours = rightAscension * RAD_TO_HOURS
    decDegrees = degrees(declination)

    return CelestialCoordinates(raHours, decDegrees)
