from math import acos
from operator import index
from typing import Union

from numpy import arange
from pyevspace import norm, dot
import matplotlib.pyplot as plt

from sattrack.structures.orbit import _mean_to_true_anomaly_newton

from sattrack.api import *
from sattrack.structures.orbit import *
import sattrack.structures.orbit as orbit

tle = getTle('zarya')
# tle = TwoLineElement("""ISS (ZARYA)
# 1 25544U 98067A   22266.84431519  .00008111  00000+0  14870-3 0  9996
# 2 25544  51.6423 207.8056 0002412 286.8120 181.5821 15.50238875360488""")
# tle = TwoLineElement(tleStr)
iss = Satellite(tle)
from sattrack.util.anomalies import timeToNextTrueAnomaly

jd = now()
geo = GeoPosition(38.0608, -97.9298)
# pc = PassController(iss, geo, jd)
# np = pc.getNextPass()
# plist = pc.getPassList(2)
# elements = OrbitalElements.fromTle(tle, jd)

def __compute_s_vector(sunPosition: EVector, raan: float, inclination: float, aop: float) -> EVector:
    S = -norm(sunPosition)
    return rotateOrderTo(ZXZ, EulerAngles(raan, inclination, aop), S)

def __compute_gamma(sVector: EVector) -> float:
    return atan2(sVector[1], -sVector[0])

Shadow = Union[int]
UMBRA = Shadow(0)
PENUMBRA = Shadow(1)

# todo: make a better name for this
EnterExit = Union[int]
ENTER = EnterExit(0)
EXIT = EnterExit(1)

def __escobal_method(R: float, sVector: EVector, phi: float, zeta: float, sma: float, ecc: float, shadow: Shadow) -> float:
    eTerm = 1 + ecc * cos(phi)
    aTerm = sma * (1 - ecc * ecc)
    sTerm = -sVector[0] * cos(phi) - sVector[1] * sin(phi)

    term1 = R * R * eTerm * eTerm
    term2 = aTerm * aTerm * sTerm * sTerm
    term3 = aTerm * aTerm * cos(zeta) * cos(zeta)
    term4 = 2 * aTerm * R * sTerm * eTerm * sin(zeta)

    if shadow is PENUMBRA:
        term4 *= -1

    return term1 + term2 - term3 + term4

def __escobal_method_derivative(R: float, sVector: EVector, phi: float, zeta: float, sma: float, ecc: float, shadow: Shadow) -> float:
    eTerm = 1 + ecc * cos(phi)
    aTerm = sma * (1 - ecc * ecc)
    sTerm = -sVector[0] * cos(phi) - sVector[1] * sin(phi)
    eTermPrime = -ecc * sin(phi)
    sTermPrime = sVector[0] * sin(phi) - sVector[1] * cos(phi)

    term1 = 2 * R * R * eTerm * eTermPrime
    term2 = 2 * aTerm * aTerm * sTerm * sTermPrime
    term4 = 2 * R * aTerm * sin(zeta) * (sTerm * eTermPrime + sTermPrime * eTerm)

    if shadow is PENUMBRA:
        term4 *= -1

    return term1 + term2 + term4

# todo: geosynchronous satellites that are in sunlight during opposition don't have zeros to this method, need to check
#  and make an exception for them
def __find_zero_newton(R, sVector, guess, zeta, sma, ecc, shadow, epsilon = 1e-5):
    phi = guess
    gi = __escobal_method(R, sVector, phi, zeta, sma, ecc, shadow)
    while abs(gi) > epsilon:
        phi = phi - (gi / __escobal_method_derivative(R, sVector, phi, zeta, sma, ecc, shadow))
        gi = __escobal_method(R, sVector, phi, zeta, sma, ecc, shadow)
    return phi % TWOPI

def __find_fast_zero(R, sVector, gamma, zeta, sma, ecc, shadow, enterOrExit, epsilon = 1e-5):
    if index(enterOrExit) is ENTER:
        return __find_zero_newton(R, sVector, 3*pi/4 - gamma, zeta, sma, ecc, shadow, epsilon)
    elif index(enterOrExit) is EXIT:
        return __find_zero_newton(R, sVector, 5*pi/4 - gamma, zeta, sma, ecc, shadow, epsilon)
    else:
        raise ValueError('enterOrExit parameter value must be ENTER or EXIT')

def __find_certain_zeros(R, sVector, zeta, sma, ecc, shadow, epsilon = 1e-5):
    zeros = []
    frac = 0.5
    while len(zeros) < 4:
        guesses = [i * pi * frac for i in range(int(2 / frac))]
        computed = [__find_zero_newton(R, sVector, guess, zeta, sma, ecc, shadow, epsilon) for guess in guesses]
        # round to the place of an epsilon value to compare equivalent floating values
        rounded = {round(c*1e5)/1e5 for c in computed}
        zeros.clear()
        for phi in computed:
            roundedPhi = round(phi*1e5)/1e5
            if roundedPhi in rounded:
                zeros.append(phi)
                rounded.remove(roundedPhi)
        frac /= 2
    return zeros

def __check_zero(phi, sVector):
    return (sVector[0] * cos(phi) + sVector[1] * sin(phi)) > 0

def __check_range(psi, enterOrExit):
    return (enterOrExit is ENTER and pi/2 <= psi <= pi) or (enterOrExit is EXIT and pi <= psi <= 3*pi/2)

def _get_zero(R, sVector, zeta, sma, ecc, shadow, enterOrExit, *, guess=None, epsilon=1e-5):
    # if guess is set use it first, then try __find_fast_zero, then use __find_certain_zeros
    gamma = __compute_gamma(sVector)
    if guess is not None:
        phi = __find_zero_newton(R, sVector, guess, zeta, sma, ecc, shadow)
        if __check_zero(phi, sVector):
            psi = (gamma + phi) % TWOPI
            # if (enterOrExit is ENTER and pi/2 <= psi <= pi) or (enterOrExit is EXIT and pi >= psi >= 3*pi/2):
            if __check_range(psi, enterOrExit):
                return phi
    phi = __find_fast_zero(R, sVector, gamma, zeta, sma, ecc, shadow, enterOrExit, epsilon)
    if __check_zero(phi, sVector):
        psi = (gamma + phi) % TWOPI
        # if (enterOrExit is ENTER and phi <= pi) or (enterOrExit is EXIT and phi >= pi):
        if __check_range(psi, enterOrExit):
            return phi
    zeros = __find_certain_zeros(R, sVector, zeta, sma, ecc, shadow, epsilon)
    for phi in zeros:
        if __check_zero(phi, sVector):
            psi = (gamma + phi) % TWOPI
            # if (enterOrExit is ENTER and phi <= pi) or (enterOrExit is EXIT and phi >= pi):
            if __check_range(psi, enterOrExit):
                return phi
    # todo: make this a custom exception type
    raise Exception('unable to find a valid zero')

def __get_radius_z_comp(sVector, sma, ecc, inc, aop, phi):
    term1 = (sma * (1 - ecc*ecc)) / (1 + ecc * cos(phi))
    term2 = sin(aop) * sin(inc) * cos(phi)
    term3 = cos(aop) * sin(inc) * sin(phi)
    term4 = sVector[2] * (sVector[0] * cos(phi) + sVector[1] * sin(phi))
    return term1 * (term2 + term3 - term4)

def __get_latitude_term(Rz):
    ae = EARTH_EQUITORIAL_RADIUS
    # f = EARTH_FLATTENING
    # fTerm = -f * (2 -f)
    # -f * (2 - f) where f is EARTH_FLATTENING
    fTerm = -0.006694317778266723
    aeTerm = ae * ae * (1 + fTerm)
    return (aeTerm - Rz * Rz) / (aeTerm + (Rz * Rz * fTerm))

def __get_radius_from_latitude(latitudeTerm):
    ae = EARTH_EQUITORIAL_RADIUS
    # f = EARTH_FLATTENING
    # fTerm = f * (2 - f)
    # f * (2 - f) where f is EARTH_FLATTENING
    fTerm = 0.006694317778266723
    return (ae * sqrt(1 - fTerm)) / sqrt(1 - fTerm * latitudeTerm)

def __get_aperture_angle(rs, Re, shadow):
    if shadow is UMBRA:
        cosZeta = sqrt(rs * rs - (SUN_RADIUS - Re) ** 2) / rs
    elif shadow is PENUMBRA:
        cosZeta = sqrt(rs * rs - (SUN_RADIUS + Re) ** 2) / rs
    else:
        raise ValueError('shadow parameter must be either UMBRA or PENUMBRA')
    return acos(cosZeta)

def __get_refraction_angle(altitudeAngle):
    numerator = 0.009928887226387075 + altitudeAngle * (0.06995 + altitudeAngle * 0.004087098938599872)
    denominator = 1 + altitudeAngle * (28.934368654106574 + altitudeAngle * 277.39713657599236)
    return numerator / denominator

def __get_corrected_refraction_angle(rs, Re, shadow):
    semiApertureAngle = __get_aperture_angle(rs, Re, shadow)
    refractionAngle = __get_refraction_angle(semiApertureAngle)
    if shadow is UMBRA:
        correctedAngle = semiApertureAngle + refractionAngle
    elif shadow is PENUMBRA:
        correctedAngle = refractionAngle - semiApertureAngle
    else:
        raise ValueError('shadow parameter must be either UMBRA or PENUMBRA')
    return correctedAngle

# todo: make the sat parameter an object that an OrbitalElements object can be derived from
# idea behind the interface, use an orbitable object as the parameter instead of a TLE, so we can
# compute values from objects not constructed by a tle
def _get_shadow_anomaly(jd: JulianDate, sat: Satellite, shadow: Shadow, zeroEpsilon=1e-5,
                        radiusEpsilon=1e-5) -> (float, JulianDate):
    Re = 6371 # average radius of the earth
    tle = sat.getTle()
    time = jd
    sunPosition = getSunPosition(time)
    elements = OrbitalElements.fromTle(tle, time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
    apertureAngle = __get_corrected_refraction_angle(sunPosition.mag(), Re, shadow)

    phi0 = elements.trueAnomalyAt(jd)
    approxPhi1 = _get_zero(Re, sVector, apertureAngle, elements.sma, elements.ecc, shadow, ENTER, epsilon=zeroEpsilon)
    approxPhi2 = _get_zero(Re, sVector, apertureAngle, elements.sma, elements.ecc, shadow, EXIT, epsilon=zeroEpsilon)

    # maximum difference between the ranges of Re seems to be about 0.015, so take a very conservative value of 0.1 rad
    errorBuffer = 0.1
    # were close enough to the exit we don't know if we're before or after it with errors
    if (abs(phi0 - approxPhi2) < errorBuffer) or (abs(phi0 + TWOPI - approxPhi2) < errorBuffer) \
            or (abs(approxPhi2 + TWOPI - phi0) < errorBuffer):
        dt = (jd - elements.timeToPrevTrueAnomaly(approxPhi1, jd)) / 2
        referenceTime = jd.future(-dt)
    else:
        phi2Time = elements.timeToNextTrueAnomaly(approxPhi2, jd)
        phi1Time = elements.timeToPrevTrueAnomaly(approxPhi1, phi2Time)
        dt = (phi2Time - phi1Time) / 2
        referenceTime = phi1Time.future(dt)

    enterPhi, enterTime = __compute_anomaly_loop(jd, referenceTime, sat, shadow, ENTER, zeroEpsilon, radiusEpsilon)
    exitPhi, exitTime = __compute_anomaly_loop(jd, referenceTime, sat, shadow, EXIT, zeroEpsilon, radiusEpsilon)

    if enterTime < exitTime < jd:
        gamma = __compute_gamma(sVector)
        # todo: make sure this is at conjunction
        updatedJd = elements.timeToNextTrueAnomaly(gamma, jd)
        return _get_shadow_anomaly(updatedJd, sat, shadow, zeroEpsilon)

    return (enterPhi, enterTime), (exitPhi, exitTime)

def __compute_anomaly_loop(startTime, referenceTime, sat, shadow, enterOrExit, zeroEpsilon=1e-5,
                           radiusEpsilon=1e-5):
    Re = 6371
    tle = sat.getTle()
    time = startTime
    sunPosition = getSunPosition(time)
    elements = OrbitalElements.fromTle(tle, time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)

    # need a do-while structure here
    while True:
        apertureAngle = __get_corrected_refraction_angle(sunPosition.mag(), Re, shadow)
        phi = _get_zero(Re, sVector, apertureAngle, elements.sma, elements.ecc, shadow, enterOrExit, epsilon=zeroEpsilon)
        if enterOrExit is ENTER:
            time = elements.timeToPrevTrueAnomaly(phi, referenceTime)
        elif enterOrExit is EXIT:
            time = elements.timeToNextTrueAnomaly(phi, referenceTime)
        elements = OrbitalElements.fromTle(tle, time)
        sunPosition = getSunPosition(time)
        sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
        Rz = __get_radius_z_comp(sVector, elements.sma, elements.ecc, elements.inc, elements.aop, phi)
        latitudeTerm = __get_latitude_term(Rz)
        previousRe = Re
        Re = __get_radius_from_latitude(latitudeTerm)
        if (abs(Re - previousRe) > radiusEpsilon):
            break

    return phi, time

if __name__ == "__main__":
    starlink_numbers = [
        1007, 1011, 1022, 1029,
        1073, 1071, 1110, 1085,
        1132, 1141, 1173, 1138,
        1201, 1234, 1279, 1301,
        1306, 1294, 1346, 1332,
        1441, 1393, 1461, 1475,
        1391, 1523, 1524, 1515,
        1585, 1588, 1654, 1686,
        1550, 1644, 1687, 1774
    ]
    R = 6371
    zeta = 0.013609786618843114
    passedCount = 0

    if False:
        printString = '{:^10}|{:^10}|{:^10}|{:^10}{:^32}'\
            .format('phi1P', 'phi2P', 'phi1U', 'phi2U', 'range results')
        print(printString)
        for i, starlink in enumerate(starlink_numbers):
            print(f'testing starlink-{starlink}', end='\r')
            tle = getTle(f'starlink-{starlink}')
            sat = Satellite(tle)
            jd = now()
            sunPos = getSunPosition(jd)
            elements = OrbitalElements.fromTle(tle, jd)
            sVector = __compute_s_vector(jd, elements.raan, elements.inc, elements.aop)
            phi1P = __find_zero_newton(R, sVector, 3*pi/4, zeta, elements.sma, elements.ecc, PENUMBRA)
            phi2P = __find_zero_newton(R, sVector, 5*pi/4, zeta, elements.sma, elements.ecc, PENUMBRA)
            phi1U = __find_zero_newton(R, sVector, 3*pi/4, zeta, elements.sma, elements.ecc, UMBRA)
            phi2U = __find_zero_newton(R, sVector, 5*pi/4, zeta, elements.sma, elements.ecc, UMBRA)
            phi1PResult = pi/2 <= phi1P <= pi
            phi1UResult = pi/2 <= phi1U <= pi
            phi2PResult = pi <= phi2P <= 3*pi/2
            phi2UResult = pi <= phi2U <= 3*pi/2
            printString = '{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^8}{:^8}{:^8}{:^8}' \
                .format(phi1P, phi2P, phi1U, phi2U, bool(phi1PResult), bool(phi2PResult), bool(phi1UResult),
                        bool(phi2UResult), w=10, p=6)
            print(printString)
            if phi1PResult and phi1UResult and phi2PResult and phi2UResult:
                passedCount += 1
        print(f"RESULTS: {passedCount}/36")

    if False:
        printString = '{:^10}|{:^10}|{:^10}|{:^10}|{:^10}|{:^40}' \
            .format('number', 'phi0', 'phi1', 'phi2', 'phi3', 'escobal method values')
        print(printString)
        for i, starlink in enumerate(starlink_numbers):
            print(f'testing starlink-{starlink}', end='\r')
            tle = getTle(f'starlink-{starlink}')
            sat = Satellite(tle)
            jd = now()
            sunPos = getSunPosition(jd)
            elements = OrbitalElements.fromTle(tle, jd)
            sVector = __compute_s_vector(jd, elements.raan, elements.inc, elements.aop)
            phis = __find_certain_zeros(R, sVector, zeta, elements.sma, elements.ecc, PENUMBRA)
            phis.sort()
            vals = [round(__escobal_method(R, sVector, phi, zeta, elements.sma, elements.ecc, PENUMBRA)*1e5)/1e5 for phi in phis]
            printString = '{:^10}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}'\
                .format(starlink, phis[0], phis[1], phis[2], phis[3], vals[0], vals[1], vals[2], vals[3], w=10, p=6)
            print(printString)
            if vals[0] <= 1e-5 and vals[1] <= 1e-5 and vals[2] <= 1e-5 and vals[3] <= 1e-5:
                passedCount += 1
        print(f'RESULTS: {passedCount}/36')

    if False:
        printString = '{:^10}|{:^10}|{:^10}|{:^10}|{:^10}'.format('number', 'phi1P', 'phi1U', 'phi2U', 'phi2P')
        printString2 = '{:-^10} {:-^10} {:-^10} {:-^10} {:-^10}'.format('','','','','')
        print(printString, printString2, sep='\n')
        for i, starlink in enumerate(starlink_numbers):
            print(f'testing starlink-{starlink}', end='\r')
            tle = getTle(f'starlink-{starlink}')
            sat = Satellite(tle)
            jd = now()
            sunPos = getSunPosition(jd)
            elements = OrbitalElements.fromTle(tle, jd)
            # todo: this interface changed
            sVector = __compute_s_vector(jd, elements.raan, elements.inc, elements.aop)
            phi1P = __get_zero(R, sVector, zeta, elements.sma, elements.ecc, PENUMBRA, ENTER)
            phi1U = __get_zero(R, sVector, zeta, elements.sma, elements.ecc, UMBRA, ENTER)
            phi2P = __get_zero(R, sVector, zeta, elements.sma, elements.ecc, PENUMBRA, EXIT)
            phi2U = __get_zero(R, sVector, zeta, elements.sma, elements.ecc, UMBRA, EXIT)
            printString = '{:^10}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}'\
                .format(starlink, phi1P, phi1U, phi2U, phi2P, w=10, p=6)
            print(printString)
        print("FINISHED")

    if False:
        def approxTA(M, ecc):
            return M + 2 * ecc * sin(M) + 1.25 * ecc * ecc * sin(2 * M)
        ecc = 0.017
        print('eccentricity: ', ecc)
        printString = '{:^10}|{:^10}|{:^10}|{:^10}'.format('M', 'direct v', 'approx v', 'error')
        printString2 = '{:-^10} {:-^10} {:-^10} {:-^10}'.format('', '', '', '')
        print(printString, printString2, sep='\n')
        ma = [i * 0.1 for i in range(62)]
        for m in ma:
            directV = _mean_to_true_anomaly_newton(m, ecc)
            approxV = approxTA(m, ecc)
            error = degrees(directV - approxV) * 60
            printString = '{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}|{:^{w}.{p}f}'.format(m, directV, approxV, error, w=10, p=6)
            print(printString)
        print("FINISHED")


    if True:
        def getConst(r, v):
            return r*v*v/EARTH_MU

        def getEcc(r, v, gamma):
            c = getConst(r, v)
            return sqrt(((c-1)**2) * (sin(gamma)**2) + (cos(gamma)**2))

        x = [i for i in range(1440)]

        # states = [iss.getState(jd.future(i/1440)) for i in x]
        # r = [s[0].mag() for s in states]
        # v = [s[1].mag() for s in states]
        # gammas = [vang(s[0], s[1]) for s in states]
        # y = [getEcc(r, v, gamma) for r, v, gamma in zip(r, v, gammas)]
        y = [Elements.fromTle(tle, jd.future(i/1440)).ecc for i in x]

        plt.plot(x, y)
        plt.show()

        # elements = Elements.fromTle(tle, jd)
        # pos, vel = iss.getState(jd)
        # gamma = vang(pos, vel)
        # r = pos.mag()
        # v = vel.mag()
        # c = (r*v*v/EARTH_MU)
        # ecc = sqrt((c-1)*(c-1)*sin(gamma)*sin(gamma) + cos(gamma)*cos(gamma))
        # ta = atan3(c * sin(gamma) * cos(gamma), c * sin(gamma)*sin(gamma) - 1)

        # print('from elements: e: ', elements.ecc, ' v: ', orbit._mean_to_true_anomaly(elements.meanAnomaly, elements.ecc))
        # print('by hand: e: ', ecc, ' v: ', ta)


class SatelliteShadow:

    __slots__ = '_phi', '_elements', '_s', '_Re', '_cosTerm', '_zeta', '_'

    def __init__(self, sat: Satellite, jd: JulianDate):
        pass


Rs = 6.957e8
Re = 6371000
rs = 149597870700

def getRefractionAngle(apparentAltitude: float):
    h = degrees(apparentAltitude)
    return (0.009928887226387075 + 0.0012208578117700335*h + 1.2450015330892884e-06*h*h) \
                      / (1 + 0.505*h + 0.0845*h*h)

def getUmbraAperture(rs, Rs, Re):
    zetaU = acos(sqrt(rs*rs - (Rs-Re)**2) / rs)
    deltaZetaU = getRefractionAngle(zetaU)
    return zetaU + deltaZetaU
    # sinZeta = (Rs - Re) / rs

def getPenumbraAperture(rs, Rs, Re):
    zetaP = acos(sqrt(rs*rs - (Rs+Re)**2) / rs)
    deltaZetaP = getRefractionAngle(zetaP)
    return deltaZetaP - zetaP
    # sinZeta = (Rs + Re) / rs

# THIS ONLY WORKS FOR PENUMBRA
def gPenumbra(Re, e, phi, a, sx, sy, z):
    eTerm = 1 + e * cos(phi)
    aTerm = a * (1 - e*e)
    sTerm = -sx*cos(phi) - sy*sin(phi)
    cosZeta = cos(z)

    term1 = Re * Re * eTerm * eTerm
    term2 = aTerm * aTerm * sTerm * sTerm
    term3 = aTerm * aTerm * cosZeta * cosZeta
    term4 = 2 * aTerm * Re * sTerm * eTerm * sin(z)

    return term1 + term2 - term3 - term4

def updateElements(sat, jd):
    elements = OrbitalElements.fromTle(sat.getTle(), jd)
    return elements.ecc, elements.sma, elements.raan, elements.inc, elements.aop

def getSVector(jd, raan, inc, aop):
    S = -norm(getSunPosition(jd))
    return rotateOrderTo(ZXZ, EulerAngles(raan, inc, aop), S)

def getRe(Re, phi, ecc, a, inc, aop, s):
    f = EARTH_FLATTENING
    ae = EARTH_EQUITORIAL_RADIUS

    Rz = _getRz(a, ecc, inc, aop, phi, s)
    fTerm = 2 * f - f * f
    mainTerm = ae * ae * (1 - fTerm)
    cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
    return ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)

def _getRz(a, ecc, inc, aop, phi, s):
    return (a * (1 - ecc * ecc) / (1 + ecc * cos(phi))) * (
            sin(aop) * sin(inc) * cos(phi) + cos(aop) * sin(inc) * sin(phi) - s[2] * (
            s[0] * cos(phi) + s[1] * sin(phi)))
