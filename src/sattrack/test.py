from math import acos

from numpy import arange
from pyevspace import norm, dot

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import rotateOrderTo
from sattrack.satpass import nextPassMax
from sattrack.structures.coordinates import *
from sattrack.structures.elements import OrbitalElements
from sattrack.structures.satellite import Satellite
from sattrack.structures.tle import *
from sattrack.sun import getSunPosition, getSunPosition2
from sattrack.topos import getAltitude
from sattrack.util.anomalies import trueToMean
from sattrack.util.constants import SUN_RADIUS
from sattrack.util.conversions import smaToMeanMotion

from sattrack.spacetime.juliandate import now
from sattrack.structures.tle import getTle

tle = getTle('zarya')
# tle = TwoLineElement("""ISS (ZARYA)
# 1 25544U 98067A   22266.84431519  .00008111  00000+0  14870-3 0  9996
# 2 25544  51.6423 207.8056 0002412 286.8120 181.5821 15.50238875360488""")
iss = Satellite(tle)
#tle = getTLE('ixpe')
#ixpe = Satellite(tle)
# jd = JulianDate(7, 20, 2022, 1, 0, 0, -5)
jd = now()
geo = GeoPosition(38.0608, -97.9298)
elements = OrbitalElements.fromTle(tle, jd)

def f(a, e, r, psi, gamma, delta):
    """psi, gamma, delta in radians"""
    return 1 + e * cos(psi-gamma) - (a/r)*(1-e*e)*sqrt(1-(cos(psi)**2)*(cos(delta)**2))


def fprime(a, e, r, psi, gamma, delta):
    """psi, gamma, delta in radians"""
    term1 = -e*sin(psi- gamma)
    term2 = (a/r)*(1-e*e)
    term3 = sin(2*psi) * (cos(delta)**2)
    term4 = 2 * sqrt(1-(cos(psi)**2)*(cos(delta)**2))
    return term1 - term2 * term4 / term3


def g(a, ecc, Re, phi, s):
    """phi in radians"""
    c = (a * (1 - ecc * ecc)) ** 2
    term1 = Re * Re * ((1 + ecc * cos(phi)) ** 2)
    term2 = (-s[0]*cos(phi) - s[1]*sin(phi)) ** 2
    return term1 + c * term2 - c


def gprime(a, ecc, Re, phi, s):
    """phi in radians"""
    term1 = 2 * (1 + ecc * cos(phi)) * (-ecc * sin(phi))
    term2 = 2 * (-s[0]*cos(phi)-s[1]*sin(phi)) * (s[0]*sin(phi)-s[1]*cos(phi))
    return Re * Re * term1 + ((a * (1 - ecc * ecc)) ** 2) * term2


def gi(a, e, r, phiU, phi, s, rs, enterExit):
    """phiU, phi in radians"""
    # enterExit is 1 for enter and 2 for exit
    cosZeta = sqrt((rs-(SUN_RADIUS-r))*(rs+(SUN_RADIUS-r))) / rs
    #cosZeta = 30.932 / 60
    sinZeta = (SUN_RADIUS - r) / rs
    #sinZeta = 14.729 / 60
    term1 = r*r * (1 + e*cos(phiU))**2
    term2 = a * (1 - e*e)
    term3 = -s[0] * cos(phi) - s[1] * sin(phi)
    if enterExit == 1:
        rhs = -2 * term2 * r * term3 * (1 + e*cos(phi)) * sinZeta
    elif enterExit == 2:
        rhs = 2 * term2 * r * term3 * (1 + e*cos(phi)) * sinZeta
    else:
        raise ValueError()
    return term1 + term2 * term2 * term3 * term3 - term2 * term2 * cosZeta * cosZeta + rhs


def giprime(a, e, r, phi, s, rs, enterExit):
    """phi in radians"""
    sinZeta = (SUN_RADIUS - r) / rs
    term1 = a * (1 - e*e)
    term2 = -s[0] * cos(phi) - s[1] * sin(phi)
    term3 = s[0] * sin(phi) - s[1] * cos(phi)
    if enterExit == 1:
        rhs = -2 * term1 * r * sinZeta
    elif enterExit == 2:
        rhs = 2 * term1 * r * sinZeta
    else:
        raise ValueError()
    return 2*term1*term1*term2*term3 + rhs*(term2*(e*sin(phi)) + (1+e*cos(phi))*term3)


def getPhii(a, e, r, phiU, phi, s, rs, enterExit):
    """phiU, phi in radians"""
    nextPhi = phi - (gi(a, e, r, phiU, phi, s, rs, enterExit) / giprime(a, e, r, phi, s, rs, enterExit))
    while abs(phi-nextPhi) > 1e-7:
        phi = nextPhi
        nextPhi = phi - (gi(a, e, r, phiU, phi, s, rs, enterExit) / giprime(a, r, r, phi, s, rs, enterExit))
    return phi


def getPhi(a, ecc, Re, guess, s):
    """guess in radians"""
    phi = guess
    nextPhi = phi - (g(a, ecc, Re, phi, s) / gprime(a, ecc, Re, phi, s))
    #print('phi', phi, 'next', nextPhi)
    while abs(phi-nextPhi) > 1e-8:
        phi = nextPhi
        nextPhi = phi - (g(a, ecc, Re, phi, s) / gprime(a, ecc, Re, phi, s))
        #print('phi', phi, 'next', nextPhi)
    return phi

def getPhi2(a, ecc, Re, guess, s):
    gi = g(a, ecc, Re, guess, s)
    phi = guess - gi / gprime(a, ecc, Re, guess, s)
    while abs(gi) > 1e-5:
        gi = g(a, ecc, Re, phi, s)
        phi = phi - gi / gprime(a, ecc, Re, phi, s)
    return phi


def getRzFromAnomaly(a, e, i, aop, phi, s):
    """i, aop in radians"""
    '''i = radians(i)
    aop = radians(aop)'''
    term1 = (a * (1-e*e)) / (1 + e*cos(phi))
    term2 = sin(aop) * sin(i) * cos(phi)
    term3 = cos(aop) * sin(i) * sin(phi)
    term4 = s[2] * (s[0]*cos(phi) + s[1]*sin(phi))
    return term1 * (term2 + term3 - term4)


def getRzFromCos(cosTerm):
    fTerm = 2 * EARTH_FLATTENING - EARTH_FLATTENING*EARTH_FLATTENING
    term1 = EARTH_EQUITORIAL_RADIUS * sqrt(1 - fTerm)
    term2 = sqrt(1 - fTerm * cosTerm)
    return term1 / term2


def getCosTerm(Rz):
    fTerm = 2 * EARTH_FLATTENING - EARTH_FLATTENING*EARTH_FLATTENING
    term1 = EARTH_EQUITORIAL_RADIUS*EARTH_EQUITORIAL_RADIUS * (1 - fTerm)
    return (term1 - Rz*Rz) / (term1 - Rz*Rz*fTerm)


def getZeros(a, e, s, R):
    frac = 0.5
    twoPi = 2 * pi
    ls = [getPhi(a, e, R, pi*i*frac, s) for i in range(int(2/frac))]
    for i in range(len(ls)):
        if ls[i] >= twoPi:
            ls[i] %= twoPi
    frac /= 2
    while len(ls) < 4:
        tmp = [getPhi(a, e, R, pi * i * frac, s) for i in range(int(2/frac))]
        for i in range(len(tmp)):
            if tmp[i] >= twoPi:
                tmp[i] %= twoPi
        for i in tmp:
            for j in ls:
                if abs(i-j) > 1e-5:
                    ls.append(i)
    return tuple(ls)


def getPsi(a, e, r, guess, gamma, delta):
    psi = guess
    nextPsi = psi - (f(a, e, r, psi, gamma, delta) / fprime(a, e, r, psi, gamma, delta))
    while abs(psi-nextPsi) > 1e-6:
        psi = nextPsi
        nextPsi = psi - (f(a, e, r, psi, gamma, delta) / fprime(a, e, r, psi, gamma, delta))
    return psi


def getSunPos(elements, time: JulianDate):
    sunPos = getSunPosition(time)
    S = -norm(sunPos)
    return rotateOrderTo(ZXZ, EulerAngles(elements.getRaan(), elements.getInc(), elements.getAop()), S)


def getDelta(s: EVector):
    return atan((s[2]*s[2]) / sqrt(s[0]*s[0] + s[1]*s[1]))


def getGamma(s: EVector):
    if s[0] < 0:
        if s[1] == 0:
            return 0
        elif s[1] > 0:
            return atan(-s[1]/s[0])
        else:
            return atan(s[1]/s[0])
    elif s[0] == 0:
            return pi/2
    elif s[0] > 0:
        if s[1] > 0:
            return pi - atan(s[1]/s[0])
        elif s[1] == 0:
            return pi
        else:
            return pi - atan(-s[1]/s[0])


def getShadowAnomalies2(sat: Satellite, time: JulianDate, sunPos: EVector) -> tuple[float]:
    elements = OrbitalElements.fromTle(sat.getTle(), time)
    s = rotateOrderTo(ZXZ, EulerAngles(elements.getRaan(), elements.getInc(), elements.getAop()), -norm(sunPos))
    zeros = getZeros(elements.getSma(), elements.getEcc(), s, 6371)
    checks =[s[0]*cos(z) + s[1]*sin(z) for z in zeros]
    trueAnomalies = [zeros[i] for i in range(4) if checks[i] > 0]
    gamma = getGamma(s)
    psi1, psi2 = trueAnomalies[0] + gamma, trueAnomalies[1] + gamma
    if psi1 > pi:
        return trueAnomalies[1], trueAnomalies[0]
    else:
        return trueAnomalies[0], trueAnomalies[1]


def adjustForNonSphericalEarth(sat: Satellite, time: JulianDate, anomaly: float):
    elements = OrbitalElements.fromTle(sat.getTle(), time)
    a, e, i, aop = elements.getSma(), elements.getEcc(), elements.getInc(), elements.getAop()
    s = getSunPos(elements, time)
    Rz = getRzFromAnomaly(a, e, i, aop, anomaly, s)
    cosTerm = getCosTerm(Rz)
    R = 6371
    nextR = getRzFromCos(cosTerm)
    phi = anomaly
    nextPhi = getPhi(a, e, nextR, phi, s)
    while abs(R-nextR) > 1e-4:
        phi = nextPhi
        Rz = getRzFromAnomaly(a, e, i, aop, phi, s)
        cosTerm = getCosTerm(Rz)
        R = nextR
        nextR = getRzFromCos(cosTerm)
        nextPhi = getPhi(a, e, R, phi, s)
    return nextPhi


def getShadowTimes2(sat: Satellite, time: JulianDate, sunPos: EVector) -> tuple[JulianDate]:
    elements = OrbitalElements.fromTle(sat.getTle(), time)

    ta1, ta2 = getShadowAnomalies2(sat, time, sunPos)
    ta1, ta2 = adjustForNonSphericalEarth(sat, jd, ta1), adjustForNonSphericalEarth(sat, jd, ta2)
    m0 = radians(elements.getMeanAnomaly())
    m1 = radians(trueToMean(degrees(ta1), elements.getEcc()))
    m2 = radians(trueToMean(degrees(ta2), elements.getEcc()))
    twoPi = 2 * pi
    if m1 <= m2:
        if m0 <= m1:
            dm1 = m1 - m0
            dm2 = m2 - m0
        elif m0 <= m2:
            dm1 = m1 + twoPi - m0
            dm2 = m2 - m0
        else:
            dm1 = m1 + twoPi - m0
            dm2 = m2 + twoPi - m0
    else:
        if m0 <= m2:
            dm1 = m1 - m0
            dm2 = m2 - m0
        elif m0 <= m1:
            dm1 = m1 - m0
            dm2 = m2 + twoPi - m0
        else:
            dm1 = m1 + twoPi - m0
            dm2 = m2 + twoPi - m0

    n = smaToMeanMotion(elements.getSma())
    dt1 = dm1 / n / 86400.0
    dt2 = dm2 / n / 86400.0
    return time.future(dt1), time.future(dt2)


def getShadowAnomalies(sat: Satellite, time: JulianDate, sunPos: EVector) -> tuple[float]:
    elements = OrbitalElements.fromTle(sat.getTle(), time)
    S = -norm(sunPos)
    s = rotateOrderTo(ZXZ, EulerAngles(elements.getRaan(), elements.getInc(), elements.getAop()), S)
    gamma = getGamma(s)
    delta = atan((s[2]*s[2]) / sqrt(s[0]*s[0] + s[1]*s[1]))
    twoPi = 2 * pi
    psi1 = getPsi(
        elements.getSma(),
        elements.getEcc(),
        6371.0,
        (3 * pi / 4),
        gamma,
        delta
    )
    psi2 = getPsi(
        elements.getSma(),
        elements.getEcc(),
        6371.0,
        (5 * pi / 4),
        gamma,
        delta
    )
    print('psi1', psi1, 'psi2', psi2)
    print('gamma', gamma)
    ta1 = (psi1 - gamma) % twoPi
    ta2 = (psi2 - gamma) % twoPi
    print('ta1', ta1, 'ta2', ta2)
    if ta2 < ta1:
    #if psi2 - gamma < psi1 - gamma:
        tmp = ta1
        ta1 = ta2
        ta2 = tmp
    return ta1 % twoPi, ta2 % twoPi
    #return ((psi1 - gamma) % (2 * pi), (psi2 - gamma) % (2 * pi))


def getShadowTimes(sat: Satellite, time: JulianDate, sunPos: EVector) -> tuple[JulianDate]:
    """jd2 is the next time the satellite leaves the earth's shadow, jd1 is the preceding time the
    satellite enters it.
    """
    tAnom1, tAnom2 = getShadowAnomalies(sat, time, sunPos)
    elements = OrbitalElements.fromTle(sat.getTle(), time)
    m1 = radians(trueToMean(degrees(tAnom1), elements.getEcc()))
    m2 = radians(trueToMean(degrees(tAnom2), elements.getEcc()))
    m0 = radians(elements.getMeanAnomaly())
    twoPi = 2 * pi
    if m1 <= m2:
        if m0 <= m1:
            dm1 = m1 - m0
            dm2 = m2 - m0
        elif m0 <= m2:
            dm1 = m1 + twoPi - m0
            dm2 = m2 - m0
        else:
            dm1 = m1 + twoPi - m0
            dm2 = m2 + twoPi - m0
    else:
        if m0 <= m2:
            dm1 = m1 - m0
            dm2 = m2 - m0
        elif m0 <= m1:
            dm1 = m1 - m0
            dm2 = m2 + twoPi - m0
        else:
            dm1 = m1 + twoPi - m0
            dm2 = m2 + twoPi - m0

    print('dm1', dm1, 'dm2', dm2)
    n = smaToMeanMotion(elements.getSma())
    dt1 = dm1 / n / 86400.0
    dt2 = dm2 / n / 86400.0
    print('dt1', dt1, 'dt2', dt2)
    return time.future(dt1), time.future(dt2)

    """jd1 is the next time to enter the shadow, jd2 is the next time to leave the shadow"""
    '''tAnom1, tAnom2 = getShadowAnomalies(sat, time, sunPos)
    elements = OrbitalElements.fromTle(sat.tle(), time)
    dt1 = timeToMeanAnomaly(elements, trueToMean(degrees(tAnom1), elements.getEcc()))
    dt2 = timeToMeanAnomaly(elements, trueToMean(degrees(tAnom2), elements.getEcc()))
    return time.future(dt1), time.future(dt2)'''


def findAnomalies(elements: OrbitalElements, s) -> tuple[float]:
    zeros = getZeros(elements.getSma(), elements.getEcc(), s, 6371)
    checks = [s[0] * cos(z) + s[1] * sin(z) for z in zeros]
    trueAnomalies = [zeros[i] for i in range(4) if checks[i] > 0]
    gamma = getGamma(s)
    if trueAnomalies[0] + gamma > pi:
        phi1, phi2 = trueAnomalies[1], trueAnomalies[0]
    else:
        phi1, phi2 = trueAnomalies[0], trueAnomalies[1]

    # find correct R's
    R1 = getReFlattening(elements, s, phi1)
    R2 = getReFlattening(elements, s, phi2)

    zeros1 = getZeros(elements.getSma(), elements.getEcc(), s, R1)
    print(zeros1)
    zeros2 = getZeros(elements.getSma(), elements.getEcc(), s, R2)
    dif1 = [abs(phi1 - z) for z in zeros1]
    print(dif1)
    dif2 = [abs(phi2 - z) for z in zeros2]
    phi1Adjust = zeros1[dif1.index(min(dif1))]
    phi2Adjust = zeros2[dif2.index(min(dif2))]
    return phi1Adjust, phi2Adjust

def getReFlattening(elements: OrbitalElements, s,  phi) -> float:
    a = elements.getSma()
    e = elements.getEcc()
    i = elements.getInc()
    aop = elements.getAop()
    f = EARTH_FLATTENING
    ae = EARTH_EQUITORIAL_RADIUS
    R = 6371
    Rz = (a*(1-e*e) / (1 + e*cos(phi))) * (sin(aop)*sin(i)*cos(phi) + cos(aop)*sin(i)*sin(phi)
                                           -s[2]*(s[0]*cos(phi) + s[1]*sin(phi)))
    fTerm = 2 * f - f * f
    mainTerm = ae * ae * (1 - fTerm)
    cosTerm = (mainTerm - Rz*Rz) / (mainTerm - Rz*Rz*fTerm)
    nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
    while abs(nextR - R) > 1e-5:
        R = nextR
        phi = getPhi(a, e, R, phi, s)
        print('phi:', phi)
        Rz = (a * (1 - e * e) / (1 + e * cos(phi))) * (sin(aop) * sin(i) * cos(phi) + cos(aop) * sin(i) * sin(phi)
                                                       - s[2] * (s[0] * cos(phi) + s[1] * sin(phi)))
        fTerm = 2 * f - f * f
        mainTerm = ae * ae * (1 - fTerm)
        cosTerm = (mainTerm - Rz * Rz) / (mainTerm - Rz * Rz * fTerm)
        nextR = ae * sqrt(1 - fTerm) / sqrt(1 - fTerm * cosTerm)
    print('phi:', getPhi(a, e, R, phi, s))
    return nextR


def findUmbra(sat, time, Re):
    sunPos = getSunPosition2(time)
    S = -norm(sunPos)
    s = rotateOrderTo(ZXZ, EulerAngles(elements.getRaan(), elements.getInc(), elements.getAop()), S)

    rs = sunPos.mag()
    cosZ = sqrt(rs*rs - ((SUN_RADIUS - Re) ** 2)) / rs
    sinZ = (SUN_RADIUS - Re) / rs

    return s, cosZ, sinZ


def gi2(elements, s, Re, cosZ, sinZ, phi):
    a = elements.getSma()
    e = elements.getEcc()


    term1 = Re * Re * (1 + e * cos(phi)) ** 2
    term2 = a * (1 - e * e)
    term3 = -s[0] * cos(phi) - s[1] * sin(phi)
    rhs = -2 * term2 * Re * term3 * (1 + e * cos(phi)) * sinZ
    
    return term1 + term2 * term2 * term3 * term3 - term2 * term2 * cosZ * cosZ + rhs


def getGeometry(Dp, Ds, deltaSP):
    Xu = (Dp * deltaSP) / (Ds - Dp)
    Xp = (Dp * deltaSP) / (Ds + Dp)
    alphaU = asin(Dp / (2 * Xu))
    alphaP = asin(Dp / (2 * Xp))
    return Xu, alphaU, Xp, alphaP


def getTerminatorPoints(Xu, alphaU, Xp, alphaP, rs):
    ksi = (Xu - rs) * tan(alphaU)
    kappa = (Xp + rs) * tan(alphaP)
    return ksi, kappa


if __name__ == '__main__':
    nextjd = nextPassMax(iss, geo, jd)
    for i in arange(-10/86400, 10/86400, 1/86400):
        print('time:', nextjd.future(i), 'alt:', getAltitude(iss, nextjd.future(i), geo))
