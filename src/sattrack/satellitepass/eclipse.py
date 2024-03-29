from enum import Enum
from math import atan2, cos, sin, pi, sqrt, acos, asin

from pyevspace import Vector, norm, vang, ZXZ, Angles, rotateEulerTo, cross, dot

from sattrack.bodies.sun import Sun
from sattrack.satellitepass.exceptions import NoSatelliteEclipseException, NoFunctionRootFound
from sattrack.util.constants import TWOPI, EARTH_EQUITORIAL_RADIUS, SUN_RADIUS

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate
    from sattrack.orbit.satellite import Orbitable


class Shadow(Enum):
    UMBRA = 0
    PENUMBRA = 1
    ANNULAR = 2


class Eclipse(Enum):
    ENTER = 0
    EXIT = 1


def __compute_s_vector(sunPosition: Vector, raan: float, inclination: float, aop: float) -> Vector:
    S = -norm(sunPosition)

    return rotateEulerTo(ZXZ, Angles(raan, inclination, aop), S)


def __compute_gamma(sVector: Vector) -> float:
    return atan2(sVector[1], -sVector[0])


def __escobal_method(R: float, sVector: Vector, phi: float, zeta: float, sma: float, ecc: float,
                     shadow: Shadow) -> float:
    # Escobal's re-defined shadow function from Ref[1]
    cosPhi = cos(phi)
    cosZeta = cos(zeta)
    eTerm = 1 + ecc * cosPhi
    aTerm = sma * (1 - ecc * ecc)
    sTerm = -sVector[0] * cosPhi - sVector[1] * sin(phi)

    term1 = R * R * eTerm * eTerm
    term2 = aTerm * aTerm * sTerm * sTerm
    term3 = aTerm * aTerm * cosZeta * cosZeta
    term4 = 2 * aTerm * R * sTerm * eTerm * sin(zeta)

    # subtract term4 for penumbra
    if shadow is Shadow.PENUMBRA:
        term4 *= -1

    return term1 + term2 - term3 + term4


def __escobal_method_derivative(R: float, sVector: Vector, phi: float, zeta: float, sma: float, ecc: float,
                                shadow: Shadow) -> float:
    # derivative of Escobal's re-defined shadow function from Ref[1]
    cosPhi = cos(phi)
    sinPhi = sin(phi)
    eTerm = 1 + ecc * cosPhi
    aTerm = sma * (1 - ecc * ecc)
    sTerm = -sVector[0] * cosPhi - sVector[1] * sinPhi
    eTermPrime = -ecc * sinPhi
    sTermPrime = sVector[0] * sinPhi - sVector[1] * cosPhi

    term1 = 2 * R * R * eTerm * eTermPrime
    term2 = 2 * aTerm * aTerm * sTerm * sTermPrime
    term4 = 2 * R * aTerm * sin(zeta) * (sTerm * eTermPrime + sTermPrime * eTerm)

    if shadow is Shadow.PENUMBRA:
        term4 *= -1

    return term1 + term2 + term4


# todo: a satellite can have 0, 2, or 4 solutions, need to find solution to this issue
#   idea: test if guessing 8 always hits the unique solutions
#   later if there is only 2 solutions and they don't checkout they're invalid
def __find_zero_newton(R: float, sVector: Vector, guess: float, zeta: float, sma: float, ecc: float, shadow: Shadow,
                       epsilon: float = 1e-5) -> float:
    phi = guess
    gi = __escobal_method(R, sVector, phi, zeta, sma, ecc, shadow)
    while abs(gi) > epsilon:
        phi = phi - (gi / __escobal_method_derivative(R, sVector, phi, zeta, sma, ecc, shadow))
        gi = __escobal_method(R, sVector, phi, zeta, sma, ecc, shadow)
    return phi % TWOPI


def __find_fast_zero(R: float, sVector: Vector, gamma: float, zeta: float, sma: float, ecc: float, shadow: Shadow,
                     enterOrExit: Eclipse, epsilon: float = 1e-5):
    if enterOrExit is Eclipse.ENTER:
        return __find_zero_newton(R, sVector, 3*pi/4 - gamma, zeta, sma, ecc, shadow, epsilon)
    elif enterOrExit is Eclipse.EXIT:
        return __find_zero_newton(R, sVector, 5*pi/4 - gamma, zeta, sma, ecc, shadow, epsilon)
    else:
        raise ValueError('enterOrExit parameter value must be ENTER or EXIT')


def __find_certain_zeros(R: float, sVector: Vector, zeta: float, sma: float, ecc: float, shadow: Shadow,
                         epsilon: float = 1e-5) -> float:
    zeros = []
    frac = 0.5
    while len(zeros) < 4:
        guesses = [i * pi * frac for i in range(int(2 / frac))]
        computed = [__find_zero_newton(R, sVector, guess, zeta, sma, ecc, shadow, epsilon) for guess in guesses]
        # round to the place of an epsilon value to compare for duplicate floating values
        rounded = {round(c * 1e5) / 1e5 for c in computed}
        zeros.clear()
        for phi in computed:
            roundedPhi = round(phi * 1e5) / 1e5
            if roundedPhi in rounded:
                zeros.append(phi)
                rounded.remove(roundedPhi)
        frac /= 2
    return zeros


def __check_zero(phi: float, sVector: Vector):
    return (sVector[0] * cos(phi) + sVector[1] * sin(phi)) > 0


def __check_range(psi: float, enterOrExit: Eclipse) -> bool:
    cmp1 = (enterOrExit is Eclipse.ENTER and pi / 2 <= psi <= pi)
    cmp2 = (enterOrExit is Eclipse.EXIT and pi <= psi <= 3 * pi / 2)
    return cmp1 or cmp2


def _get_zero(R: float, sVector: float, zeta: float, sma: float, ecc: float, shadow: Shadow, enterOrExit: Eclipse, *,
              guess: float = None, epsilon: float = 1e-5) -> float:
    # if guess is set use it first, then try __find_fast_zero, then use __find_certain_zeros
    gamma = __compute_gamma(sVector)
    if guess is not None:
        phi = __find_zero_newton(R, sVector, guess, zeta, sma, ecc, shadow)
        if __check_zero(phi, sVector):
            psi = (gamma + phi) % TWOPI
            if __check_range(psi, enterOrExit):
                return phi
    phi = __find_fast_zero(R, sVector, gamma, zeta, sma, ecc, shadow, enterOrExit, epsilon)
    if __check_zero(phi, sVector):
        psi = (gamma + phi) % TWOPI
        if __check_range(psi, enterOrExit):
            return phi
    zeros = __find_certain_zeros(R, sVector, zeta, sma, ecc, shadow, epsilon)
    for phi in zeros:
        if __check_zero(phi, sVector):
            psi = (gamma + phi) % TWOPI
            if __check_range(psi, enterOrExit):
                return phi

    raise NoFunctionRootFound('unable to find a valid zero')


def __get_radius_z_comp(sVector: Vector, sma: float, ecc: float, inc: float, aop: float, phi: float) -> float:
    cosPhi = cos(phi)
    sinPhi = sin(phi)
    sinInc = sin(inc)
    term1 = (sma * (1 - ecc * ecc)) / (1 + ecc * cosPhi)
    term2 = sin(aop) * sinInc * cosPhi
    term3 = cos(aop) * sinInc * sinPhi
    term4 = sVector[2] * (sVector[0] * cosPhi + sVector[1] * sinPhi)
    return term1 * (term2 + term3 - term4)


def __get_latitude_term(Rz: float) -> float:
    ae = EARTH_EQUITORIAL_RADIUS
    fTerm = -0.006694317778266723
    aeTerm = ae * ae * (1 + fTerm)
    return (aeTerm - Rz * Rz) / (aeTerm + (Rz * Rz * fTerm))


def __get_radius_from_latitude(latitudeTerm: float) -> float:
    ae = EARTH_EQUITORIAL_RADIUS
    fTerm = 0.006694317778266723
    return (ae * sqrt(1 - fTerm)) / sqrt(1 - fTerm * latitudeTerm)


def __get_aperture_angle(rs: float, Re: float, shadow: Shadow) -> float:
    if shadow is Shadow.UMBRA:
        cosZeta = sqrt(rs * rs - (SUN_RADIUS - Re) ** 2) / rs
    elif shadow is Shadow.PENUMBRA:
        cosZeta = sqrt(rs * rs - (SUN_RADIUS + Re) ** 2) / rs
    else:
        raise ValueError('shadow parameter must be either UMBRA or PENUMBRA')
    return acos(cosZeta)


def __get_refraction_angle(altitudeAngle: float) -> float:
    numerator = 0.009928887226387075 + altitudeAngle * (0.06995 + altitudeAngle * 0.004087098938599872)
    denominator = 1 + altitudeAngle * (28.934368654106574 + altitudeAngle * 277.39713657599236)
    return numerator / denominator


def __get_corrected_refraction_angle(rs: float, Re: float, shadow: Shadow) -> float:
    semiApertureAngle = __get_aperture_angle(rs, Re, shadow)
    refractionAngle = __get_refraction_angle(semiApertureAngle)
    if shadow is Shadow.UMBRA:
        correctedAngle = semiApertureAngle + refractionAngle
    elif shadow is Shadow.PENUMBRA:
        correctedAngle = refractionAngle - semiApertureAngle
    else:
        raise ValueError('shadow parameter must be either UMBRA or PENUMBRA')
    return correctedAngle


def _get_shadow_positions(jd: 'JulianDate', sat: 'Orbitable', shadow: Shadow, zeroEpsilon: float = 1e-5,
                          radiusEpsilon: float = 1e-5) -> ((float, 'JulianDate'), (float, 'JulianDate')):
    Re = 6371  # start with average radius of the earth
    time = jd
    sunPosition = Sun.computePosition(time)
    elements = sat.getElements(time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
    apertureAngle = __get_corrected_refraction_angle(sunPosition.mag(), Re, shadow)

    # phi0 = sat.anomalyAt(jd, Orbitable.TRUE)
    phi0 = sat.anomalyAtTime(jd, 'true')
    approxPhi1 = _get_zero(Re, sVector, apertureAngle, elements.sma,
                           elements.ecc, shadow, Eclipse.ENTER, epsilon=zeroEpsilon)
    approxPhi2 = _get_zero(Re, sVector, apertureAngle, elements.sma, elements.ecc,
                           shadow, Eclipse.EXIT, epsilon=zeroEpsilon)

    # maximum difference between the rages of Re seems to be about 0.015, so take a very conservative value of 0.1 rad
    errorBuffer = 0.1
    # we're close enough to the exit anomaly that we don't know if the approximated anomaly is before or after it
    if (abs(phi0 - approxPhi2) < errorBuffer) or (abs(phi0 + TWOPI - approxPhi2) < errorBuffer) \
            or (abs(approxPhi2 + TWOPI - phi0) < errorBuffer):
        # dt = (jd - sat.timeToAnomaly(approxPhi1, jd, Orbitable.PREVIOUS, Orbitable.TRUE))
        dt = (jd - sat.timeToPreviousAnomaly(approxPhi1, jd, 'true'))
        referenceTime = jd.future(-dt)
    else:
        # phi2Time = sat.timeToAnomaly(approxPhi2, jd, Orbitable.NEXT, Orbitable.TRUE)
        # phi1Time = sat.timeToAnomaly(approxPhi1, phi2Time, Orbitable.PREVIOUS, Orbitable.TRUE)
        phi2Time = sat.timeToNextAnomaly(approxPhi2, jd, 'true')
        phi1Time = sat.timeToPreviousAnomaly(approxPhi1, phi2Time, 'true')
        dt = (phi2Time - phi1Time) / 2
        referenceTime = phi1Time.future(dt)

    enterPhi, enterTime = __compute_anomaly_loop(jd, referenceTime, sat, shadow,
                                                 Eclipse.ENTER, zeroEpsilon, radiusEpsilon)
    exitPhi, exitTime = __compute_anomaly_loop(jd, referenceTime, sat, shadow, Eclipse.EXIT, zeroEpsilon, radiusEpsilon)

    # if our approxPhi2 ended up being
    if enterTime < exitTime < jd:
        gamma = __compute_gamma(sVector)
        updatedJd = sat.timeToNextAnomaly(gamma, jd, 'true')
        return _get_shadow_positions(updatedJd, sat, shadow, zeroEpsilon)

    return (enterPhi, enterTime), (exitPhi, exitTime)


def __compute_anomaly_loop(startTime: 'JulianDate', referenceTime: 'JulianDate', sat: 'Orbitable', shadow: Shadow,
                           enterOrExit: Eclipse, zeroEpsilon: float = 1e-5,
                           radiusEpsilon: float = 1e-5) -> (float, 'JulianDate'):
    Re = 6371
    time = startTime
    sunPosition = Sun.computePosition(time)
    elements = sat.getElements(time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)

    # imitating a do-while structure here
    while True:
        apertureAngle = __get_corrected_refraction_angle(sunPosition.mag(), Re, shadow)
        phi = _get_zero(Re, sVector, apertureAngle, elements.sma, elements.ecc, shadow, enterOrExit,
                        epsilon=zeroEpsilon)
        if enterOrExit is Eclipse.ENTER:
            # time = sat.timeToAnomaly(phi, referenceTime, Orbitable.PREVIOUS, Orbitable.TRUE)
            time = sat.timeToPreviousAnomaly(phi, referenceTime, 'true')
        elif enterOrExit is Eclipse.EXIT:
            # time = sat.timeToAnomaly(phi, referenceTime, Orbitable.NEXT, Orbitable.TRUE)
            time = sat.timeToNextAnomaly(phi, referenceTime, 'true')
        elements = sat.getElements(time)
        sunPosition = Sun.computePosition(time)
        sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
        previousRe = Re
        Re = _get_perspective_radius(sVector, elements.sma, elements.ecc, elements.inc, elements.aop, phi)
        if abs(Re - previousRe) > radiusEpsilon:
            break

    return phi, time


def _get_perspective_radius(sVector: Vector, sma: float, ecc: float, inc: float, aop: float, phi: float) -> float:
    Rz = __get_radius_z_comp(sVector, sma, ecc, inc, aop, phi)
    latitudeTerm = __get_latitude_term(Rz)
    return __get_radius_from_latitude(latitudeTerm)


def getShadowPositions(satellite: 'Orbitable', time: 'JulianDate',
                       shadowType: Shadow) -> ((float, 'JulianDate'), (float, 'JulianDate')):
    if checkEclipse(satellite, time) is False:
        raise NoSatelliteEclipseException(f"{satellite.name} is not eclipsed by Earth's shadow")
    return _get_shadow_positions(time, satellite, shadowType)


def getShadowAnomalies(satellite: 'Orbitable', time: 'JulianDate', shadowType: Shadow) -> (float, float):
    if checkEclipse(satellite, time) is False:
        raise NoSatelliteEclipseException(f"{satellite.name} is not eclipsed by Earth's shadow")
    (enterPhi, enterTime), (exitPhi, exitTime) = _get_shadow_positions(time, satellite, shadowType)
    return enterPhi, exitPhi


def getShadowTimes(satellite: 'Orbitable', time: 'JulianDate', shadowType: Shadow) -> ('JulianDate', 'JulianDate'):
    if checkEclipse(satellite, time) is False:
        raise NoSatelliteEclipseException(f"{satellite.name} is not eclipsed by Earth's shadow")
    (enterPhi, enterTime), (exitPhi, exitTime) = _get_shadow_positions(time, satellite, shadowType)
    return enterTime, exitTime


def isEclipsed(satellite: 'Orbitable', time: 'JulianDate', shadowType: Shadow = Shadow.PENUMBRA) -> bool:

    satPosition = satellite.getState(time)[0]
    sunPosition = Sun.computePosition(time)
    # convert vectors to relative to the satellite
    earthPosition = -satPosition
    relativeSunPosition = -satPosition + sunPosition

    # get the perspective radius of the earth towards the sun
    # elements = Elements.fromTle(satellite.getTle(), time)
    elements = satellite.getElements(time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
    # phi = elements.trueAnomalyAt(time)
    # phi = satellite.anomalyAt(time, Orbitable.TRUE)
    phi = satellite.anomalyAtTime(time, 'true')
    earthRadius = _get_perspective_radius(sVector, elements.sma, elements.ecc, elements.inc, elements.aop, phi)

    # semi-diameters of sun and earth relative to the satellite
    thetaE = asin(earthRadius / earthPosition.mag())
    thetaS = asin(SUN_RADIUS / relativeSunPosition.mag())
    theta = vang(earthPosition, relativeSunPosition)

    # logical eclipse values from Ref[2]
    if shadowType is Shadow.UMBRA:
        return (thetaE > thetaS) and (theta < (thetaE - thetaS))
    elif shadowType is Shadow.PENUMBRA:
        return abs(thetaE - thetaS) < theta < (thetaE + thetaS) or (thetaE > thetaS) and (theta < (thetaE - thetaS))
    elif shadowType is Shadow.ANNULAR:
        return (thetaS > thetaE) and (theta < (thetaS - thetaE))


def checkEclipse(sat, time):
    # compute delta
    state = sat.getState(time)
    h = cross(*state)
    uz = norm(h)
    sunPosition = Sun.computePosition(time)
    s = norm(sunPosition)
    tmp = acos(dot(uz, s))
    if tmp > pi / 2:
        tmp = pi - tmp
    delta = pi / 2 - tmp

    # compute delta'
    elements = sat.getElements(time)
    sVector = __compute_s_vector(sunPosition, elements.raan, elements.inc, elements.aop)
    gamma = __compute_gamma(sVector)
    rhs = EARTH_EQUITORIAL_RADIUS * (1 + elements.ecc * cos(pi - gamma))
    lhs = elements.sma * (1 - elements.ecc * elements.ecc)
    deltaP = asin(rhs / lhs)

    return delta < deltaP
