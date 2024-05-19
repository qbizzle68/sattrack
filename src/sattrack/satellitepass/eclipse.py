"""Ref[1] = SHADOW TIMES OF EARTH SATELLITES, Alessandro de Iaco Veris - Rivista Italiana di Compositi e
Nanotecnologie – Volume 9, n°1 Giugno 2014

Ref[1] defines the algorithm used for computing the shadow times of earth satellites, incorporating the
oblateness of the earth, correction for umbra/penumbra, and atmospheric refraction. These yield very
accurate results; so much so that it is worth it to compute them and compare their positions to satellite
anomalies to determine if the satellite is eclipsed, as opposed to basic geometry assuming round earth
and avoiding specifics like umbra/penumbra and atmospheric refraction."""

from math import atan2, cos, sin, pi, sqrt, acos, asin
from typing import TYPE_CHECKING

from pyevspace.core import norm, Angles, rotateEulerTo, ZXZ, cross, dot

from sattrack.bodies.sun import Sun
from sattrack.core.exceptions import SattrackException
from sattrack.satellitepass.exceptions import NoFunctionRootFound, NoSatelliteEclipseException
from sattrack.util.constants import TWOPI, SUN_RADIUS, EARTH_EQUITORIAL_RADIUS
from sattrack.util.helpers import computeAngleDifference

if TYPE_CHECKING:
    from sattrack.core.juliandate import JulianDate
    from sattrack.orbit.satellite import Orbitable
    from pyevspace import Vector

# Enumerations to distinguishing significant values in finding eclipse positions.
UMBRA = 0x0
PENUMBRA = 0x1
ANNULAR = 0x2
ENTER = 0x10
EXIT = 0x11


class EclipseFinder:
    __slots__ = '_satellite', '_sVector', '_gamma', '_R', '_elements', '_zeta',

    def __init__(self, satellite: 'Orbitable'):
        self._satellite = satellite
        self._R = 6371
        self._elements = self._sVector = self._zeta = None

    def _computeSVector(self, sunPosition: 'Vector') -> 'Vector':
        """Computes the vector s from Ref[1]."""

        sVector = -norm(sunPosition)
        elements = self._elements
        angles = Angles(elements.raan, elements.inc, elements.aop)

        return rotateEulerTo(ZXZ, angles, sVector)

    @staticmethod
    def _computeGamma(sVector: 'Vector') -> float:
        """Computes the vector γ from Ref[1]."""

        return atan2(sVector[1], -sVector[0])

    def _zeroFunction(self, phi: float, shadow: int) -> float:
        """The modified shadow zero function g(φ) for Escobal's method of finding shadow times
        from Ref[1]. The shadow parameter must be UMBRA or PENUMBRA and is used to determine
        the sign of certain terms."""

        R = self._R
        zeta = self._zeta
        cosPhi = cos(phi)
        cosZeta = cos(zeta)
        ecc = self._elements.ecc
        eTerm = 1 + ecc * cosPhi
        aTerm = self._elements.sma * (1 - ecc * ecc)
        sTerm = -self._sVector[0] * cosPhi - self._sVector[1] * sin(phi)

        term1 = R * R * eTerm * eTerm
        term2 = aTerm * aTerm * sTerm * sTerm
        term3 = aTerm * aTerm * cosZeta * cosZeta
        term4 = 2 * aTerm * R * sTerm * eTerm * sin(zeta)

        # subtract term4 for penumbra
        if shadow == PENUMBRA:
            term4 = -term4

        return term1 + term2 - term3 + term4

    def _zeroFunctionDerivative(self, phi: float, shadow: int) -> float:
        """The derivative of the modified shadow zero function g(φ) for Escobal's method of finding
        shadow times in Ref[1]."""

        ecc = self._elements.ecc
        sVector = self._sVector
        R = self._R
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        eTerm = 1 + ecc * cosPhi
        aTerm = self._elements.sma * (1 - ecc * ecc)
        sTerm = -sVector[0] * cosPhi - sVector[1] * sinPhi
        eTermPrime = -ecc * sinPhi
        sTermPrime = sVector[0] * sinPhi - sVector[1] * cosPhi

        term1 = 2 * R * R * eTerm * eTermPrime
        term2 = 2 * aTerm * aTerm * sTerm * sTermPrime
        term4 = 2 * R * aTerm * sin(self._zeta) * (sTerm * eTermPrime + sTermPrime * eTerm)

        if shadow == PENUMBRA:
            term4 *= -1

        return term1 + term2 + term4

    def _computeZeroNewton(self, guess: float, shadow: int, epsilon: float = 1e-5) -> float:
        """Find the zeros of the zero function g(φ) using Newton-Raphson method."""

        phi = guess
        gi = self._zeroFunction(phi, shadow)

        while abs(gi) > epsilon:
            giPrime = self._zeroFunctionDerivative(phi, shadow)
            phi = phi - (gi / giPrime)
            gi = self._zeroFunction(phi, shadow)

        return phi % TWOPI

    def _computeZeroFast(self, shadow: int, direction: int, epsilon: float = 1e-5) -> float:
        """Compute the zero of the function g(φ) using Newton-Raphson method with initial guess being
        (3 * pi / 4 - gamma) if direction is ENTER, and (5 * pi / 4 - gamma) if direction is EXIT. ValueError is
        raised if direction is any other value."""

        if direction is ENTER:
            guess = 3 * pi / 4 - self._gamma
        elif direction is EXIT:
            guess = 5 * pi / 4 - self._gamma
        else:
            raise ValueError('direction parameter value must be ENTER or EXIT')

        return self._computeZeroNewton(guess, shadow, epsilon)

    @staticmethod
    def checkZero(phi: float, sVector: 'Vector') -> bool:
        """Check the root (zero) phi is on the shadow side of the Earth."""
        return (sVector[0] * cos(phi) + sVector[1] * sin(phi)) > 0

    @staticmethod
    def checkAnomalyRange(psi: float, direction: int) -> bool:
        """Check if a given anomaly corresponds to the correct direction. This is used to ensure a root
        anomaly is indeed the entrance or exit anomaly as believed, and should catch any errors where the
        root finding algorithm converges on an unintended root."""

        return (direction == ENTER and pi / 2 <= psi < pi) or (direction == EXIT and pi <= psi <= 3 * pi / 2)

    def _computeAllZeros(self, shadow: int, epsilon: float = 1e-5, iterationLimit: int = 4) -> list[float]:
        """Attempts to compute all possible roots to the zero function g(φ) by continually doubling the
        number of initial evenly spaced guesses. To avoid infinite loops on cases when shadow intersections
        do not exist, only iterationLimit number of iterations are run before a SattrackException is raised.
        Any zeros returned have already been validated by the _checkZero() method."""
        # It is possible that a single root is found, but not its pair if they are very close together but rounding
        # makes them look the same, and the pair is not added to the 'rounded' set. Decreasing the epsilon
        # value should help this.

        iterationCount = 1
        fraction = 0.5
        while True:
            guesses = [i * pi * fraction for i in range(int(2 / fraction))]
            values = [self._computeZeroNewton(guess, shadow, epsilon) for guess in guesses]
            # round to the place of an epsilon value to compare for duplicate floating values ignoring rounding errors
            invertedEpsilon = 1 / epsilon
            rounded = {round(value * invertedEpsilon) / invertedEpsilon for value in values}
            zeros = []
            for phi in values:
                roundedPhi = round(phi * invertedEpsilon) / invertedEpsilon
                if roundedPhi in rounded:
                    zeros.append(phi)
                    rounded.remove(roundedPhi)

            validZeros = [phi for phi in zeros if self.checkZero(phi, self._sVector)]
            if len(validZeros) == 2:
                return validZeros
            elif iterationCount >= iterationLimit:
                raise SattrackException('satellite eclipse limit exceeded')

            fraction /= 2
            iterationCount += 1

    def _computeShadowRoot(self, shadow: int, direction: int, guess: float = None, *, epsilon: float = 1e-5,
                           iterationLimit: int = 4) -> float:
        """Contains the logic for computing a shadow root from the zero function. The parameters
        shadow and direction determine which root is considered valid. If guess is provided, it
        is used as the initial guess to _computeZeroNewton. If that value does not converge correctly,
        _computeZeroFast is used to compute the most efficient starting guess. If this also fails,
        _computeAllZeros is used as a brute force method of finding all zeros, and then discerning
        the correct root needed. The parameter iterationLimit is used in the brute force method, and
        if more than n iterations are needed a NoFunctionRootFound exception is raised."""

        # Initially either use the guess passed in, or use gamma as the initial guess.
        sVector = self._sVector
        gamma = self._gamma = self._computeGamma(sVector)
        if guess is not None:
            phi = self._computeZeroNewton(guess, shadow, epsilon=epsilon)
            if self.checkZero(phi, sVector):
                psi = (gamma + phi) % TWOPI
                if self.checkAnomalyRange(psi, direction):
                    return phi
        # If the root found wasn't correct, use _computeZeroFast to make effective guesses for us.
        phi = self._computeZeroFast(shadow, direction, epsilon)
        if self.checkZero(phi, sVector):
            psi = (gamma + phi) % TWOPI
            if self.checkAnomalyRange(psi, direction):
                return phi

        # Finally, fall back on brute force attack.
        try:
            zeros = self._computeAllZeros(shadow, epsilon, iterationLimit)
        except SattrackException as e:
            if e.args[0].startswith('satellite eclipse limit exceeded'):
                raise NoFunctionRootFound('unable to find a valid zero') from None
            raise e
        else:
            for phi in zeros:
                psi = (gamma + phi) % TWOPI
                if self.checkAnomalyRange(psi, direction):
                    return phi

        raise NoFunctionRootFound('unable to find a valid zero')

    def _computeApertureAngle(self, rs: float, shadow: int) -> float:
        """Compute the aperture angle ξ from Ref[1] in correcting for umbra-penumbra. The shadow
        parameter must be UMBRA or PENUMBRA."""

        if shadow == UMBRA:
            cosZeta = sqrt(rs * rs - (SUN_RADIUS - self._R) ** 2) / rs
        elif shadow == PENUMBRA:
            cosZeta = sqrt(rs * rs - (SUN_RADIUS + self._R) ** 2) / rs
        else:
            raise ValueError('Shadow parameter must be either UMBRA or PENUMBRA')

        return acos(cosZeta)

    @staticmethod
    def _computeRefractionAngle(altitudeAngle: float) -> float:
        """Compute the refraction angle Δξ from Ref[1] to adjust the aperture angle."""

        numerator = 0.009928887226387075 + altitudeAngle * (0.06995 + altitudeAngle * 0.004087098938599872)
        denominator = 1 + altitudeAngle * (28.934368654106574 + altitudeAngle * 277.39713657599236)
        return numerator / denominator

    def _computeCorrectedRefractionAngle(self, rs: float, shadow: int) -> float:
        """Compute the aperture angle, corrected for atmospheric refraction from Ref[1].
        The shadow parameter must be UMBRA or PENUMBRA."""

        semiApertureAngle = self._computeApertureAngle(rs, shadow)
        refractionAngle = self._computeRefractionAngle(semiApertureAngle)
        if shadow == UMBRA:
            correctedAngle = semiApertureAngle + refractionAngle
        elif shadow == PENUMBRA:
            correctedAngle = refractionAngle - semiApertureAngle
        else:
            raise ValueError('shadow parameter must be either UMBRA or PENUMBRA')

        return correctedAngle

    def _computeApproximateAnomalies(self, time: 'JulianDate', shadow: int) -> (float, float):
        """Compute the first iteration of finding the shadow anomalies. These anomalies
        still need refining of the perspective Earth radius, but give us a close approximation
        to determine how to compute the times to each anomaly (forward or backwards). The shadow
        parameter must be UMBRA or PENUMBRA."""

        self._R = 6371
        sunPosition = Sun.computePosition(time)
        self._elements = self._satellite.getElements(time)
        self._sVector = self._computeSVector(sunPosition)
        self._zeta = self._computeCorrectedRefractionAngle(sunPosition.mag(), shadow)

        approxPhi1 = self._computeShadowRoot(shadow, ENTER)
        approxPhi2 = self._computeShadowRoot(shadow, EXIT)

        return approxPhi1, approxPhi2

    def _computePerspectiveRadius(self, phi: float) -> float:
        """Execute an iteration of computing the perspective Earth radius of Escobal's method
        found in Ref[1]."""

        # Compute Z-component of radius
        cosPhi = cos(phi)
        sinPhi = sin(phi)
        sinInc = sin(self._elements.inc)
        sVector = self._sVector
        ecc = self._elements.ecc
        aop = self._elements.aop
        term1 = (self._elements.sma * (1 - ecc * ecc)) / (1 + ecc * cosPhi)
        term2 = sin(aop) * sinInc * cosPhi
        term3 = cos(aop) * sinInc * sinPhi
        term4 = sVector[2] * (sVector[0] * cosPhi + sVector[1] * sinPhi)
        Rz = term1 * (term2 + term3 - term4)

        # Compute latitude term
        ae = EARTH_EQUITORIAL_RADIUS
        fTerm = 0.006694317778266723
        aeTerm = ae * ae * (1 - fTerm)
        rzSquared = Rz * Rz
        latitudeTerm = (aeTerm - rzSquared) / (aeTerm - (rzSquared * fTerm))

        # Radius from latitude
        return (ae * sqrt(1 - fTerm)) / sqrt(1 - fTerm * latitudeTerm)

    def _computeAnomalyLoop(self, time: 'JulianDate', shadow: int, direction: int, *,
                            zeroEpsilon: float = 1e-5, radiusEpsilon: float = 1e-5) -> (float, 'JulianDate'):
        """Run the loop of iterating the algorithm until the radius value stops changing significantly.
        The shadow parameter must be UMBRA or PENUMBRA, and the direction parameter must be ENTER or EXIT."""

        self._R = 6371
        sunPosition = Sun.computePosition(time)
        self._elements = self._satellite.getElements(time)
        self._sVector = self._computeSVector(sunPosition)

        while True:
            self._zeta = self._computeCorrectedRefractionAngle(sunPosition.mag(), shadow)
            phi = self._computeShadowRoot(shadow, direction, epsilon=zeroEpsilon)
            time = self._satellite.timeToNearestAnomaly(phi, time, 'true')
            self._elements = self._satellite.getElements(time)
            sunPosition = Sun.computePosition(time)
            self._sVector = self._computeSVector(sunPosition)
            previousRe = self._R
            self._R = self._computePerspectiveRadius(phi)
            if abs(self._R - previousRe) < radiusEpsilon:
                break

        return phi, time

    def _eclipseIsValid(self, time: 'JulianDate'):
        state = self._satellite.getState(time)
        angularMomentum = cross(*state)
        uz = norm(angularMomentum)
        sunPosition = Sun.computePosition(time)
        s = norm(sunPosition)
        tmp = acos(dot(uz, s))
        if tmp > pi / 2:
            tmp = pi - tmp
        delta = pi / 2 - tmp

        elements = self._elements = self._satellite.getElements(time)
        sVector = self._computeSVector(sunPosition)
        gamma = self._computeGamma(sVector)
        ecc = elements.ecc
        numerator = EARTH_EQUITORIAL_RADIUS * (1 + ecc * cos(pi - gamma))
        denominator = elements.sma * (1 - ecc * ecc)
        deltaPrime = asin(numerator / denominator)

        return delta < deltaPrime

    def computeShadowPositions(self, time: 'JulianDate', shadow: int, *, zeroEpsilon: float = 1e-5,
                               radiusEpsilon: float = 1e-5) -> ((float, 'JulianDate'), (float, 'JulianDate')):
        """Compute the anomaly and time of the entrance and exit positions of the satellite. If the satellite
        is eclipsed at time, the previous occurrence of the entrance time and the next occurrence of the exit
        time is found. Otherwise, the next occurrence of each instant is found. Shadow must be UMBRA or PENUMBRA."""

        if not self._eclipseIsValid(time):
            raise NoSatelliteEclipseException(f'{self._satellite.name} is not eclipsed by Earth\'s shadow')

        phi0 = self._satellite.anomalyAtTime(time, 'true')
        try:
            approxPhi1, approxPhi2 = self._computeApproximateAnomalies(time, shadow)
        except NoFunctionRootFound:
            raise NoSatelliteEclipseException(f'{self._satellite.name} is not eclipsed by Earth\'s shadow') from None

        deltaAnomaly1 = computeAngleDifference(approxPhi1 - phi0)
        deltaAnomaly2 = computeAngleDifference(approxPhi2 - phi0)
        if deltaAnomaly1 < 0:
            if deltaAnomaly2 < 0:
                enterTime = self._satellite.timeToNextAnomaly(approxPhi1, time, 'true')
                exitTime = self._satellite.timeToNextAnomaly(approxPhi2, time, 'true')
            elif deltaAnomaly2 > 0:
                enterTime = self._satellite.timeToPreviousAnomaly(approxPhi1, time, 'true')
                exitTime = self._satellite.timeToNextAnomaly(approxPhi2, time, 'true')
            else:  # deltaAnomaly2 == 0
                enterTime = self._satellite.timeToPreviousAnomaly(approxPhi1, time, 'true')
                exitTime = time
        elif deltaAnomaly1 > 0:
            enterTime = self._satellite.timeToNextAnomaly(approxPhi1, time, 'true')
            exitTime = self._satellite.timeToNextAnomaly(approxPhi2, time, 'true')
        else:
            enterTime = time
            exitTime = self._satellite.timeToNextAnomaly(approxPhi2, time, 'true')

        try:
            enterPhi, enterTime = self._computeAnomalyLoop(enterTime, shadow, ENTER, zeroEpsilon=zeroEpsilon,
                                                           radiusEpsilon=radiusEpsilon)
            exitPhi, exitTime = self._computeAnomalyLoop(exitTime, shadow, EXIT, zeroEpsilon=zeroEpsilon,
                                                         radiusEpsilon=radiusEpsilon)
        except NoFunctionRootFound:
            raise NoSatelliteEclipseException(f'{self._satellite.name} is not eclipsed by Earth\'s shadow') from None

        # Shadow times should not be before the calling time
        if enterTime < exitTime < time:
            gamma = self._computeGamma(self._sVector)
            updatedTime = self._satellite.timeToNextAnomaly(gamma, time, 'true')
            return self.computeShadowPositions(updatedTime, shadow, zeroEpsilon=zeroEpsilon,
                                               radiusEpsilon=radiusEpsilon)

        return (enterPhi, enterTime), (exitPhi, exitTime)

    def computeShadowAnomalies(self, time: 'JulianDate', shadow: int, *, zeroEpsilon: float = 1e-5,
                               radiusEpsilon: float = 1e-5) -> (float, float):
        """Convenience function to only return the anomaly values of the eclipse positions. This is just a
        simple wrapper around computeShadowPositions(). The shadow parameter must be UMBRA or PENUMBRA,
        and a NoSatelliteEclipseException is raised if the satellite is not eclipsed by Earth's shadow."""

        (enterPhi, _), (exitPhi, _) = self.computeShadowPositions(time, shadow, zeroEpsilon=zeroEpsilon,
                                                                  radiusEpsilon=radiusEpsilon)

        return enterPhi, exitPhi

    def computeShadowTimes(self, time: 'JulianDate', shadow: int, *, zeroEpsilon: float = 1e-5,
                           radiusEpsilon: float = 1e-5) -> ('JulianDate', 'JulianDate'):
        """Convenience function to only return the anomaly times of the eclipse positions. This is just a
        simple wrapper around computeShadowPositions(). The shadow parameter must be UMBRA or PENUMBRA,
        and a NoSatelliteEclipseException is raised if the satellite is not eclipsed by Earth's shadow."""

        (_, enterTime), (_, exitTime) = self.computeShadowPositions(time, shadow, zeroEpsilon=zeroEpsilon,
                                                                    radiusEpsilon=radiusEpsilon)

        return enterTime, exitTime


def isEclipsed(satellite: 'Orbitable', time: 'JulianDate', shadow: int = UMBRA, *, zeroEpsilon: float = 1e-5,
               radiusEpsilon: float = 1e-5) -> bool:
    """Returns if the satellite is eclipsed by Earth's shadow at time. If a NoSatelliteEclipseException
    would ordinarily be raised, False is returned instead of propagating the exception. The shadow
    parameter must be UMBRA or PENUMBRA."""

    finder = EclipseFinder(satellite)

    try:
        enterTime, exitTime = finder.computeShadowTimes(time, shadow, zeroEpsilon=zeroEpsilon,
                                                        radiusEpsilon=radiusEpsilon)
    except NoSatelliteEclipseException:
        return False

    return enterTime <= time < exitTime
