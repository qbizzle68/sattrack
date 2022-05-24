from pyevspace import EVector

from anomalies import atan2
from spacetime import JulianDate
from tle import TwoLineElement
from math import radians, cos, sqrt, sin
from constants import CK2, CK4, E6A, Q0MS2T, S, TOTHRD, XJ3, XKE, XKMPER, XMNPDA, \
    AE, TWOPI


class _Mean_Elements:

    def __init__(self, tle: TwoLineElement):
        self.xm0 = radians(tle.meanAnomaly())
        self.xnode0 = radians(tle.raan())
        self.omega0 = radians(tle.argumentOfPeriapsis())
        self.e0 = tle.eccentricity()
        self.xincl = radians(tle.inclination())
        self.xn0 = tle.meanMotion() * TWOPI / XMNPDA
        self.xndt20 = tle.meanMotionDot() * TWOPI / (XMNPDA * XMNPDA)
        self.xndd60 = tle.meanMotionDDot() * TWOPI / (XMNPDA * XMNPDA * XMNPDA)


class _SGP4_Independent:

    def __init__(self, tle: TwoLineElement):
        self._me = _Mean_Elements(tle)
        self.initialize(tle)

    def initialize(self, tle: TwoLineElement) -> None:

        # Recover original mean motion (xn0dp) and semi-major axis (a0dp) from input elements
        a1 = pow(XKE / self._me.xn0, TOTHRD)
        self._cosi0 = cos(self._me.xincl)
        theta2 = self._cosi0 * self._cosi0
        self._x3thm1 = 3 * theta2 - 1
        e0sq = self._me.e0 * self._me.e0
        beta02 = 1 - e0sq
        beta0 = sqrt(beta02)
        del1 = 1.5 * CK2 * self._x3thm1 / (a1 * a1 * beta0 * beta02)
        a0 = a1 * (1 - del1 * (0.5 * TOTHRD + del1 * (1 + 134.0 / 81.0 * del1)))
        del0 = 1.5 * CK2 * self._x3thm1 / (a0 * a0 * beta0 * beta02)
        self._xn0dp = self._me.xn0 / (1 + del0)
        self._a0dp = a0 / (1 - del0)

        # Initialization
        if (self._a0dp * (1 - self._me.e0) / AE) < (220 / XKMPER + AE):
            self._isImp = True
        self._s4 = S
        self._q0ms24 = Q0MS2T
        perige = (self._a0dp * (1 - self._me.e0) - AE) * XKMPER
        if perige < 156:
            if perige <= 98:
                self._s4 = 20
            else:
                self._s4 = perige - 78
            self._q0ms24 = pow((120 - self._s4) * AE / XKMPER, 4)
            self._s4 /= (XKMPER + AE)
        pinvsq = 1.0 / (self._a0dp * self._a0dp * beta02 * beta02)
        self._tsi = 1 / (self._a0dp - self._s4)
        self._eta = self._a0dp * self._me.e0 * self._tsi
        etasq = self._eta * self._eta
        eeta = self._me.e0 * self._eta
        psisq = abs(1 - etasq)
        coef = self._q0ms24 * pow(self._tsi, 4)
        coef1 = pow(coef / psisq, 3.5)
        c2 = coef1 * self._xn0dp * (self._a0dp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + .75
            * CK2 * self._tsi / psisq * self._x3thm1 * (8 + 3 * etasq * (8 + etasq)))
        self._c1 = tle.bStar() * c2
        self._sini0 = sin(self._me.xincl)
        a30vk2 = -XJ3 / CK2 * pow(AE, 3)
        c3 = coef * self._tsi * a30vk2 * self._xn0dp * AE * self._sini0 / self._me.e0
        self._x1mth2 = 1 - theta2
        self._c4 = 2 * self._xn0dp * coef1 * self._a0dp * beta02 * (self._eta
                * (2 + 0.5 * etasq) + self._me.e0 * (0.5 + 2 * etasq) - 2 * CK2 * self._tsi
                / (self._a0dp * psisq) * (-3 * self._x3thm1 * (1 - 2 * eeta + etasq
                * (1.5 - 0.5 * eeta)) + 0.75 * self._x1mth2 * (2 * etasq - eeta
                * (1 + etasq)) * cos(2 * self._me.omega0)))
        self._c5 = 2 * coef1 * self._a0dp * beta02 * (1 + 2.75 * (etasq + eeta) + eeta * etasq)
        theta4 = theta2 * theta2
        temp1 = 3 * CK2 * pinvsq * self._xn0dp
        temp2 = temp1 * CK4 * pinvsq
        temp3 = 1.25 * CK4 * pinvsq * pinvsq * self._xn0dp
        self._xmdot = self._xn0dp + 0.5 * temp1 * beta0 * self._x3thm1 + 0.0625 * temp2 * beta0 \
                * (13 - 78 * theta2 + 137 * theta4)
        x1m5th = 1 - 5 * theta2
        self._omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7 - 114 * theta2
                + 395 * theta4) + temp3 * (3 - 36 * theta2 + 49 * theta4)
        xhdot1 = -temp1 * self._cosi0
        self._xn0dot = xhdot1 + (0.5 * temp2 * (4 - 19 * theta2) + 2 * temp3 * (3
                - 7 * theta2)) * self._cosi0
        self._omgcof = tle.bStar() * c3 * cos(self._me.omega0)
        self._xmcof = -TOTHRD * coef * tle.bStar() * AE / eeta
        self._xnodcf = 3.5 * beta02 * xhdot1 * self._c1
        self._t2cof = 1.5 * self._c1
        self._xlcof = 0.125 * a30vk2 * self._sini0 * (3 + 5 * self._cosi0) / (1 + self._cosi0)
        self._aycof = 0.25 * a30vk2 * self._sini0
        self._delm0 = pow(1 + self._eta * cos(self._me.xm0), 3)
        self._sinm0 = sin(self._me.xm0)
        self._x7thm1 = 7 * theta2 - 1
        if not self._isImp:
            c1sq = self._c1 * self._c1
            self._d2 = 4 * self._a0dp * self._tsi * c1sq
            temp = self._d2 * self._tsi * self._c1 / 3.0
            self._d3 = (17 * self._a0dp + self._s4) * temp
            self._d4 = 0.5 * temp * self._a0dp * self._tsi * (221 * self._a0dp + 31 * self._s4) * self._c1
            self._t3cof = self._d2 + 2 * c1sq
            self._t4cof = 0.25 * (3 * self._d3 + self._c1 * (12 * self._d2 + 10 * c1sq))
            self._t5cof = 0.2 * (3 * self._d4 + 12 * self._c1 * self._d3 + 6 * self._d2 * self._d2 + 15 * c1sq
                                 * (2 * self._d2 + c1sq))


# for now this returns a list of vectors for state vecots,
# we might want to implement a statevector object
class _SGP4_Propagator(_SGP4_Independent):

    def __init__(self, tle: TwoLineElement):
        super().__init__(tle)
        self._tle = tle
        self._tleEpoch = tle.epoch()

    def getState(self, jd: JulianDate) -> tuple[EVector, EVector]:
        dt = jd.difference(self._tleEpoch) * XMNPDA

        #   update for secular gravity and atmospheric drag
        xmdf = self._me.xm0 + self._xmdot * dt
        omgadf = self._me.omega0 + self._omgdot * dt
        xnoddf = self._me.xnode0 + self._xn0dot * dt
        omega = omgadf
        xmp = xmdf
        tsq = dt * dt
        xnode = xnoddf + self._xnodcf * tsq
        tempa = 1 - self._c1 * dt
        tempe = self._tle.bStar() * self._c4 * dt
        templ = self._t2cof * tsq
        if not self._isImp:
            delomg = self._omgcof * dt
            delm = self._xmcof * (pow( 1 + self._eta * cos(xmdf), 3) - self._delm0)
            temp = delomg + delm
            xmp = xmdf + temp
            omega = omgadf - temp
            tcube = tsq * dt
            tfour = dt * tcube
            tempa = tempa - self._d2 * tsq - self._d3 * tcube - self._d4 * tfour
            tempe = tempe + self._tle.bStar() * self._c5 * (sin(xmp) - self._sinm0)
            templ = templ + self._t3cof * tcube + tfour * (self._t4cof + dt * self._t5cof)
        a = self._a0dp * pow(tempa, 2)
        e = self._me.e0 - tempe
        xl = xmp + omega + xnode + self._xn0dp * templ
        beta = sqrt(1 - e * e)
        xn = XKE / pow(a, 1.5)

        # long period periodics
        axn = e * cos(omega)
        temp = 1 / (a * beta * beta)
        xll = temp * self._xlcof * axn
        aynl = temp * self._aycof
        xlt = xl + xll
        ayn = e * sin(omega) + aynl

        # solve keplers equation
        # WANT TO TRY DOING THIS WITH OUT KEPLER SOLVER
        capu = (xlt - xnode) % TWOPI
        if capu < 0:
            capu += TWOPI
        temp2 = capu
        sinepw = 0
        cosepw = 0
        temp3 = 0
        temp4 = 0
        temp5 = 0
        temp6 = 0
        for i in range(10):
            sinepw = sin(temp2)
            cosepw = cos(temp2)
            temp3 = axn * sinepw
            temp4 = ayn * cosepw
            temp5 = axn * cosepw
            temp6 = ayn * sinepw
            epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2
            if abs(epw - temp2) <= E6A:
                break
            temp2 = epw

        # short period preliminary quantities
        ecose = temp5 + temp6
        esine = temp3 - temp4
        elsq = axn * axn + ayn * ayn
        temp = 1 - elsq
        pl = a * temp
        r = a * (1 - ecose)
        temp1 = 1.0 / r
        rdot = XKE * sqrt(a) * esine * temp1
        rfdot = XKE * sqrt(pl) * temp1
        temp2 = a * temp1
        betal = sqrt(temp)
        temp3 = 1.0 / (1 + betal)
        cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
        sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
        u = atan2(sinu, cosu)
        sin2u = 2 * sinu * cosu
        cos2u = 2 * cosu * cosu - 1
        temp = 1 / pl
        temp1 = CK2 * temp
        temp2 = temp1 * temp

        # update for short periodics
        rk = r * (1 - 1.5 * temp2 * betal * self._x3thm1) + 0.5 * temp1 * self._x1mth2 * cos2u
        uk = u - 0.25 * temp2 * self._x7thm1 * sin2u
        xnodek = xnode + 1.5 * temp2 * self._cosi0 * sin2u
        xinck = self._me.xincl + 1.5 * temp2 * self._cosi0 * self._sini0 * cos2u
        rdotk = rdot - xn * temp1 * self._x1mth2 * sin2u
        rfdotk = rfdot + xn * temp1 * (self._x1mth2 * cos2u + 1.5 * self._x3thm1)

        # orientation vectors
        sinuk = sin(uk)
        cosuk = cos(uk)
        sinik = sin(xinck)
        cosik = cos(xinck)
        sinnok = sin(xnodek)
        cosnok = cos(xnodek)
        xmx = -sinnok * cosik
        xmy = cosnok * cosik
        ux = xmx * sinuk + cosnok * cosuk
        uy = xmy * sinuk + sinnok * cosuk
        uz = sinik * sinuk
        vx = xmx * cosuk - cosnok * sinuk
        vy = xmy * cosuk - sinnok * sinuk
        vz = sinik * cosuk

        # position and velocity
        position = EVector(rk * ux, rk * uy, rk * uz) * 1000.0 * XKMPER
        velocity = EVector(
            rdotk * ux + rfdotk * vx,
            rdotk * uy + rfdotk * vy,
            rdotk * uz + rfdotk * vz
        ) * 1000.0 * XKMPER / 60.0
        return (position, velocity)
