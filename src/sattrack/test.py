from math import acos

from numpy import arange
from pyevspace import norm, dot

from sattrack.rotation.order import ZXZ
from sattrack.rotation.rotation import rotateOrderTo
from sattrack.satpass import nextPassMax, PassController
from sattrack.structures.coordinates import *
from sattrack.structures.elements import OrbitalElements
from sattrack.structures.satellite import Satellite
from sattrack.structures.tle import *
from sattrack.sun import getSunPosition, SunPositionController2
from sattrack.topos import getAltitude
from sattrack.util.anomalies import trueToMean
from sattrack.util.constants import SUN_RADIUS, TWOPI
from sattrack.util.conversions import smaToMeanMotion

from sattrack.spacetime.juliandate import now
from sattrack.structures.tle import getTle

# tle = getTle('zarya')
# tle = TwoLineElement("""ISS (ZARYA)
# 1 25544U 98067A   22266.84431519  .00008111  00000+0  14870-3 0  9996
# 2 25544  51.6423 207.8056 0002412 286.8120 181.5821 15.50238875360488""")
tleStr = """STARLINK-1332
1 45579U 20025BA  22267.40990775  .00000340  00000+0  41736-4 0  9990
2 45579  53.0552 303.3671 0001112  89.5038 270.6079 15.06396587135812"""
tle = TwoLineElement(tleStr)
time = JulianDate(9, 25, 2022, 11, 47, 15.468)
iss = Satellite(tle)
jd = now()
geo = GeoPosition(38.0608, -97.9298)
pc = PassController(iss, geo, jd)
np = pc.getNextPass()
plist = pc.getPassList(1)
# elements = OrbitalElements.fromTle(tle, jd)

# def foo(time: JulianDate):
#     tzOffset = time.getTimeZone() / 24.0
#     localVal = time.value() + tzOffset
#     if localVal - int(localVal) < 0.5:
#         # num = (int(localVal) - 0.5)
#         num = int(localVal) - 1
#     # elif localVal - int(localVal) > 0.5:
#     else:
#         # num = (int(localVal) + 0.5)
#         num = int(localVal)
#     return JulianDate.fromNumber(num, 0.5 - tzOffset, time.getTimeZone())
#
# time = JulianDate(9, 29, 2022, 0, 0, 0, -5)
# target = -0.8333333333
# DELTAT = 72.6
#
# sc = SunPositionController2()
# # figure the details here later, for now use middle of the day, not noon
# # A.2.1
# # jd = JulianDate.fromNumber(time.number())
# jd = time
# v = sc.getApparentSiderealTime(jd)
# # A.2.2
# dt = DELTAT / 86400.0
# time_m1 = jd.future(dt - 1)
# time_0 = jd.future(dt)
# time_p1 = jd.future(dt + 1)
# alpha_m1 = sc.getSolarRightAscension(time_m1)
# alpha_0 = sc.getSolarRightAscension(time_0)
# alpha_p1 = sc.getSolarRightAscension(time_p1)
# delta_m1 = sc.getSolarDeclination(time_m1)
# delta_0 = sc.getSolarDeclination(time_0)
# delta_p1 = sc.getSolarDeclination(time_p1)
# # A.2.3
# sigma = radians(geo.getLongitude())
# m0 = (alpha_0 - sigma - v) / TWOPI
# # A.2.4
# hp0 = radians(target)
# phi = radians(geo.getLatitude())
# H0 = acos((sin(hp0) - sin(phi)*sin(delta_0)) / (cos(phi) * cos(delta_0))) % pi
# # A.2.5
# m1 = (m0 - (H0 / TWOPI)) % 1.0
# # A.2.6
# m2 = (m0 + (H0 / TWOPI)) % 1.0
# # A.2.8
# v0 = v + radians(360.985647) * m0
# v1 = v + radians(360.985647) * m1
# v2 = v + radians(360.985647) * m2
# # A.2.9
# n0 = m0 + dt
# n1 = m1 + dt
# n2 = m2 + dt
# # A.2.10
# a = degrees(alpha_0 - alpha_m1)
# if abs(a) > 2:
#     a %= 1.0
# a = radians(a)
# b = degrees(alpha_p1 - alpha_0)
# if abs(b) > 2:
#     b %= 1.0
# b = radians(b)
# ap = degrees(delta_0 - delta_m1)
# if abs(ap) > 2:
#     ap %= 1.0
# ap = radians(ap)
# bp = degrees(delta_p1 - delta_0)
# if abs(bp) > 2:
#     bp %= 1.0
# bp = radians(bp)
# c = b - a
# cp = bp - ap
# alphap0 = alpha_0 + ((n0 * (a + b + c * n0)) / 2.0)
# alphap1 = alpha_0 + ((n1 * (a + b + c * n1)) / 2.0)
# alphap2 = alpha_0 + ((n2 * (a + b + c * n2)) / 2.0)
# deltap0 = delta_0 + ((n0 * (ap + bp + cp * n0)) / 2.0)
# deltap1 = delta_0 + ((n1 * (ap + bp + cp * n1)) / 2.0)
# deltap2 = delta_0 + ((n2 * (ap + bp + cp * n2)) / 2.0)
# # A.2.11
# Hp0 = (v0 + sigma - alphap0) % TWOPI
# if Hp0 >= pi:
#     Hp0 += -TWOPI
# Hp1 = (v1 + sigma - alphap1) % TWOPI
# if Hp1 >= pi:
#     Hp1 += -TWOPI
# Hp2 = (v2 + sigma - alphap2) % TWOPI
# if Hp2 >= pi:
#     Hp2 += -TWOPI
# # A.2.12
# h0 = asin(sin(phi)*sin(deltap0) + cos(phi)*cos(deltap0)*cos(Hp0)) # sun altitude at transit time
# h1 = asin(sin(phi)*sin(deltap1) + cos(phi)*cos(deltap1)*cos(Hp1))
# h2 = asin(sin(phi)*sin(deltap2) + cos(phi)*cos(deltap2)*cos(Hp2))
# # A.2.13
# T = m0 - (Hp0 / TWOPI)
# # A.2.14
# R = m1 + ((h1 - hp0) / (TWOPI*cos(deltap1)*cos(phi)*sin(Hp1)))
# # A.2.15
# S = m2 + ((h2 - hp0) / (TWOPI*cos(deltap2)*cos(phi)*sin(Hp2)))
