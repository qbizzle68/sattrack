from enum import Enum
from functools import total_ordering
from math import sin, cos, radians, atan2, tan, asin, atan, pi, acos, sqrt, degrees

from pyevspace import Vector

from sattrack.util.constants import TWOPI, AU
from sattrack.util.conversions import atan3

# DELTAT = 72.6


# class Variable:
#     # self computing variable
#     __slots__ = '_args', '_value', '_callback', '_parent'
#
#     def __init__(self, callback, args, parent):
#         self._value = None
#         self._callback = callback
#         self._args = args
#         self._parent = parent
#
#     @property
#     def value(self):
#         if self._value is None:
#             arguments = (self._parent[arg] for arg in self._args)
#             self._value = self._callback(*arguments)
#         return self._value
#
#
# class Constant:
#     # class with the same interface as sampa.Variable but has constant value and doesn't need a callback
#     __slots__ = '_value'
#
#     def __init__(self, value):
#         self._value = value
#
#     @property
#     def value(self):
#         return self._value


# class Register:
#     __slots__ = '_internalState'
#
#     def __init__(self, jd):
#         # self._geo = geo
#         JD, JDE, JC, JCE, JME = _generate_times(jd)
#         # list all variables involved in the sampa
#         self._internalState = {
#             'JD': JD, 'JDE': JDE, 'JC': JC, 'JCE': JCE, 'JME': JME,
#             # 'geoLatitude': Constant(radians(geo.latitude)), 'geoLongitude': Constant(radians(geo.longitude)),
#             # 'elevation': Constant(geo.elevation),
#             # MPA: 3.2
#             'moonMeanLongitude': Variable(_mpa_moon_mean_longitude, ('JCE',), self),
#             'moonMeanElongation': Variable(_mpa_moon_mean_elongation, ('JCE',), self),
#             'sunMeanAnomaly': Variable(_mpa_sun_mean_anomaly, ('JCE',), self),
#             'moonMeanAnomaly': Variable(_mpa_moon_mean_anomaly, ('JCE',), self),
#             'moonArgumentLatitude': Variable(_mpa_moon_argument_latitude, ('JCE',), self),
#             'ETerm': Variable(_mpa_E_term, ('JCE',), self),
#             'lrProductTable': Variable(_mpa_lr_table, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
#                                                        'moonArgumentLatitude', 'ETerm'), self),
#             'lTerm': Variable(_mpa_l_term, ('lrProductTable',), self),
#             'rTerm': Variable(_mpa_r_term, ('lrProductTable',), self),
#             'bTerm': Variable(_mpa_b_term, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
#                                             'moonArgumentLatitude', 'ETerm'), self),
#             'a1Term': Variable(_mpa_a1_term, ('JCE',), self),
#             'a2Term': Variable(_mpa_a2_term, ('JCE',), self),
#             'a3Term': Variable(_mpa_a3_term, ('JCE',), self),
#             'deltal': Variable(_mpa_delta_l, ('a1Term', 'a2Term', 'moonMeanLongitude', 'moonArgumentLatitude'), self),
#             'deltab': Variable(_mpa_delta_b, ('a1Term', 'a3Term', 'moonMeanLongitude', 'moonMeanAnomaly',
#                                               'moonArgumentLatitude'), self),
#             'moonLongitude': Variable(_mpa_moon_longitude, ('moonMeanLongitude', 'lTerm', 'deltal'), self),
#             'moonLatitude': Variable(_mpa_moon_latitude, ('bTerm', 'deltab'), self),
#             'moonDistance': Variable(_mpa_moon_distance, ('rTerm',), self),
#             'moonParallax': Variable(_mpa_moon_parallax, ('moonDistance',), self),
#             'xValues': Variable(_x_values, ('JCE',), self),
#             'xyProductTable': Variable(_xy_table, ('xValues',), self),
#             'nutationLongitude': Variable(_nutation_longitude, ('JCE', 'xyProductTable'), self),
#             'nutationObliquity': Variable(_nutation_obliquity, ('JCE', 'xyProductTable'), self),
#             'meanObliquity': Variable(_mean_obliquity, ('JME',), self),
#             'trueObliquity': Variable(_true_obliquity, ('meanObliquity', 'nutationObliquity'), self),
#             'apparentMoonLongitude': Variable(_mpa_apparent_moon_longitude, ('moonLongitude', 'nutationLongitude'),
#                                               self),
#             'meanSiderealTime': Variable(_mean_sidereal_time, ('JD', 'JC'), self),
#             'apparentSiderealTime': Variable(_apparent_sidereal_time, ('meanSiderealTime', 'nutationLongitude',
#                                                                        'nutationObliquity'), self),
#             'moonRightAscension': Variable(_right_ascension, ('moonLongitude', 'moonLatitude',
#                                                               'trueObliquity'), self),
#             'moonDeclination': Variable(_declination, ('moonLongitude', 'moonLatitude', 'trueObliquity'),
#                                         self),
#             'earthHeliocentricLongitude': Variable(_spa_earth_heliocentric_longitude, ('JME',), self),
#             'earthHeliocentricLatitude': Variable(_spa_earth_heliocentric_latitude, ('JME',), self),
#             'sunDistance': Variable(_spa_earth_heliocentric_radius, ('JME',), self),
#             'sunLongitude': Variable(_spa_geocentric_longitude, ('earthHeliocentricLongitude',), self),
#             'sunLatitude': Variable(_spa_geocentric_latitude, ('earthHeliocentricLatitude',), self),
#             'sunAberrationCorrection': Variable(_spa_aberration_correction, ('sunDistance',), self),
#             'apparentSunLongitude': Variable(_spa_apparent_sun_longitude, ('sunLongitude', 'nutationLongitude',
#                                                                            'sunAberrationCorrection'), self),
#             'sunRightAscension': Variable(_right_ascension, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
#             'sunDeclination': Variable(_declination, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
#             'equationOfTime': Variable(_equation_of_time,
#                                        ('JME', 'sunRightAscension', 'nutationLongitude', 'nutationObliquity'), self),
#         }
#
#     def __getitem__(self, item):
#         return self._internalState[item].value


# class RegisterTopocentric(Register):
#
#     def __init__(self, jd, geo):
#         super().__init__(jd)
#         self._internalState['geoLatitude'] = Constant(radians(geo.latitude))
#         self._internalState['geoLongitude'] = Constant(radians(geo.longitude))
#         self._internalState['elevation'] = Constant(geo.elevation)
#         # todo: implement these correctly when able to
#         self._internalState['pressure'] = Constant(0)
#         self._internalState['temperature'] = Constant(10)
#         self._internalState['slopeSurface'] = Constant(0)
#         self._internalState['surfaceAzimuth'] = Constant(0)
#         self._internalState['moonLocalHourAngle'] = Variable(_local_hour_angle, ('apparentSiderealTime', 'geoLongitude',
#                                                                                  'moonRightAscension'), self)
#         self._internalState['uTerm'] = Variable(_u_term, ('geoLatitude',), self)
#         self._internalState['xTerm'] = Variable(_x_term, ('uTerm', 'geoLatitude', 'elevation'), self)
#         self._internalState['yTerm'] = Variable(_y_term, ('uTerm', 'geoLatitude', 'elevation'), self)
#         self._internalState['moonParallaxRightAscension'] = Variable(_parallax_right_ascension,
#                                                                      ('xTerm', 'moonParallax', 'moonLocalHourAngle',
#                                                                       'moonDeclination'), self)
#         self._internalState['moonTopocentricRightAscension'] = Variable(_topocentric_right_ascension,
#                                                                         ('moonRightAscension',
#                                                                          'moonParallaxRightAscension'), self)
#         self._internalState['moonTopocentricDeclination'] = Variable(_topocentric_declination,
#                                                                      ('yTerm', 'moonDeclination', 'moonParallax',
#                                                                       'moonParallaxRightAscension',
#                                                                       'moonLocalHourAngle'), self)
#         self._internalState['sunLocalHourAngle'] = Variable(_local_hour_angle, ('apparentSiderealTime', 'geoLongitude',
#                                                                                 'sunRightAscension'), self)
#         self._internalState['sunParallax'] = Variable(_spa_equitorial_parallax_sun, ('sunDistance',), self)
#         self._internalState['sunParallaxRightAscension'] = Variable(_parallax_right_ascension,
#                                                                     ('xTerm', 'sunParallax', 'sunLocalHourAngle',
#                                                                      'sunDeclination'), self)
#         self._internalState['sunTopocentricRightAscension'] = Variable(_topocentric_right_ascension,
#                                                                        ('sunRightAscension',
#                                                                         'sunParallaxRightAscension'), self)
#         self._internalState['sunTopocentricDeclination'] = Variable(_topocentric_declination,
#                                                                     ('yTerm', 'sunDeclination', 'sunParallax',
#                                                                      'sunParallaxRightAscension', 'sunLocalHourAngle'),
#                                                                     self)
#         self._internalState['moonTopocentricHourAngle'] = Variable(_topocentric_local_hour_angle,
#                                                                    ('moonLocalHourAngle', 'moonParallaxRightAscension'),
#                                                                    self)
#         self._internalState['sunTopocentricHourAngle'] = Variable(_topocentric_local_hour_angle,
#                                                                   ('sunLocalHourAngle', 'sunParallaxRightAscension'),
#                                                                   self)
#         self._internalState['moonElevationAngleWithout'] = Variable(_topocentric_elevation_angle_without,
#                                                                     ('geoLatitude', 'moonTopocentricDeclination',
#                                                                      'moonTopocentricLocalHour'), self)
#         self._internalState['sunElevationAngleWithout'] = Variable(_topocentric_elevation_angle_without,
#                                                                    ('geoLatitude', 'sunTopocentricDeclination',
#                                                                     'sunTopocentricHourAngle'), self)
#         self._internalState['moonAtmosphericRefraction'] = Variable(_atmospheric_refraction_correction,
#                                                                     ('pressure', 'temperature',
#                                                                      'moonElevationAngleWithout'), self)
#         self._internalState['sunAtmosphericRefraction'] = Variable(_atmospheric_refraction_correction,
#                                                                    ('pressure', 'temperature',
#                                                                     'sunElevationAngleWithout'), self)
#         self._internalState['moonElevationAngle'] = Variable(_topocentric_elevation_angle,
#                                                              ('moonElevationAngleWithout', 'moonAtmosphericRefraction'),
#                                                              self)
#         self._internalState['sunElevationAngle'] = Variable(_topocentric_elevation_angle,
#                                                             ('sunElevationAngleWithout', 'sunAtmosphericRefraction'),
#                                                             self)
#         self._internalState['moonZenithAngle'] = Variable(_topocentric_zenith_angle, ('moonElevationAngle',), self)
#         self._internalState['sunZenithAngle'] = Variable(_topocentric_zenith_angle, ('sunElevationAngle',), self)
#         self._internalState['moonAstronomersAzimuth'] = Variable(_topocentric_astronomers_azimuth_angle,
#                                                                  ('moonTopocentricHourAngle', 'geoLatitude',
#                                                                   'moonTopocentricDeclination'), self)
#         self._internalState['sunAstronomersAzimuth'] = Variable(_topocentric_astronomers_azimuth_angle,
#                                                                 ('sunTopocentricHourAngle', 'geoLatitude',
#                                                                  'sunTopocentricDeclination'), self)
#         self._internalState['moonAzimuthAngle'] = Variable(_topocentric_azimuth_angle, ('moonAstronomersAzimuth',),
#                                                            self)
#         self._internalState['sunAzimuthAngle'] = Variable(_topocentric_azimuth_angle, ('sunAstronomersAzimuth',),
#                                                           self)
#         self._internalState['sunIncidenceAngle'] = Variable(_spa_incidence_angle,
#                                                             ('sunZenithAngle', 'sunAstronomersAzimuth', 'slopeSurface',
#                                                              'surfaceAzimuth'), self)


'''
     ------------------------------------------------------------------------------------
    |       description of the variable names and their relation to the SPA and MPA      |
     ------------------------------- -------- ----------- ---------------- --------------
    |     Register Variable Name    | Symbol | Algorithm |    Section     |    Units     |
     ------------------------------- -------- ----------- ---------------- --------------
    |      moonMeanLongitude        |   L'   |    MPA    |     3.2.1      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |      moonMeanElongation       |   D    |    MPA    |     3.2.2      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        sunMeanAnomaly         |   M    |    MPA    |     3.2.3      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       moonMeanAnomaly         |   M'   |    MPA    |     3.2.4      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    moonArgumentLatitude       |   F    |    MPA    |     3.2.5      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |           ETerm               |   E    |    MPA    |    3.2.[6-8]   |   unit-less  |
     ------------------------------- -------- ----------- ---------------- --------------
    |       lrProductTable          |  ----  |    MPA    |    3.2.[6-7]   |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            lTerm              |   l    |    MPA    |     3.2.6      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            rTerm              |   r    |    MPA    |     3.2.7      |  kilometers  |
     ------------------------------- -------- ----------- ---------------- --------------
    |            bTerm              |   b    |    MPA    |     3.2.8      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            a1Term             |   a1   |    MPA    |     3.2.8*     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            a2Term             |   a2   |    MPA    |     3.2.9      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            a3Term             |   a3   |    MPA    |     3.2.10     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            deltal             |   Δl   |    MPA    |     3.2.11     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            deltab             |   Δb   |    MPA    |     3.2.12     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        moonLongitude          |   λ'   |    MPA    |     3.2.13     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         moonLatitude          |   β    |    MPA    |     3.2.14     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         moonDistance          |   Δ    |    MPA    |     3.2.15     |  kilometers  |
     ------------------------------- -------- ----------- ---------------- --------------
    |         moonParallax          |   π    |    MPA    |      3.3       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |           xValues             |   Xi   |  MPA, SPA |    3.4[1-5]    |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        xyProductTable         | Xi*Yij |  MPA, SPA |     3.4.6      | 1e-6 radians |
     ------------------------------- -------- ----------- ---------------- --------------
    |       nutationLongitude       |   Δψ   |  MPA, SPA |     3.4.7      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       nutationObliquity       |   Δε   |  MPA, SPA |     3.4.8      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         meanObliquity         |   ε0   |  MPA, SPA |     3.5.1      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         trueObliquity         |   ε    |  MPA, SPA |     3.5.2      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |     apparentMoonLongitude     |   λ    |    MPA    |      3.6       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       meanSiderealTime        |   ν0   |    MPA    |     3.7.1      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |     apparentSiderealTime      |   ν    |    MPA    |     3.7.2      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |      moonRightAscension       |   α    |    MPA    |      3.8       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       moonDeclination         |   δ    |    MPA    |      3.9       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |      moonLocalHourAngle       |   H    |    MPA    |      3.10      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            uTerm              |   u    |  MPA, SPA | 3.11.1, 3.12.2 |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            xTerm              |   x    |  MPA, SPA | 3.11.2, 3.12.3 |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |            yTerm              |   y    |  MPA, SPA | 3.11.3, 3.12.4 |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   moonParallaxRightAscension  |   Δα   |    MPA    |     3.11.4     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    | moonTopocentricRightAscension |   α'   |    MPA    |     3.11.5     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |  moonTopocentricDeclination   |   δ'   |    MPA    |     3.11.6     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |  earthHeliocentricLongitude   |   L    |    SPA    |    3.2.[4-6]   |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   earthHeliocentricLatitude   |   B    |    SPA    |     3.2.7      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         sunDistance           |   R    |    SPA    |     3.2.8      |      AU      |
     ------------------------------- -------- ----------- ---------------- --------------
    |         sunLongitude          |   Θ    |    SPA    |     3.3.1      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         sunLatitude           |   β    |    SPA    |     3.3.3      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    sunAberrationCorrection    |   Δτ   |    SPA    |      3.6       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |      apparentSunLongitude     |   λ    |    SPA    |      3.7       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       sunRightAscension       |   α    |    SPA    |      3.9       |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        sunDeclination         |   δ    |    SPA    |      3.10      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       sunLocalHourAngle       |   H    |    SPA    |      3.11      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |          sunParallax          |   ξ    |    SPA    |     3.12.1     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   sunParallaxRightAscension   |   Δα   |    SPA    |     3.12.5     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |  sunTopocentricRightAscension |   α'   |    SPA    |     3.12.6     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   sunTopocentricDeclination   |   δ'   |    SPA    |     3.12.7     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    moonTopocentricHourAngle   |   H'   |    MPA    |      3.12      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   moonElevationAngleWithout   |   e0   |    MPA    |     3.13.1     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |   moonAtmosphericRefraction   |   Δe   |    MPA    |     3.13.2     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |      moonElevationAngle       |   e    |    MPA    |     3.13.3     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        moonZenithAngle        |   θm   |    MPA    |     3.13.4     |   radians    
     ------------------------------- -------- ----------- ---------------- --------------
    |    moonAstronomersAzimuth     |   Γ    |    MPA    |     3.14.1     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       moonAzimuthAngle        |   Φm   |    MPA    |     3.14.2     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    sunTopocentricHourAngle    |   H'   |    SPA    |      3.13      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    sunElevationAngleWithout   |   e0   |    SPA    |     3.14.1     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |    sunAtmosphericRefraction   |   Δe   |    SPA    |     3.14.2     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |       sunElevationAngle       |   e    |    SPA    |     3.14.3     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         sunZenithAngle        |   θ    |    SPA    |     3.14.4     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |     sunAstronomersAzimuth     |   Γ    |    SPA    |     3.15.1     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         sunAzimuthAngle       |   Φ    |    SPA    |     3.15.2     |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |        sunIncidenceAngle      |   I    |    SPA    |      3.16      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
    |         equationOfTime        |   E    |    SPA    |       A.1      |   radians    |
     ------------------------------- -------- ----------- ---------------- --------------
'''


'''class SampaComputer:
    __slots__ = '_registry'

    def __init__(self, jd: JulianDate = None):
        if jd is not None and not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type')

        if jd is None:
            self._registry = {}
        else:
            self._registry = {jd: Register(jd)}

    def __call__(self, jd: JulianDate, variableStr: str):
        if not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type')

        if jd not in self._registry:
            self._registry[jd] = Register(jd)

        return self._registry[jd][variableStr]

    def __len__(self):
        return len(self._registry)

    def clear(self, jd=None):
        if jd is not None and _check_jd(jd):
            if jd in self._registry:
                self._registry.pop(jd)
        else:
            self._registry = {}

    def getSunPosition(self, jd: JulianDate):
        _check_jd(jd)
        rightAscension = self.__call__(jd, 'sunRightAscension')
        declination = self.__call__(jd, 'sunDeclination')
        sunDistance = self.__call__(jd, 'sunDistance')

        return _celestial_coordinates_to_position_vector(rightAscension, declination, sunDistance)

    def getSunCoordinates(self, jd: JulianDate):
        _check_jd(jd)
        rightAscension = self.__call__(jd, 'sunRightAscension') * 12 / pi
        declination = degrees(self.__call__(jd, 'sunDeclination'))
        return CelestialCoordinates(rightAscension, declination)

    def getSunTopocentricCoordinates(self, jd: JulianDate):
        _check_jd(jd)
        toposRightAscension = self.__call__(jd, 'sunTopocentricRightAscension') * 12 / pi
        toposDeclination = degrees(self.__call__(jd, 'sunTopocentricDeclination'))
        return CelestialCoordinates(toposRightAscension, toposDeclination)

    def getSunPositionTimes(self, jd: JulianDate, target=None):
        pass

    def getPrecisePositionTimes(self, jd: JulianDate, target=None, epsilon=1e-5):
        pass

    def getMoonPosition(self, jd: JulianDate):
        _check_jd(jd)
        rightAscension = self.__call__(jd, 'moonRightAscension')
        declination = self.__call__(jd, 'moonDeclination')
        moonDistance = self.__call__(jd, 'moonDistance')

        return _celestial_coordinates_to_position_vector(rightAscension, declination, moonDistance)

    # broad getters here (position, rise/set time)

    def nutationLongitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'nutationLongitude')

    def nutationObliquity(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'nutationObliquity')

    def meanObliquity(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'meanObliquity')

    def trueObliquity(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'trueObliquity')

    def meanSiderealTime(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'meanSiderealTime')

    def apparentSiderealTime(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'apparentSiderealTime')

    def moonLongitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonLongitude')

    def moonLatitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonLatitude')

    def moonDistance(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonDistance')

    def moonParallax(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonParallax')

    def moonApparentLongitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'apparentMoonLongitude')

    def moonRightAscension(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonRightAscension')

    def moonDeclination(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonDeclination')

    def sunDistance(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunDistance')

    def sunLongitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunLongitude')

    def sunLatitude(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunLatitude')

    def sunAberrationCorrection(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunAberrationCorrection')

    def sunRightAscension(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunRightAscension')

    def sunDeclination(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunDeclination')

    def sunParallax(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunParallax')

    def equationOfTime(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'equationOfTime')


class SampaTopocentricComputer(SampaComputer):
    __slots__ = '_geo'

    def __init__(self, geo: GeoPosition, jd: JulianDate = None):
        if not isinstance(geo, GeoPosition):
            raise TypeError('geo parameter must be GeoPosition type')
        self._geo = geo
        super().__init__(jd)
        if jd is not None:
            self._registry = {jd: RegisterTopocentric(jd, geo)}

    def __call__(self, jd: JulianDate, variableStr: str):
        if not isinstance(jd, JulianDate):
            raise TypeError('jd parameter must be a JulianDate type')

        if jd not in self._registry:
            self._registry[jd] = Register(jd)

        return self._registry[jd][variableStr]

    def moonLocalHourAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonLocalHourAngle')

    def moonTopocentricRightAscension(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonTopocentricRightAscension')

    def moonTopocentricDeclination(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonTopocentricDeclination')

    def moonTopocentricHourAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonTopocentricHourAngle')

    def moonElevationAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonElevationAngle')

    def moonZenithAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonZenithAngle')

    def moonAzimuthAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'moonAzimuthAngle')

    def sunLocalHourAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunLocalHourAngle')

    def sunTopocentricRightAscension(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunTopocentricRightAscension')

    def sunTopocentricDeclination(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunTopocentricDeclination')

    def sunTopocentricHourAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunTopocentricHourAngle')

    def sunElevationAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunElevationAngle')

    def sunZenithAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunZenithAngle')

    def sunAzimuthAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunAzimuthAngle')

    def sunIncidenceAngle(self, jd: JulianDate):
        _check_jd(jd)
        return self.__call__(jd, 'sunIncidenceAngle')'''


@total_ordering
class TwilightType(Enum):
    Day = 0
    Civil = 1
    Nautical = 2
    Astronomical = 3
    Night = 4

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented

    def __eq__(self, other):
        if self.__class__ is other.__class__:
            return self.value == other.value
        return NotImplemented


# def getTwilightType(jd: JulianDate, geo: GeoPosition):
#     # todo: get sun position in the most efficient way
#     computer = SampaComputer(jd)



# def _generate_times(jd: JulianDate):
#     JD = Constant(jd.value)
#     JDE = Constant(JD.value + DELTAT / 86400)
#     JC = Constant((JD.value - 2451545) / 36525)
#     JCE = Constant((JDE.value - 2451545) / 36525)
#     return JD, JDE, JC, JCE, Constant(JCE.value / 10)


# def _check_jd(jd):
#     if not isinstance(jd, JulianDate):
#         raise TypeError('jd parameter must be JulianDate type')
#     return True


def _mpa_moon_mean_longitude(JCE):
    # radians
    longitude = 3.8103408236 \
                + JCE * (8399.709111634 + JCE * (-2.75517675712e-5 + JCE * (3.2390431537e-8 + JCE * -2.6771317176e-10)))
    return longitude % TWOPI


def _mpa_moon_mean_elongation(JCE):
    # radians
    elongation = 5.19846652984 \
                 + JCE * (7771.37714483 + JCE * (
                  -3.284535119328e-5 + JCE * (3.197346706519e-8 + JCE * 1.5436512200896e-10)))
    return elongation % TWOPI


def _mpa_sun_mean_anomaly(JCE):
    # radians
    anomaly = 6.24006012726 + JCE * (628.301955167 + JCE * (-2.68082573106e-6 + JCE * 7.126701723129e-10))
    return anomaly % TWOPI


def _mpa_moon_mean_anomaly(JCE):
    # 3.2.4 (M'), radians
    anomaly = 2.355555636854 \
              + JCE * (8328.691424759 + JCE * (0.000152566211 + JCE * (2.50409511183e-7 + JCE * -1.1863303779189e-9)))
    return anomaly % TWOPI


def _mpa_moon_argument_latitude(JCE):
    # 3.2.5 (F), radians
    return 1.62790515798 \
           + JCE * (8433.46615806 + JCE * (-6.37725855386e-5 + JCE * (-4.949884435605e-9 + JCE * 2.021671533973e-11)))


def _mpa_E_term(JCE):
    # 3.2.6 (E)
    return 1 + JCE * (-0.002516 + JCE * -0.0000074)


def _mpa_lr_table(moonMeanElongation, sunMeanAnomaly, moonMeanAnomaly, moonLatitude, E):
    # todo: not sure if the lr table values are unit-less, if not change to radians
    # parameters are D, M, M', F
    productTable = []
    for di, mi, miPrime, fi, l, r in _MPA_LR_TERM_TABLE:
        i = di * moonMeanElongation + mi * sunMeanAnomaly + miPrime * moonMeanAnomaly + fi * moonLatitude
        lrAdjust = E ** abs(mi)
        productTable.append([i, lrAdjust, l, r])
    return productTable


def _mpa_l_term(productTable):
    # 3.2.6 (l), in radians
    lTerm = 0.0
    for rowNumber, (product, lAdjust, l, _) in enumerate(productTable):
        lTerm += l * lAdjust * sin(product)
    return radians(lTerm / 1000000)


def _mpa_r_term(productTable):
    # 3.2.7 (r), in radians
    rTerm = 0.0
    for rowNumber, (product, rAdjust, _, r) in enumerate(productTable):
        rTerm += r * rAdjust * cos(product)
    return rTerm / 1000


def _mpa_b_term(moonMeanElongation, sunMeanAnomaly, moonMeanAnomaly, moonLatitude, E):
    # 3.2.8 first (b), in radians
    bTerm = 0.0
    for di, mi, miPrime, fi, bi in _MPA_B_TERM_TABLE:
        bAdjust = E ** abs(mi)
        product = di * moonMeanElongation + mi * sunMeanAnomaly + miPrime * moonMeanAnomaly + fi * moonLatitude
        bTerm += bi * bAdjust * sin(product)
    return radians(bTerm / 1000000)


def _mpa_a1_term(JCE):
    # 3.2.8 second
    return 2.0900317792632 + 2.301199165462 * JCE


def _mpa_a2_term(JCE):
    # 3.2.9
    return 0.926595299884 + 8364.7398477329 * JCE


def _mpa_a3_term(JCE):
    # 3.2.10
    return 5.470734540376 + 8399.6847252966 * JCE


def _mpa_delta_l(a1, a2, moonMeanLongitude, moonLatitude):
    # 3.2.11 (deltal) in radians
    # todo: divide each constant by 1000000
    return (69.0801317939 * sin(a1) + 34.2433599241 * sin(moonMeanLongitude - moonLatitude) + 5.550147021342
            * sin(a2)) / 1000000


def _mpa_delta_b(a1, a3, moonMeanLongitude, moonMeanAnomaly, moonLatitude):
    # 3.2.12 (deltab) in radians
    # todo: divide each constant by 1000000
    deltab = -39.00810878207 * sin(moonMeanLongitude) + 6.667157742618 * sin(a3) + 3.05432619099 \
             * (sin(a1 - moonLatitude) + sin(a1 + moonLatitude)) + 2.21656815003 \
             * sin(moonMeanLongitude - moonMeanAnomaly) - 2.00712863979 * sin(moonMeanLongitude + moonMeanAnomaly)
    return deltab / 1000000


def _mpa_moon_longitude(moonMeanLongitude, lTerm, deltal):
    # 3.2.13 (lambda') in radians
    return (moonMeanLongitude + lTerm + deltal) % TWOPI


def _mpa_moon_latitude(bTerm, deltab):
    # 3.2.14 (beta) in radians
    return (bTerm + deltab) % TWOPI


def _mpa_moon_distance(rTerm):
    # 3.2.16 (delta) in kilometers
    return 385000.56 + rTerm


def _mpa_moon_parallax(moonDistance):
    # 3.3 in radians
    return asin(6378.14 / moonDistance)


def _x_values(JCE):
    # 3.4 in radians
    xTerms = []
    for a, b, c, d in _X_TABLE:
        i = a + JCE * (b + JCE * (c + JCE * d))
        xTerms.append(i)
    return xTerms


def _xy_table(xValues):
    productTable = []
    # for row in _Nutation_Table:
    for row in _Y_TABLE:
        iSum = 0
        # for xj, yij in zip(xValues, row[:4]):
        for xj, yij in zip(xValues, row):
            iSum += xj * yij
        productTable.append(iSum)
    return productTable


def _nutation_longitude(JCE, xyTable):
    # 3.4.7 in radians
    dPsi = 0
    for (a, b, *_), xyProduct in zip(_NUTATION_TABLE, xyTable):
        dPsi += (a + b * JCE) * sin(xyProduct)
    return radians(dPsi / 36000000)


def _nutation_obliquity(JCE, xyTable):
    # 3.4.8 in radians
    dEpsilon = 0
    for (_, _, c, d), xyProduct in zip(_NUTATION_TABLE, xyTable):
        dEpsilon += (c + d * JCE) * cos(xyProduct)
    return radians(dEpsilon / 36000000)


def _mean_obliquity(JME):
    # 3.5.1 (e0) in radians
    U = JME / 10
    return 0.4090928042223 + U \
        * (-0.022693789043 + U
           * (-7.5146120571978e-6 + U
              * (0.00969263751958 + U
                 * (-0.00024909726935 + U
                    * (-0.0012104343176 + U
                       * (-0.00018931974247 + U
                          * (3.4518734095e-5 + U
                             * (0.0001351175729 + U * (2.80707121362e-5 + U * 1.187793518718e-5)))))))))


def _true_obliquity(meanObliquity, nutationObliquity):
    # 3.5.2 (epsilon) in radians
    return meanObliquity + nutationObliquity


def _mpa_apparent_moon_longitude(moonLongitude, nutationLongitude):
    # 3.6 (lambda) in radians
    return moonLongitude + nutationLongitude


def _mean_sidereal_time(JD, JC):
    # 3.7.1 (v0) in radians
    return 4.8949612127358 + 6.300388098985 * (JD - 2451545) + JC * JC * (6.770708127139e-6 + JC * -4.50872966157e-10)


def _apparent_sidereal_time(meanSiderealTime, nutationLongitude, nutationObliquity):
    # 3.7.2 (v) in radians
    return (meanSiderealTime + nutationLongitude * cos(nutationObliquity)) % TWOPI


def _right_ascension(longitude, latitude, trueObliquity):
    # 3.8.1 (alpha) in radians (0, TWOPI)
    alpha = atan2(sin(longitude) * cos(trueObliquity) - tan(latitude) * sin(trueObliquity), cos(longitude))
    return alpha if alpha >= 0 else alpha + TWOPI


def _declination(longitude, latitude, trueObliquity):
    # 3.9 (delta) in radians (-90, 90)
    return asin(sin(latitude) * cos(trueObliquity) + cos(latitude) * sin(trueObliquity) * sin(longitude))


def _local_hour_angle(apparentSiderealTime, geoLongitude, rightAscension):
    # 3.10 (H) in radians
    return apparentSiderealTime + geoLongitude - rightAscension


def _u_term(geoLatitude):
    # 3.11.1 (u) in radians
    return atan(0.99664719 * tan(geoLatitude))


def _x_term(uTerm, geoLatitude, elevation):
    # 3.11.2 (x) in radians
    return cos(uTerm) + (elevation / 6378140) * cos(geoLatitude)


def _y_term(uTerm, geoLatitude, elevation):
    # 3.11.3 (y), in radians
    return 0.99664719 * sin(uTerm) + (elevation / 6378140) * sin(geoLatitude)


def _parallax_right_ascension(xTerm, parallax, hourAngle, declination):
    # 3.11.4 (deltaAlpha) in radians
    return atan2(-xTerm * sin(parallax) * sin(hourAngle), cos(declination) - xTerm * sin(parallax) * cos(hourAngle))


def _topocentric_right_ascension(rightAscension, parallax):
    # 3.11.5 (alpha') in radians
    return rightAscension + parallax


def _topocentric_declination(yTerm, declination, parallax, rightAscensionParallax, hourAngle):
    # 3.11.6 (delta') in radians
    return atan2((sin(declination) - yTerm * sin(parallax)) * cos(rightAscensionParallax),
                 cos(declination) - yTerm * sin(parallax) * cos(hourAngle))


def _topocentric_hour_angle(hourAngle, rightAscensionParallax):
    # 3.12 (H') in radians
    return hourAngle - rightAscensionParallax


def _spa_expand_tables(table, JME):
    # 3.2.1 - 3.2.3
    computedTable = []
    for subTable in table:
        rowSum = 0
        for ai, bi, ci in subTable:
            rowSum += ai * cos(bi + ci * JME)
        computedTable.append(rowSum)
    return computedTable


def _spa_sum_expanded_table(table, JME):
    # 3.2.4 - 3.2.5
    expandedTable = _spa_expand_tables(table, JME)
    rowSum = 0
    for power, i in enumerate(expandedTable):
        rowSum += i * JME ** power
    return rowSum


def _spa_earth_heliocentric_longitude(JME):
    # 3.2.6 (L), in radians
    return (_spa_sum_expanded_table(_SPA_L_TABLE, JME) / 1e8) % TWOPI


def _spa_earth_heliocentric_latitude(JME):
    # 3.2.7 (B), in radians
    B = (_spa_sum_expanded_table(_SPA_B_TABLE, JME) / 1e8) % TWOPI
    return B if B < pi else B - TWOPI


def _spa_earth_heliocentric_radius(JME):
    # 3.2.8 (R), in AU
    return _spa_sum_expanded_table(_SPA_R_TABLE, JME) / 1e8


def _spa_geocentric_longitude(helioLongitude):
    # 3.3.1 (THETA), in radians
    return (helioLongitude + pi) % TWOPI


def _spa_geocentric_latitude(helioLatitude):
    # 3.3.3 (beta), in radians
    return -helioLatitude


def _spa_aberration_correction(sunDistance):
    # 3.6 (deltaTau), in radians
    return -0.357614473075 / (3600 * sunDistance)


def _spa_apparent_sun_longitude(geocentricLongitude, nutationLongitude, aberrationCorrection):
    # 3.7 (lambda), in radians
    return geocentricLongitude + nutationLongitude + aberrationCorrection


def _spa_equitorial_parallax_sun(sunDistance):
    # 3.12.1 (ksi), in radians
    return 0.15348425442038 / (3600 * sunDistance)


def _topocentric_local_hour_angle(hourAngle, parallaxRightAscension):
    # 3.12, 3.12.5 (H'), in radians
    # todo: valid range here?
    return hourAngle - parallaxRightAscension


def _topocentric_elevation_angle_without(geoLatitude, topoDeclination, topoHourAngle):
    # 3.13.1, 3.14.1 (e0), in radians
    return asin(sin(geoLatitude)*sin(topoDeclination) + cos(geoLatitude)*cos(topoDeclination)*cos(topoHourAngle))


def _atmospheric_refraction_correction(pressure, temperature, e0):
    # 3.13.2, 3.14.2 (deltae), in radians
    denominator = 60 * tan(e0 + (10.3 / (e0 + 5.11)))
    return (pressure / 1010) * (283 / (273 + temperature)) * (0.01780235837034 / denominator)


def _topocentric_elevation_angle(e0, atmosphericCorrection):
    # 3.13.3, 3.14.3 (e), in radians
    return e0 + atmosphericCorrection


def _topocentric_zenith_angle(elevationAngle):
    # 3.13.4, 3.14.4 (THETAm), in radians
    return pi/2 - elevationAngle


def _topocentric_astronomers_azimuth_angle(topoHourAngle, geoLatitude, topoDeclination):
    # 3.14.1, 3.15.1 (GAMMA), in radians
    return atan3(sin(topoHourAngle), cos(topoHourAngle)*sin(geoLatitude) - tan(topoDeclination)*cos(geoLatitude))


def _topocentric_azimuth_angle(astronomersAngle):
    # 3.14.2, 3.15.2 (PHIm), in radians
    return (astronomersAngle + pi) % TWOPI


def _spa_incidence_angle(zenithAngle, astronomersAzimuth, slopeSurface, surfaceAzimuth):
    # 3.16 (I), in radians
    return acos(cos(zenithAngle)*cos(slopeSurface)
                + sin(slopeSurface)*sin(zenithAngle)*cos(astronomersAzimuth-surfaceAzimuth))


def _equation_of_time(JME, sunRightAscension, nutationLongitude, nutationObliquity):
    # A.1 (E), in radians
    sunMeanLongitude = 4.895063110817 + JME \
                       * (6283.3196674757 + JME
                          * (0.000529188716 + JME
                             * (3.49548226952e-7 + JME * (-1.1407380731989e-6 + JME * -8.726646259972e-9))))
    E = sunMeanLongitude - 9.980316262e-5 - sunRightAscension + nutationLongitude * cos(nutationObliquity)
    return E % TWOPI


def _celestial_coordinates_to_position_vector(rightAscension, declination, sunDistance):
    # radians and kilometers, return kilometers
    zComp = sin(declination)
    xComp = sqrt((cos(declination) ** 2) / (1 + (tan(rightAscension) ** 2)))
    if pi / 2 < rightAscension < 3 * pi / 2:
        xComp = -xComp
    yComp = xComp * tan(rightAscension)
    yComp = abs(yComp) if rightAscension < pi else -abs(yComp)

    return Vector((xComp, yComp, zComp)) * sunDistance * AU


def _get_twilight_type(topocentricPositionVector):
    # vector in kilometers
    sunAngle = degrees(asin(topocentricPositionVector[2] / topocentricPositionVector.mag()))

    if sunAngle < -18:
        return TwilightType.Night
    elif sunAngle < -12:
        return TwilightType.Astronomical
    elif sunAngle < -6:
        return TwilightType.Nautical
    elif sunAngle < -(5.0 / 6.0):
        return TwilightType.Civil
    else:
        return TwilightType.Day


_SPA_L_TABLE = [
    [
        [175347046, 0, 0],
        [3341656, 4.6692568, 6283.07585],
        [34894, 4.6261, 12566.1517],
        [3497, 2.7441, 5753.3849],
        [3418, 2.8289, 3.5231],
        [3136, 3.6277, 77713.7715],
        [2676, 4.4181, 7860.4194],
        [2343, 6.1352, 3930.2097],
        [1324, 0.7425, 11506.7698],
        [1273, 2.0371, 529.691],
        [1199, 1.1096, 1577.3435],
        [990, 5.233, 5884.927],
        [902, 2.045, 26.298],
        [857, 3.508, 398.149],
        [780, 1.179, 5223.694],
        [753, 2.533, 5507.553],
        [505, 4.583, 18849.228],
        [492, 4.205, 775.523],
        [357, 2.92, 0.067],
        [317, 5.849, 11790.629],
        [284, 1.899, 796.298],
        [271, 0.315, 10977.079],
        [243, 0.345, 5486.778],
        [206, 4.806, 2544.314],
        [205, 1.869, 5573.143],
        [202, 2.458, 6069.777],
        [156, 0.833, 213.299],
        [132, 3.411, 2942.463],
        [126, 1.083, 20.775],
        [115, 0.645, 0.98],
        [103, 0.636, 4694.003],
        [102, 0.976, 15720.839],
        [102, 4.267, 7.114],
        [99, 6.21, 2146.17],
        [98, 0.68, 155.42],
        [86, 5.98, 161000.69],
        [85, 1.3, 6275.96],
        [85, 3.67, 71430.7],
        [80, 1.81, 17260.15],
        [79, 3.04, 12036.46],
        [75, 1.76, 5088.63],
        [74, 3.5, 3154.69],
        [74, 4.68, 801.82],
        [70, 0.83, 9437.76],
        [62, 3.98, 8827.39],
        [61, 1.82, 7084.9],
        [57, 2.78, 6286.6],
        [56, 4.39, 14143.5],
        [56, 3.47, 6279.55],
        [52, 0.19, 12139.55],
        [52, 1.33, 1748.02],
        [51, 0.28, 5856.48],
        [49, 0.49, 1194.45],
        [41, 5.37, 8429.24],
        [41, 2.4, 19651.05],
        [39, 6.17, 10447.39],
        [37, 6.04, 10213.29],
        [37, 2.57, 1059.38],
        [36, 1.71, 2352.87],
        [36, 1.78, 6812.77],
        [33, 0.59, 17789.85],
        [30, 0.44, 83996.85],
        [30, 2.74, 1349.87],
        [25, 3.16, 4690.48]
    ],
    [
        [628331966747.0, 0, 0],
        [206059, 2.678235, 6283.07585],
        [4303, 2.6351, 12566.1517],
        [425, 1.59, 3.523],
        [119, 5.769, 26.298],
        [109, 2.966, 1577.344],
        [93, 2.59, 18849.23],
        [72, 1.14, 529.69],
        [68, 1.87, 398.15],
        [67, 4.41, 5507.55],
        [59, 2.89, 5223.69],
        [56, 2.17, 155.42],
        [45, 0.4, 796.3],
        [36, 0.47, 775.52],
        [29, 2.65, 7.11],
        [21, 5.34, 0.98],
        [19, 1.85, 5486.78],
        [19, 4.97, 213.3],
        [17, 2.99, 6275.96],
        [16, 0.03, 2544.31],
        [16, 1.43, 2146.17],
        [15, 1.21, 10977.08],
        [12, 2.83, 1748.02],
        [12, 3.26, 5088.63],
        [12, 5.27, 1194.45],
        [12, 2.08, 4694],
        [11, 0.77, 553.57],
        [10, 1.3, 6286.6],
        [10, 4.24, 1349.87],
        [9, 2.7, 242.73],
        [9, 5.64, 951.72],
        [8, 5.3, 2352.87],
        [6, 2.65, 9437.76],
        [6, 4.67, 4690.48]
    ],
    [
        [52919, 0, 0],
        [8720, 1.0721, 6283.0758],
        [309, 0.867, 12566.152],
        [27, 0.05, 3.52],
        [16, 5.19, 26.3],
        [16, 3.68, 155.42],
        [10, 0.76, 18849.23],
        [9, 2.06, 77713.77],
        [7, 0.83, 775.52],
        [5, 4.66, 1577.34],
        [4, 1.03, 7.11],
        [4, 3.44, 5573.14],
        [3, 5.14, 796.3],
        [3, 6.05, 5507.55],
        [3, 1.19, 242.73],
        [3, 6.12, 529.69],
        [3, 0.31, 398.15],
        [3, 2.28, 553.57],
        [2, 4.38, 5223.69],
        [2, 3.75, 0.98]
    ],
    [
        [289, 5.844, 6283.076],
        [35, 0, 0],
        [17, 5.49, 12566.15],
        [3, 5.2, 155.42],
        [1, 4.72, 3.52],
        [1, 5.3, 18849.23],
        [1, 5.97, 242.73]
    ],
    [
        [114, 3.142, 0],
        [8, 4.13, 6283.08],
        [1, 3.84, 12566.15]
    ],
    [
        [1, 3.14, 0]
    ]
]

_SPA_B_TABLE = [
    [
        [280, 3.199, 84334.662],
        [102, 5.422, 5507.553],
        [80, 3.88, 5223.69],
        [44, 3.7, 2352.87],
        [32, 4, 1577.34]
    ],
    [
        [9, 3.9, 5507.55],
        [6, 1.73, 5223.69]
    ]
]

_SPA_R_TABLE = [
    [
        [100013989, 0, 0],
        [1670700, 3.0984635, 6283.07585],
        [13956, 3.05525, 12566.1517],
        [3084, 5.1985, 77713.7715],
        [1628, 1.1739, 5753.3849],
        [1576, 2.8469, 7860.4194],
        [925, 5.453, 11506.77],
        [542, 4.564, 3930.21],
        [472, 3.661, 5884.927],
        [346, 0.964, 5507.553],
        [329, 5.9, 5223.694],
        [307, 0.299, 5573.143],
        [243, 4.273, 11790.629],
        [212, 5.847, 1577.344],
        [186, 5.022, 10977.079],
        [175, 3.012, 18849.228],
        [110, 5.055, 5486.778],
        [98, 0.89, 6069.78],
        [86, 5.69, 15720.84],
        [86, 1.27, 161000.69],
        [65, 0.27, 7260.15],
        [63, 0.92, 529.69],
        [57, 2.01, 83996.85],
        [56, 5.24, 71430.7],
        [49, 3.25, 2544.31],
        [47, 2.58, 775.52],
        [45, 5.54, 9437.76],
        [43, 6.01, 6275.96],
        [39, 5.36, 4694],
        [38, 2.39, 8827.39],
        [37, 0.83, 19651.05],
        [37, 4.9, 12139.55],
        [36, 1.67, 12036.46],
        [35, 1.84, 2942.46],
        [33, 0.24, 7084.9],
        [32, 0.18, 5088.63],
        [32, 1.78, 398.15],
        [28, 1.21, 6286.6],
        [28, 1.9, 6279.55],
        [26, 4.59, 10447.39]
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
        [4359, 5.7846, 6283.0758],
        [124, 5.579, 12566.152],
        [12, 3.14, 0],
        [9, 3.63, 77713.77],
        [6, 1.87, 5573.14],
        [3, 5.47, 18849.23]
    ],
    [
        [145, 4.273, 6283.076],
        [7, 3.92, 12566.15]
    ],
    [
        [4, 2.56, 6283.08]
    ]
]

_X_TABLE = [
    [5.1984694602504185, 7771.377146170642, -3.340909254167546e-05, 9.211444588673535e-08],
    [6.240036788719593, 628.3019560241842, -2.79776279094691e-06, -5.817764173314432e-08],
    [2.355548369303256, 8328.691422882925, 0.00015179477570445083, 3.10280755910103e-07],
    [1.6279019291238244, 8433.466158317484, -6.427174970469119e-05, 5.332994933829344e-08],
    [2.1824385855759, -33.757045936662394, 3.614227815029857e-05, 3.8785094488762875e-08]
]

_Y_TABLE = [
    [0, 0, 0, 0, 1],
    [-2, 0, 0, 2, 2],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [-2, 1, 0, 2, 2],
    [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2],
    [-2, -1, 0, 2, 2],
    [-2, 0, 1, 0, 0],
    [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2],
    [2, 0, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1],
    [0, 0, 1, 2, 1],
    [-2, 0, 2, 0, 0],
    [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2],
    [0, 0, 2, 2, 2],
    [0, 0, 2, 0, 0],
    [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0],
    [-2, 0, 0, 2, 0],
    [0, 0, -1, 2, 1],
    [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1],
    [-2, 2, 0, 2, 2],
    [0, 1, 0, 0, 1],
    [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1],
    [0, 0, 2, -2, 0],
    [2, 0, -1, 2, 1],
    [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2],
    [-2, 1, 1, 0, 0],
    [0, -1, 0, 2, 2],
    [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0],
    [-2, 0, 2, 2, 2],
    [-2, 0, 1, 2, 1],
    [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1],
    [0, -1, 1, 0, 0],
    [-2, -1, 0, 2, 1],
    [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1],
    [-2, 0, 2, 0, 1],
    [-2, 1, 0, 2, 1],
    [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0],
    [-2, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2],
    [-1, -1, 1, 0, 0],
    [0, 1, 1, 0, 0],
    [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2],
    [0, 0, 3, 2, 2],
    [2, -1, 0, 2, 2]
]

_NUTATION_TABLE = [
    [-171996, -174.2, 92025, 8.9],
    [-13187, -1.6, 5736, -3.1],
    [-2274, -0.2, 977, -0.5],
    [2062, 0.2, -895, 0.5],
    [1426, -3.4, 54, -0.1],
    [712, 0.1, -7, 0],
    [-517, 1.2, 224, -0.6],
    [-386, -0.4, 200, 0],
    [-301, 0, 129, -0.1],
    [217, -0.5, -95, 0.3],
    [-158, 0, 0, 0],
    [129, 0.1, -70, 0],
    [123, 0, -53, 0],
    [63, 0, 0, 0],
    [63, 0.1, -33, 0],
    [-59, 0, 26, 0],
    [-58, -0.1, 32, 0],
    [-51, 0, 27, 0],
    [48, 0, 0, 0],
    [46, 0, -24, 0],
    [-38, 0, 16, 0],
    [-31, 0, 13, 0],
    [29, 0, 0, 0],
    [29, 0, -12, 0],
    [26, 0, 0, 0],
    [-22, 0, 0, 0],
    [21, 0, -10, 0],
    [17, -0.1, 0, 0],
    [16, 0, -8, 0],
    [-16, 0.1, 7, 0],
    [-15, 0, 9, 0],
    [-13, 0, 7, 0],
    [-12, 0, 6, 0],
    [11, 0, 0, 0],
    [-10, 0, 5, 0],
    [-8, 0, 3, 0],
    [7, 0, -3, 0],
    [-7, 0, 0, 0],
    [-7, 0, 3, 0],
    [-7, 0, 3, 0],
    [6, 0, 0, 0],
    [6, 0, -3, 0],
    [6, 0, -3, 0],
    [-6, 0, 3, 0],
    [-6, 0, 3, 0],
    [5, 0, 0, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0]
]

_MPA_LR_TERM_TABLE = [
    # d m  m' f      l         r
    [0, 0, 1, 0, 6288774, -20905355],
    [2, 0, -1, 0, 1274027, -3699111],
    [2, 0, 0, 0, 658314, -2955968],
    [0, 0, 2, 0, 213618, -569925],
    [0, 1, 0, 0, -185116, 48888],
    [0, 0, 0, 2, -114332, -3149],
    [2, 0, -2, 0, 58793, 246158],
    [2, -1, -1, 0, 57066, -152138],
    [2, 0, 1, 0, 53322, -170733],
    [2, -1, 0, 0, 45758, -204586],
    [0, 1, -1, 0, -40923, -129620],
    [1, 0, 0, 0, -34720, 108743],
    [0, 1, 1, 0, -30383, 104755],
    [2, 0, 0, -2, 15327, 10321],
    [0, 0, 1, 2, -12528, 0],
    [0, 0, 1, -2, 10980, 79661],
    [4, 0, -1, 0, 10675, -34782],
    [0, 0, 3, 0, 10034, -23210],
    [4, 0, -2, 0, 8548, -21636],
    [2, 1, -1, 0, -7888, 24208],
    [2, 1, 0, 0, -6766, 30824],
    [1, 0, -1, 0, -5163, -8379],
    [1, 1, 0, 0, 4987, -16675],
    [2, -1, 1, 0, 4036, -12831],
    [2, 0, 2, 0, 3994, -10445],
    [4, 0, 0, 0, 3861, -11650],
    [2, 0, -3, 0, 3665, 14403],
    [0, 1, -2, 0, -2689, -7003],
    [2, 0, -1, 2, -2602, 0],
    [2, -1, -2, 0, 2390, 10056],
    [1, 0, 1, 0, -2348, 6322],
    [2, -2, 0, 0, 2236, -9884],
    [0, 1, 2, 0, -2120, 5751],
    [0, 2, 0, 0, -2069, 0],
    [2, -2, -1, 0, 2048, -4950],
    [2, 0, 1, -2, -1773, 4130],
    [2, 0, 0, 2, -1595, 0],
    [4, -1, -1, 0, 1215, -3958],
    [0, 0, 2, 2, -1110, 0],
    [3, 0, -1, 0, -892, 3258],
    [2, 1, 1, 0, -810, 2616],
    [4, -1, -2, 0, 759, -1897],
    [0, 2, -1, 0, -713, -2117],
    [2, 2, -1, 0, -700, 2354],
    [2, 1, -2, 0, 691, 0],
    [2, -1, 0, -2, 596, 0],
    [4, 0, 1, 0, 549, -1423],
    [0, 0, 4, 0, 537, -1117],
    [4, -1, 0, 0, 520, -1571],
    [1, 0, -2, 0, -487, -1739],
    [2, 1, 0, -2, -399, 0],
    [0, 0, 2, -2, -381, -4421],
    [1, 1, 1, 0, 351, 0],
    [3, 0, -2, 0, -340, 0],
    [4, 0, -3, 0, 330, 0],
    [2, -1, 2, 0, 327, 0],
    [0, 2, 1, 0, -323, 1165],
    [1, 1, -1, 0, 299, 0],
    [2, 0, 3, 0, 294, 0],
    [2, 0, -1, -2, 0, 8752]
]

_MPA_B_TERM_TABLE = [
    # d  m  m' f     b
    [0, 0, 0, 1, 5128122],
    [0, 0, 1, 1, 280602],
    [0, 0, 1, -1, 277693],
    [2, 0, 0, -1, 173237],
    [2, 0, -1, 1, 55413],
    [2, 0, -1, -1, 46271],
    [2, 0, 0, 1, 32573],
    [0, 0, 2, 1, 17198],
    [2, 0, 1, -1, 9266],
    [0, 0, 2, -1, 8822],
    [2, -1, 0, -1, 8216],
    [2, 0, -2, -1, 4324],
    [2, 0, 1, 1, 4200],
    [2, 1, 0, -1, -3359],
    [2, -1, -1, 1, 2463],
    [2, -1, 0, 1, 2211],
    [2, -1, -1, -1, 2065],
    [0, 1, -1, -1, -1870],
    [4, 0, -1, -1, 1828],
    [0, 1, 0, 1, -1794],
    [0, 0, 0, 3, -1749],
    [0, 1, -1, 1, -1565],
    [1, 0, 0, 1, -1491],
    [0, 1, 1, 1, -1475],
    [0, 1, 1, -1, -1410],
    [0, 1, 0, -1, -1344],
    [1, 0, 0, -1, -1335],
    [0, 0, 3, 1, 1107],
    [4, 0, 0, -1, 1021],
    [4, 0, -1, 1, 833],
    [0, 0, 1, -3, 777],
    [4, 0, -2, 1, 671],
    [2, 0, 0, -3, 607],
    [2, 0, 2, -1, 596],
    [2, -1, 1, -1, 491],
    [2, 0, -2, 1, -451],
    [0, 0, 3, -1, 439],
    [2, 0, 2, 1, 422],
    [2, 0, -3, -1, 421],
    [2, 1, -1, 1, -366],
    [2, 1, 0, 1, -351],
    [4, 0, 0, 1, 331],
    [2, -1, 1, 1, 315],
    [2, -2, 0, -1, 302],
    [0, 0, 1, 3, -283],
    [2, 1, 1, -1, -229],
    [1, 1, 0, -1, 223],
    [1, 1, 0, 1, 223],
    [0, 1, -2, -1, -220],
    [2, 1, -1, -1, -220],
    [1, 0, 1, 1, -185],
    [2, -1, -2, -1, 181],
    [0, 1, 2, 1, -177],
    [4, 0, -2, -1, 176],
    [4, -1, -1, -1, 166],
    [1, 0, 1, -1, -164],
    [4, 0, 1, -1, 132],
    [1, 0, -1, -1, -119],
    [4, -1, 0, -1, 115],
    [2, -2, 0, 1, 107]
]
