from sattrack._sampa import *
from sattrack.util.constants import DELTAT

from ._sampa import _mpa_moon_mean_longitude, _mpa_moon_mean_elongation, _mpa_sun_mean_anomaly, \
    _mpa_moon_mean_anomaly, _mpa_moon_argument_latitude, _mpa_E_term, _mpa_lr_table, _mpa_l_term, _mpa_r_term, \
    _mpa_b_term, _mpa_a1_term, _mpa_a2_term, _mpa_a3_term, _mpa_delta_l, _mpa_delta_b, _mpa_moon_longitude, \
    _mpa_moon_latitude, _mpa_moon_distance, _mpa_moon_parallax, _x_values, _xy_table, _nutation_longitude, \
    _nutation_obliquity, _mean_obliquity, _true_obliquity, _mpa_apparent_moon_longitude, _mean_sidereal_time, \
    _apparent_sidereal_time, _right_ascension, _declination, _spa_earth_heliocentric_longitude, \
    _spa_earth_heliocentric_latitude, _spa_geocentric_longitude, _spa_geocentric_latitude, _spa_aberration_correction, \
    _spa_apparent_sun_longitude, _equation_of_time, _spa_earth_heliocentric_radius, _local_hour_angle, _u_term, \
    _x_term, _y_term, _parallax_right_ascension, _topocentric_right_ascension, _topocentric_declination, \
    _spa_equitorial_parallax_sun, _topocentric_local_hour_angle, _topocentric_elevation_angle_without, \
    _atmospheric_refraction_correction, _topocentric_elevation_angle, _topocentric_zenith_angle, \
    _topocentric_astronomers_azimuth_angle, _topocentric_azimuth_angle, _spa_incidence_angle


# todo: put useful values here like true obliquity etc.


''' saving this so we dont lose it in case we want to use it'''
class Variable:
    # self computing variable
    __slots__ = '_args', '_value', '_callback', '_parent'

    def __init__(self, callback, args, parent):
        self._value = None
        self._callback = callback
        self._args = args
        self._parent = parent

    @property
    def value(self):
        if self._value is None:
            arguments = (self._parent[arg] for arg in self._args)
            self._value = self._callback(*arguments)
        return self._value


class Constant:
    # class with the same interface as sampa.Variable but has constant value and doesn't need a callback
    __slots__ = '_value'

    def __init__(self, value):
        self._value = value

    @property
    def value(self):
        return self._value


def _generate_times(jd: JulianDate):
    JD = jd.value
    JDE = JD + DELTAT / 86400
    JC = (JD - 2451545) / 36525
    JCE = (JDE - 2451545) / 36525
    return JD, JDE, JC, JCE, JCE / 10


def _check_jd(jd):
    if not isinstance(jd, JulianDate):
        raise TypeError('jd parameter must be JulianDate type')
    return True


class Register:
    __slots__ = '_internalState'

    def __init__(self, jd):
        # self._geo = geo
        JD, JDE, JC, JCE, JME = _generate_times(jd)
        # list all variables involved in the sampa
        self._internalState = {
            'JD': Constant(JD), 'JDE': Constant(JDE), 'JC': Constant(JC), 'JCE': Constant(JCE), 'JME': Constant(JME),
            # 'geoLatitude': Constant(radians(geo.latitude)), 'geoLongitude': Constant(radians(geo.longitude)),
            # 'elevation': Constant(geo.elevation),
            # MPA: 3.2
            'moonMeanLongitude': Variable(_mpa_moon_mean_longitude, ('JCE',), self),
            'moonMeanElongation': Variable(_mpa_moon_mean_elongation, ('JCE',), self),
            'sunMeanAnomaly': Variable(_mpa_sun_mean_anomaly, ('JCE',), self),
            'moonMeanAnomaly': Variable(_mpa_moon_mean_anomaly, ('JCE',), self),
            'moonArgumentLatitude': Variable(_mpa_moon_argument_latitude, ('JCE',), self),
            'ETerm': Variable(_mpa_E_term, ('JCE',), self),
            'lrProductTable': Variable(_mpa_lr_table, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
                                                       'moonArgumentLatitude', 'ETerm'), self),
            'lTerm': Variable(_mpa_l_term, ('lrProductTable',), self),
            'rTerm': Variable(_mpa_r_term, ('lrProductTable',), self),
            'bTerm': Variable(_mpa_b_term, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
                                            'moonArgumentLatitude', 'ETerm'), self),
            'a1Term': Variable(_mpa_a1_term, ('JCE',), self),
            'a2Term': Variable(_mpa_a2_term, ('JCE',), self),
            'a3Term': Variable(_mpa_a3_term, ('JCE',), self),
            'deltal': Variable(_mpa_delta_l, ('a1Term', 'a2Term', 'moonMeanLongitude', 'moonArgumentLatitude'), self),
            'deltab': Variable(_mpa_delta_b, ('a1Term', 'a3Term', 'moonMeanLongitude', 'moonMeanAnomaly',
                                              'moonArgumentLatitude'), self),
            'moonLongitude': Variable(_mpa_moon_longitude, ('moonMeanLongitude', 'lTerm', 'deltal'), self),
            'moonLatitude': Variable(_mpa_moon_latitude, ('bTerm', 'deltab'), self),
            'moonDistance': Variable(_mpa_moon_distance, ('rTerm',), self),
            'moonParallax': Variable(_mpa_moon_parallax, ('moonDistance',), self),
            'xValues': Variable(_x_values, ('JCE',), self),
            'xyProductTable': Variable(_xy_table, ('xValues',), self),
            'nutationLongitude': Variable(_nutation_longitude, ('JCE', 'xyProductTable'), self),
            'nutationObliquity': Variable(_nutation_obliquity, ('JCE', 'xyProductTable'), self),
            'meanObliquity': Variable(_mean_obliquity, ('JME',), self),
            'trueObliquity': Variable(_true_obliquity, ('meanObliquity', 'nutationObliquity'), self),
            'apparentMoonLongitude': Variable(_mpa_apparent_moon_longitude, ('moonLongitude', 'nutationLongitude'),
                                              self),
            'meanSiderealTime': Variable(_mean_sidereal_time, ('JD', 'JC'), self),
            'apparentSiderealTime': Variable(_apparent_sidereal_time, ('meanSiderealTime', 'nutationLongitude',
                                                                       'nutationObliquity'), self),
            'moonRightAscension': Variable(_right_ascension, ('moonLongitude', 'moonLatitude',
                                                              'trueObliquity'), self),
            'moonDeclination': Variable(_declination, ('moonLongitude', 'moonLatitude', 'trueObliquity'),
                                        self),
            'earthHeliocentricLongitude': Variable(_spa_earth_heliocentric_longitude, ('JME',), self),
            'earthHeliocentricLatitude': Variable(_spa_earth_heliocentric_latitude, ('JME',), self),
            'sunDistance': Variable(_spa_earth_heliocentric_radius, ('JME',), self),
            'sunLongitude': Variable(_spa_geocentric_longitude, ('earthHeliocentricLongitude',), self),
            'sunLatitude': Variable(_spa_geocentric_latitude, ('earthHeliocentricLatitude',), self),
            'sunAberrationCorrection': Variable(_spa_aberration_correction, ('sunDistance',), self),
            'apparentSunLongitude': Variable(_spa_apparent_sun_longitude, ('sunLongitude', 'nutationLongitude',
                                                                           'sunAberrationCorrection'), self),
            'sunRightAscension': Variable(_right_ascension, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
            'sunDeclination': Variable(_declination, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
            'equationOfTime': Variable(_equation_of_time,
                                       ('JME', 'sunRightAscension', 'nutationLongitude', 'nutationObliquity'), self),
        }

    def __getitem__(self, item):
        return self._internalState[item].value


class RegisterTopocentric(Register):

    def __init__(self, jd, geo):
        super().__init__(jd)
        self._internalState['geoLatitude'] = Constant(radians(geo.latitude))
        self._internalState['geoLongitude'] = Constant(radians(geo.longitude))
        self._internalState['elevation'] = Constant(geo.elevation)
        # todo: implement these correctly when able to
        self._internalState['pressure'] = Constant(0)
        self._internalState['temperature'] = Constant(10)
        self._internalState['slopeSurface'] = Constant(0)
        self._internalState['surfaceAzimuth'] = Constant(0)
        self._internalState['moonLocalHourAngle'] = Variable(_local_hour_angle, ('apparentSiderealTime', 'geoLongitude',
                                                                                 'moonRightAscension'), self)
        self._internalState['uTerm'] = Variable(_u_term, ('geoLatitude',), self)
        self._internalState['xTerm'] = Variable(_x_term, ('uTerm', 'geoLatitude', 'elevation'), self)
        self._internalState['yTerm'] = Variable(_y_term, ('uTerm', 'geoLatitude', 'elevation'), self)
        self._internalState['moonParallaxRightAscension'] = Variable(_parallax_right_ascension,
                                                                     ('xTerm', 'moonParallax', 'moonLocalHourAngle',
                                                                      'moonDeclination'), self)
        self._internalState['moonTopocentricRightAscension'] = Variable(_topocentric_right_ascension,
                                                                        ('moonRightAscension',
                                                                         'moonParallaxRightAscension'), self)
        self._internalState['moonTopocentricDeclination'] = Variable(_topocentric_declination,
                                                                     ('yTerm', 'moonDeclination', 'moonParallax',
                                                                      'moonParallaxRightAscension',
                                                                      'moonLocalHourAngle'), self)
        self._internalState['sunLocalHourAngle'] = Variable(_local_hour_angle, ('apparentSiderealTime', 'geoLongitude',
                                                                                'sunRightAscension'), self)
        self._internalState['sunParallax'] = Variable(_spa_equitorial_parallax_sun, ('sunDistance',), self)
        self._internalState['sunParallaxRightAscension'] = Variable(_parallax_right_ascension,
                                                                    ('xTerm', 'sunParallax', 'sunLocalHourAngle',
                                                                     'sunDeclination'), self)
        self._internalState['sunTopocentricRightAscension'] = Variable(_topocentric_right_ascension,
                                                                       ('sunRightAscension',
                                                                        'sunParallaxRightAscension'), self)
        self._internalState['sunTopocentricDeclination'] = Variable(_topocentric_declination,
                                                                    ('yTerm', 'sunDeclination', 'sunParallax',
                                                                     'sunParallaxRightAscension', 'sunLocalHourAngle'),
                                                                    self)
        self._internalState['moonTopocentricHourAngle'] = Variable(_topocentric_local_hour_angle,
                                                                   ('moonLocalHourAngle', 'moonParallaxRightAscension'),
                                                                   self)
        self._internalState['sunTopocentricHourAngle'] = Variable(_topocentric_local_hour_angle,
                                                                  ('sunLocalHourAngle', 'sunParallaxRightAscension'),
                                                                  self)
        self._internalState['moonElevationAngleWithout'] = Variable(_topocentric_elevation_angle_without,
                                                                    ('geoLatitude', 'moonTopocentricDeclination',
                                                                     'moonTopocentricLocalHour'), self)
        self._internalState['sunElevationAngleWithout'] = Variable(_topocentric_elevation_angle_without,
                                                                   ('geoLatitude', 'sunTopocentricDeclination',
                                                                    'sunTopocentricHourAngle'), self)
        self._internalState['moonAtmosphericRefraction'] = Variable(_atmospheric_refraction_correction,
                                                                    ('pressure', 'temperature',
                                                                     'moonElevationAngleWithout'), self)
        self._internalState['sunAtmosphericRefraction'] = Variable(_atmospheric_refraction_correction,
                                                                   ('pressure', 'temperature',
                                                                    'sunElevationAngleWithout'), self)
        self._internalState['moonElevationAngle'] = Variable(_topocentric_elevation_angle,
                                                             ('moonElevationAngleWithout', 'moonAtmosphericRefraction'),
                                                             self)
        self._internalState['sunElevationAngle'] = Variable(_topocentric_elevation_angle,
                                                            ('sunElevationAngleWithout', 'sunAtmosphericRefraction'),
                                                            self)
        self._internalState['moonZenithAngle'] = Variable(_topocentric_zenith_angle, ('moonElevationAngle',), self)
        self._internalState['sunZenithAngle'] = Variable(_topocentric_zenith_angle, ('sunElevationAngle',), self)
        self._internalState['moonAstronomersAzimuth'] = Variable(_topocentric_astronomers_azimuth_angle,
                                                                 ('moonTopocentricHourAngle', 'geoLatitude',
                                                                  'moonTopocentricDeclination'), self)
        self._internalState['sunAstronomersAzimuth'] = Variable(_topocentric_astronomers_azimuth_angle,
                                                                ('sunTopocentricHourAngle', 'geoLatitude',
                                                                 'sunTopocentricDeclination'), self)
        self._internalState['moonAzimuthAngle'] = Variable(_topocentric_azimuth_angle, ('moonAstronomersAzimuth',),
                                                           self)
        self._internalState['sunAzimuthAngle'] = Variable(_topocentric_azimuth_angle, ('sunAstronomersAzimuth',),
                                                          self)
        self._internalState['sunIncidenceAngle'] = Variable(_spa_incidence_angle,
                                                            ('sunZenithAngle', 'sunAstronomersAzimuth', 'slopeSurface',
                                                             'surfaceAzimuth'), self)


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
