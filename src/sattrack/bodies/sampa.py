from math import radians

from sattrack.core.juliandate import JulianDate
from sattrack.util.constants import DELTAT

from sattrack.bodies._sampa import _mpaMoonMeanLongitude, _mpaMoonMeanElongation, _mpaSunMeanAnomaly, \
    _mpaMoonMeanAnomaly, _mpaMoonArgumentLatitude, _mpaETerm, _mpaLRTable, _mpaLTerm, _mpaRTerm, \
    _mpaBTerm, _mpaA1Term, _mpaA2Term, _mpaA3Term, _mpaDeltaL, _mpaDeltaB, _mpaMoonLongitude, \
    _mpaMoonLatitude, _mpaMoonDistance, _mpaMoonParallax, _xValues, _xyTable, _nutationLongitude, \
    _nutationObliquity, _meanObliquity, _trueObliquity, _mpaApparentMoonLongitude, _meanSiderealTime, \
    _apparentSiderealTime, _right_ascension, _declination, _spaEarthHeliocentricLongitude, \
    _spaEarthHeliocentricLatitude, _spaGeocentricLongitude, _spaGeocentricLatitude, _spaAberrationCorrection, \
    _spaApparentSunLongitude, _equationOfTime, _spaEarthHeliocentricRadius, _localHourAngle, _uTerm, \
    _xTerm, _yTerm, _parallaxRightAscension, _topocentricRightAscension, _topocentricDeclination, \
    _spaEquitorialParallaxSun, _topocentricLocalHourAngle, _topocentricElevationAngleWithout, \
    _atmosphericRefractionCorrection, _topocentricElevationAngle, _topocentricZenithAngle, \
    _topocentricAstronomersAzimuthAngle, _topocentricAzimuthAngle, _spaIncidenceAngle


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
            'moonMeanLongitude': Variable(_mpaMoonMeanLongitude, ('JCE',), self),
            'moonMeanElongation': Variable(_mpaMoonMeanElongation, ('JCE',), self),
            'sunMeanAnomaly': Variable(_mpaSunMeanAnomaly, ('JCE',), self),
            'moonMeanAnomaly': Variable(_mpaMoonMeanAnomaly, ('JCE',), self),
            'moonArgumentLatitude': Variable(_mpaMoonArgumentLatitude, ('JCE',), self),
            'ETerm': Variable(_mpaETerm, ('JCE',), self),
            'lrProductTable': Variable(_mpaLRTable, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
                                                       'moonArgumentLatitude', 'ETerm'), self),
            'lTerm': Variable(_mpaLTerm, ('lrProductTable',), self),
            'rTerm': Variable(_mpaRTerm, ('lrProductTable',), self),
            'bTerm': Variable(_mpaBTerm, ('moonMeanElongation', 'sunMeanAnomaly', 'moonMeanAnomaly',
                                            'moonArgumentLatitude', 'ETerm'), self),
            'a1Term': Variable(_mpaA1Term, ('JCE',), self),
            'a2Term': Variable(_mpaA2Term, ('JCE',), self),
            'a3Term': Variable(_mpaA3Term, ('JCE',), self),
            'deltal': Variable(_mpaDeltaL, ('a1Term', 'a2Term', 'moonMeanLongitude', 'moonArgumentLatitude'), self),
            'deltab': Variable(_mpaDeltaB, ('a1Term', 'a3Term', 'moonMeanLongitude', 'moonMeanAnomaly',
                                              'moonArgumentLatitude'), self),
            'moonLongitude': Variable(_mpaMoonLongitude, ('moonMeanLongitude', 'lTerm', 'deltal'), self),
            'moonLatitude': Variable(_mpaMoonLatitude, ('bTerm', 'deltab'), self),
            'moonDistance': Variable(_mpaMoonDistance, ('rTerm',), self),
            'moonParallax': Variable(_mpaMoonParallax, ('moonDistance',), self),
            'xValues': Variable(_xValues, ('JCE',), self),
            'xyProductTable': Variable(_xyTable, ('xValues',), self),
            'nutationLongitude': Variable(_nutationLongitude, ('JCE', 'xyProductTable'), self),
            'nutationObliquity': Variable(_nutationObliquity, ('JCE', 'xyProductTable'), self),
            'meanObliquity': Variable(_meanObliquity, ('JME',), self),
            'trueObliquity': Variable(_trueObliquity, ('meanObliquity', 'nutationObliquity'), self),
            'apparentMoonLongitude': Variable(_mpaApparentMoonLongitude, ('moonLongitude', 'nutationLongitude'),
                                              self),
            'meanSiderealTime': Variable(_meanSiderealTime, ('JD', 'JC'), self),
            'apparentSiderealTime': Variable(_apparentSiderealTime, ('meanSiderealTime', 'nutationLongitude',
                                                                       'nutationObliquity'), self),
            'moonRightAscension': Variable(_right_ascension, ('moonLongitude', 'moonLatitude',
                                                              'trueObliquity'), self),
            'moonDeclination': Variable(_declination, ('moonLongitude', 'moonLatitude', 'trueObliquity'),
                                        self),
            'earthHeliocentricLongitude': Variable(_spaEarthHeliocentricLongitude, ('JME',), self),
            'earthHeliocentricLatitude': Variable(_spaEarthHeliocentricLatitude, ('JME',), self),
            'sunDistance': Variable(_spaEarthHeliocentricRadius, ('JME',), self),
            'sunLongitude': Variable(_spaGeocentricLongitude, ('earthHeliocentricLongitude',), self),
            'sunLatitude': Variable(_spaGeocentricLatitude, ('earthHeliocentricLatitude',), self),
            'sunAberrationCorrection': Variable(_spaAberrationCorrection, ('sunDistance',), self),
            'apparentSunLongitude': Variable(_spaApparentSunLongitude, ('sunLongitude', 'nutationLongitude',
                                                                           'sunAberrationCorrection'), self),
            'sunRightAscension': Variable(_right_ascension, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
            'sunDeclination': Variable(_declination, ('sunLongitude', 'sunLatitude', 'trueObliquity'), self),
            'equationOfTime': Variable(_equationOfTime,
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
        self._internalState['moonLocalHourAngle'] = Variable(_localHourAngle, ('apparentSiderealTime', 'geoLongitude',
                                                                                 'moonRightAscension'), self)
        self._internalState['uTerm'] = Variable(_uTerm, ('geoLatitude',), self)
        self._internalState['xTerm'] = Variable(_xTerm, ('uTerm', 'geoLatitude', 'elevation'), self)
        self._internalState['yTerm'] = Variable(_yTerm, ('uTerm', 'geoLatitude', 'elevation'), self)
        self._internalState['moonParallaxRightAscension'] = Variable(_parallaxRightAscension,
                                                                     ('xTerm', 'moonParallax', 'moonLocalHourAngle',
                                                                      'moonDeclination'), self)
        self._internalState['moonTopocentricRightAscension'] = Variable(_topocentricRightAscension,
                                                                        ('moonRightAscension',
                                                                         'moonParallaxRightAscension'), self)
        self._internalState['moonTopocentricDeclination'] = Variable(_topocentricDeclination,
                                                                     ('yTerm', 'moonDeclination', 'moonParallax',
                                                                      'moonParallaxRightAscension',
                                                                      'moonLocalHourAngle'), self)
        self._internalState['sunLocalHourAngle'] = Variable(_localHourAngle, ('apparentSiderealTime', 'geoLongitude',
                                                                                'sunRightAscension'), self)
        self._internalState['sunParallax'] = Variable(_spaEquitorialParallaxSun, ('sunDistance',), self)
        self._internalState['sunParallaxRightAscension'] = Variable(_parallaxRightAscension,
                                                                    ('xTerm', 'sunParallax', 'sunLocalHourAngle',
                                                                     'sunDeclination'), self)
        self._internalState['sunTopocentricRightAscension'] = Variable(_topocentricRightAscension,
                                                                       ('sunRightAscension',
                                                                        'sunParallaxRightAscension'), self)
        self._internalState['sunTopocentricDeclination'] = Variable(_topocentricDeclination,
                                                                    ('yTerm', 'sunDeclination', 'sunParallax',
                                                                     'sunParallaxRightAscension', 'sunLocalHourAngle'),
                                                                    self)
        self._internalState['moonTopocentricHourAngle'] = Variable(_topocentricLocalHourAngle,
                                                                   ('moonLocalHourAngle', 'moonParallaxRightAscension'),
                                                                   self)
        self._internalState['sunTopocentricHourAngle'] = Variable(_topocentricLocalHourAngle,
                                                                  ('sunLocalHourAngle', 'sunParallaxRightAscension'),
                                                                  self)
        self._internalState['moonElevationAngleWithout'] = Variable(_topocentricElevationAngleWithout,
                                                                    ('geoLatitude', 'moonTopocentricDeclination',
                                                                     'moonTopocentricLocalHour'), self)
        self._internalState['sunElevationAngleWithout'] = Variable(_topocentricElevationAngleWithout,
                                                                   ('geoLatitude', 'sunTopocentricDeclination',
                                                                    'sunTopocentricHourAngle'), self)
        self._internalState['moonAtmosphericRefraction'] = Variable(_atmosphericRefractionCorrection,
                                                                    ('pressure', 'temperature',
                                                                     'moonElevationAngleWithout'), self)
        self._internalState['sunAtmosphericRefraction'] = Variable(_atmosphericRefractionCorrection,
                                                                   ('pressure', 'temperature',
                                                                    'sunElevationAngleWithout'), self)
        self._internalState['moonElevationAngle'] = Variable(_topocentricElevationAngle,
                                                             ('moonElevationAngleWithout', 'moonAtmosphericRefraction'),
                                                             self)
        self._internalState['sunElevationAngle'] = Variable(_topocentricElevationAngle,
                                                            ('sunElevationAngleWithout', 'sunAtmosphericRefraction'),
                                                            self)
        self._internalState['moonZenithAngle'] = Variable(_topocentricZenithAngle, ('moonElevationAngle',), self)
        self._internalState['sunZenithAngle'] = Variable(_topocentricZenithAngle, ('sunElevationAngle',), self)
        self._internalState['moonAstronomersAzimuth'] = Variable(_topocentricAstronomersAzimuthAngle,
                                                                 ('moonTopocentricHourAngle', 'geoLatitude',
                                                                  'moonTopocentricDeclination'), self)
        self._internalState['sunAstronomersAzimuth'] = Variable(_topocentricAstronomersAzimuthAngle,
                                                                ('sunTopocentricHourAngle', 'geoLatitude',
                                                                 'sunTopocentricDeclination'), self)
        self._internalState['moonAzimuthAngle'] = Variable(_topocentricAzimuthAngle, ('moonAstronomersAzimuth',),
                                                           self)
        self._internalState['sunAzimuthAngle'] = Variable(_topocentricAzimuthAngle, ('sunAstronomersAzimuth',),
                                                          self)
        self._internalState['sunIncidenceAngle'] = Variable(_spaIncidenceAngle,
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
