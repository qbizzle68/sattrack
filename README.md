# Sattrack 0.4.1

Sattrack is a Python package that simplifies computations of celestial mechanics. While Sattrack
supports user defined satellites and celestial bodies, it was designed for real satellites in Earth orbit,
specifically those with a recently generated [two-line element set][tle wiki]. Sattrack makes it easy to
produce an array of data, from SGP4 state vectors, to orbital elements and anomalies at a given time. You
can generate the exact data to solve your specific problem, or use Sattrack's `PassFinder` class to find
the time and details of the next overhead pass of the International Space Station, as well as a list of
all overhead passes for the next several days.

Sattrack uses [Pyevspace], a fast and lightweight Euclidean vector space package for representing vector
quantities and rotating them between reference frames. We implement our own time and geo-position types
with simple interfaces for ease of use and accurate results. The `JulianDate` class allows quick and easy
computations of time differences, while abstracting away the complexities of using Gregorian calendar dates,
as well as timezone differences. There is also support for use with Python's built-in `datetime` module. The
`GeoPosition` class enables painless rotation of vector quantities to and from a topocentric reference frame.

Sattrack Web App
----------------

An online version of Sattrack can be found [here][sattrack web app]. It is a Python Flask app that wraps
around Sattrack, and allows you to use some of Sattrack's tools wherever you have access to the internet.
This is handy for those who want the convenience of not needing to use the package yourselves, or as an
example of the possibilities the package can do for you. 

Install
-------

You can install Sattrack from [PyPi] on the command line with:
``` python
pip install sattrack
```

The source code can be found on [github], and can be cloned using:
``` bash
git clone https://github.com/qbizzle68/sattrack.git
```
The SGP4 module in sattrack is an extension module that uses the Python C-API to wrap around an SGP4
library written in C. This means if you're using an operating system that doesn't have an already built
wheel on PyPi, or you make changes to the source code and (re-)install the package, you'll need a C
compiler installed on your machine (most people won't need to worry about this).

Usage
-----

### Time

Timing is handled using the Julian date. According to the [wikipedia page on the Julian day][JD wiki]:

> The Julian date (JD) of any instant is the Julian day number plus the fraction of a
> day since the preceding noon in Universal Time.

where

> The Julian day is the continuous count of days since the beginning of the Julian period.

The lower level details of the Julian date are out of the scope of these docs, but the key takeaways are:
- Each Julian date has a Julian day number, an integer, that increments every day at 12:00:00 UT.
- The current Julian period started counting on November 24, 4714 BC, in the [Gregorian calendar][gc wiki].
- The Julian date is the Julian day number plus the fraction of the day following 12:00:00 UT.
- The Julian date is not relative to each time zone, but to [Universal Time][ut wiki], however time zones
  are still significant when converting between Julian dates and Gregorian calendar dates.
- January 1, 2000 12:00:00 UT has a Julian date of 2451545.0, and is known as the J2000 epoch.

The `JulianDate` class handles all the complexities of the Julian date computations for you. It will convert
between Gregorian calendar elements and Julian date values, compute future/past dates, printing various
formats, all with varying time zones. The `JulianDate` class also supports interfacing with Python's `datetime`
module via the `JulianDate.fromDatetime()` class method, and the `toDateTime()` instance method.

#### Examples

A `JulianDate` object can be instantiated in various ways:
``` python
from sattrack.api import JulianDate, now
from datetime import datetime, timezone, timedelta

# Assuming this was in a Python script that was ran on February 17, 2024 at 21:35:00 -6 UTC, these would
# all produce approximately the same Julian dates, give or take a few milliseconds.

# Explicit construction of Julian dates:
jdFromComponents = JulianDate(2024, 2, 17, 21, 35, 0, -6)
jdFromValue = JulianDate.fromNumber(2460358.642361111, -6)

# Via Python's datetime module:
offset = timedelta(hours=-6)
tz = timezone(offset)
datetimeInstance = datetime(2024, 2, 17, 21, 35, 0, 0, tzinfo=tz)
jdFromDatetime = JulianDate.fromDatetime(datetimeInstance)

# Get current time and timezone (implicitly via time.localtime().tm_gmtoff, or explicitly):
currentJDImplicit = now()
currentJDExplicit = now(-6)
```

The complexities in computing future/past dates as well as differences in times are completely abstracted away:
``` python
from sattrack.api import JulianDate

# Differences in time measured in solar days.
oneHour = 1 / 24
oneMinute = 1 / 1440
oneSecond = 1/86400

currentTime = JulianDate(2024, 2, 17, 21, 35, 0, -6)
dt = 5 * oneHour + 4 * oneMinute + 3 * oneSecond
futureTime = currentTime.future(dt)
print(futureTime)
# prints '2460358.860451 --- 2024/02/18 02:39:03.0 -6 UTC'

# Negative values move into past.
pastTime = currentTime.future(-dt)
print(pastTime)
# prints '2460358.43816 --- 2024/02/17 16:30:57.0 -6 UTC'

timeSincePast = futureTime - pastTime
print(timeSincePast)
# prints '0.4222916666666667'
print(timeSincePast == (2 * dt))
# prints 'True'

# Order here matters, the sign of the time difference is measured relative to the second parameter.
timeSinceFuture = pastTime - futureTime
print(timeSinceFuture)
# prints '-0.4222916666666667'
```

The `JulianDate` class also supports all comparison operators:
``` python
from sattrack.api import JulianDate

firstTime =  JulianDate(2024, 2, 17, 21, 35, 0, -6)
secondTime = JulianDate(2024, 2, 17, 22, 35, 0, -6)
thirdTime =  JulianDate(2024, 2, 17, 23, 35, 0, -6)

# These all return True:
firstTime < secondTIme
thirdTime >= firstTime
secondTime != thirdTime

# These all return False:
firstTime > secondTime
thirdTime < firstTime
firstTime == secondTIme

# Time zones are just different representations of the same moment in time:
timeInUSA =     JulianDate(2024, 2, 17, 12, 0, 0, -6)
timeInGermany = JulianDate(2024, 2, 17, 19, 0, 0, +1)

# Notice the Julian date number is the same.
print(timeInUSA)
# prints '2460358.25 --- 2024/02/17 12:00:00.0 -6 UTC'
print(timeInGermany)
# prints '2460358.25 --- 2024/02/17 19:00:00.0 +1 UTC'

# Returns True:
timeInUSA == timeInGermany
```

Since times usually represent unique events, it makes sense they may be usable as keys values, meaning
they should be hashable:

``` python
from sattrack.api import JulianDate

firstTime =  JulianDate(2024, 2, 17, 21, 35, 0, -6)
secondTime = JulianDate(2024, 2, 17, 22, 35, 0, -6)
thirdTime =   JulianDate(2024, 2, 17, 23, 35, 0, -6)

# Since JulianDate objects are hashable, they can be used as keys in dicts, or in sets.
events = {firstTime: 'first event', secondTime: 'seconds event', thirdTime: 'third event'}
earliestTime = min(events)

print('The earliest event is', events[earliestTime])
# prints 'The earliest event is first event'
```

### Coordinates

There are several coordinate types in Sattrack, most inheriting from the abstract base class `Coordinates`,
the most useful of which is probably the `GeoPosition` class. This class represents the geo-position of a user, 
and is used for computing [topocentric][topos wiki] values. It takes a latitude and longitude angle as constructor
arguments, and also supports an elevation value, which can be used while computing a geo-position's radius or 
geocentric position vector.

Almost all angle values in Sattrack are in radians. This is because, as a computationally intensive and data based
package, radians are the most useful to use. This is a rule that is kept as consistent across the package as
possible for ease of use, however geo-positions are almost exclusively described using degrees, therefore the only
consistent part of the coordinates objects is how they break this rule. Internally all `Coordinates` objects use
radians, but abstract that away from the user, who will instantiate them and see them printed with degree units.

``` python
from sattrack.api import GeoPosition

# Arguments are latitude then longitude, think of the colloquialism lat/long.
detroit = GeoPosition(42.331429, -83.045753)
# Elevation is in kilometers, Detroit's is 200 meters.
detroitWithElevation = GeoPosition(42.331429, -83.045753, 0.200)
```

The `GeoPosition` class also supports several methods, specifically for computing position data with regard to a
geo-position. For example, since the Earth is rotating relative to the celestial reference frame, the vector
pointing from the center of Earth to a geo-position is dependent on time.

``` python
from sattrack.api import GeoPosition, JulianDate, SIDEREAL_PER_SOLAR

detroit = GeoPosition(42.331429, -83.045753)
time = JulianDate(2024, 2, 17, 21, 35, 0, -6)

positionVector = detroit.getPositionVector(time)
print(positionVector)
# prints '[-2239.8, 4157.31, 4272.89]'

# Time exactly one earth rotation later (one sidereal day).
futureTime = time.future(SIDEREAL_PER_SOLAR)
print(futureTime.date())
# prints '2024/02/18 21:31:04.091 -6 UTC'
futurePositionVector = detroit.getPositionVector(futureTime)
print(futurePositionVector)
# prints '[-2239.8, 4157.31, 4272.89]'
```

### Orbitables

In Sattrack, an orbitable object is an object that can orbit a celestial body. The `Orbitable` class is an
abstract base class that provides an interface for different types of orbiting objects. In Sattrack terms, a
satellite that is described using a set of orbital elements, via the `Elements` class, is implemented with
the `Orbit` class. An orbitable that is described by a TLE, via the `TwoLineElement` class, is implemented
with the `Satellite` class. This may create some confusion, as outside Sattrack, an `Orbit` object would
also be called a satellite, however this is how they are differentiated within Sattrack.

This distinction should be clear enough to avoid any ambiguities in this documentation. To restate it simply:
- An `Orbit` is a satellite using only orbital elements, instantiated with an `Elements` object.
- A `Satellite` is a satellite using a TLE, instantiated with a `TwoLineElement` object.

Both classes implement most of the same methods (via the `Orbitable` interface) so that other modules/classes
in the package do not need to distinguish between the two types.

#### Orbit

The simplest and most customizable satellite in Sattrack is an `Orbit` object. This is why the `Orbit` class was
created, to allow us to explicitly create specific satellite conditions whenever we want. To create a satellite
from a set of orbital elements, first create an `Elements` object, then instantiate an `Orbit` from the elements.
``` python
from sattrack.api import Elements, Orbit, now

# Create a JulianDate object from the current time
time = now()

# The regular constructor takes it's arguments in radians, not degrees
elements = Elements.fromDegrees(90,     # right-ascension of the ascending node
                                51.6,   # inclination
                                0,      # argument of periapsis
                                0.001,  # eccentricity
                                6700,   # semi-major axis in kilometers
                                180,    # mean anomaly, currently at apoapsis
                                time)   # time satellite is at mean anomaly
                                
orbit = Orbit(elements, 'custom-satellite')
```

#### Satellite

The more complex satellite type is a `Satellite` object. The complexity comes from the SGP4 library needed to
propagate the satellite, however this distinction is hidden away by the `TwoLineElement` class and Sattrack's
`sgp4` extension module. Thanks to the `Orbitable` interface, a `Satellite` object can be used, for the most
part, just like an `Orbit` object.

``` python
from sattrack.api import TwoLineElement, Satellite

tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   24048.87310782  .00022368  00000+0  39798-3 0  9993
2 25544  51.6400 191.6020 0001723 274.9313 221.5487 15.50057742439897''')
satellite = Satellite(tle)
```

#### Accuracy
It should be quickly noted here that no real satellite in orbit maintains constant orbital elements due
to the myriad of perturbing forces including but not limited to atmospheric drag, oblateness of the parent
body, solar radiation pressure and gravitational effects from the sun, moon, and other planetary bodies.
This means using orbital elements over a non-zero period of time with the `Orbit` class, should really
only be used when an idealized case is sufficient, such as demonstrating orbital principles or when a
patched conic approximation is needed such as while playing Kerbal Space Program.

A two-line element set is one of the most accurate ways to approximate a satellite's position using the
SGP4 algorithm. With that being said, it should also be noted that the use of two-line element sets and
the SGP4 algorithm is, by definition, an approximation. Therefore, without direct measurement of a satellite,
there is no way to know its exact position or velocity. A TLE will stay with a few kilometers of accuracy
within several days, so you should use a TLE generated close to the time of the data you wish to compute.
An error of several kilometers may sound large, however low earth satellites usually have a radius over
6,700 kilometers, so the error is relatively small and reasonable for most cases.

### Generating Data

At this point whether you have an `Orbit` or `Satellite` object is irrelevant for the most part outside of
accuracy, as all operations shown here can be done with either class.

``` python
from sattrack.api import GeoPosition, JulianDate, TwoLineElement, Satellite, toTopocentric, toTopocentricOffset

geo = GeoPosition(42.331429, -83.045753)
time = JulianDate(2024, 2, 17, 21, 35, 0, -6)

# Most current TLE for the International Space Station for the time used here.
tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   24048.87310782  .00022368  00000+0  39798-3 0  9993
2 25544  51.6400 191.6020 0001723 274.9313 221.5487 15.50057742439897''')
satellite = Satellite(tle)

position, velocity = satellite.getState(time)
print(position, velocity)
# prints '[2825.37, 4178.21, -4557.74] [-6.87903, 1.25872, -3.11139]'

# Transform state vectors to the SEZ reference frame.
topocentricPosition = toTopocentricOffset(position, geo, time)
topocentricVelocity = toTopocentric(velocity, geo, time)
print(topocentricPosition, topocentricVelocity)
# prints '[4922.7, -4469.09, -7709.13] [5.24356, 5.45901, 1.13596]'
```

### Overhead Satellite Passes

#### Overview

One of the most intriguing parts of celestial mechanics that even a lay person can enjoy is when a satellite
visibly passes overhead. Sattrack shields the complexities for quickly and accurately finding the next time
a satellite will pass overhead with the `PassFinder` class. The `PassFinder` methods `computeNextPass()` and
`computePassList()` methods return either a single or `list` of `SatellitePass` objects respectively.

#### Satellite Pass

The `SatellitePass` is container of `PositionInfo` instances. Each `PositionInfo` object contains information
about a specific event that occurs during a pass such as rise time, maximum altitude time and set time (all
passes have these three events). A printed `SatellitePass` looks like the following:

```
                 Pass details for ISS (ZARYA), at 2024/02/18 04:30:59.791 -6 UTC
     instance      |     time     | altitude |   azimuth    | illuminated | unobscured | visible
-------------------------------------------------------------------------------------------------
       rise        | 04:25:48.838 |   0.00   | 211.73 (SSW) |    False    |    True    |  False
-------------------------------------------------------------------------------------------------
 first illuminated | 04:30:26.189 |  29.34   | 159.59 (SSE) |    True     |    True    |  True
-------------------------------------------------------------------------------------------------
        max        | 04:30:59.791 |  31.64   | 138.45 (SE)  |    True     |    True    |  True
-------------------------------------------------------------------------------------------------
        set        | 04:36:16.96  |   0.00   | 62.78  (ENE) |    True     |    True    |  True
```

Each column corresponds to a different event described by the instance value, and has a time, altitude angle,
and azimuth angle. The azimuth angle refers to a compass heading to look towards at ground level, and the
altitude angle refers to how high above that point the satellite is at that time (90 degrees is straight up).
Each event also has visibility information: illuminated, unobscured and visible. Illuminated refers to whether
the satellite is **_illuminated_** by sunlight at that time. Satellites do not give off their own light, and
can only be seen when they reflect sunlight (just like the moon). Unobscured refers to the fact that the Sun
is too bright for us to see any other celestial objects in the sky (except for the bright Moon). Therefore,
celestial objects can almost always only be seen after the Sun has set and before it rises, when they are
**_unobscured_** by the Sun. In Sattrack terms, a satellite is visible if and only if it is illuminated
**_and_** unobscured.

``` python
# other code used to produce our example SatellitePass above

riseInfo = satellitePass.riseInfo
print(riseInfo.illuminated)
# prints False
print(riseInfo.unobscured)
# prints True
print(riseInfo.visible)
# prints False
```

There are also visibility values for an entire `SatellitePass`. If any event of a pass is illuminated, the
illuminated attribute of the pass is True, otherwise it is False. The same logic follows for the unobscured
and visible attributes.

#### Finding Passes

Passes can be easily found using the `PassFinder` class, which required an `Orbitable` and `GeoPosition` to
instantiate, and a list of `SatellitePass` objects can be computed with the `computePassList()` method.

``` python 
tle = TwoLineElement('''ISS (ZARYA)             
1 25544U 98067A   24048.87310782  .00022368  00000+0  39798-3 0  9993
2 25544  51.6400 191.6020 0001723 274.9313 221.5487 15.50057742439897''')
satellite = Satellite(tle)
geo = GeoPosition(42.331429, -83.045753)
time = JulianDate(2024, 2, 17, 21, 35, 0, -6)

finder = PassFinder(satellite, geo)
# Compute passes for the next seven days.
passList = finder.computePassList(time, 7)
```

You can either print and visually inspect each pass, or refine the list however you wish.

``` python 
for satellitePass in passList:
  print(satellitePass)
# prints a long list of SatellitePass objects
  
visibleList = [satellitePass for satellitePass in passList if satellitePass.visible]
print(len(visibleLlist))
# prints 17

# Compute a list of passes that are visible and are higher than 70 degrees.
highPasses = [satellitePass for satellitePass in visibleList if satellitePass.maxInfo.altitude > 70]
print(len(highPasses))
# prints 2
```

#### Pass Events

Each satellite pass can have several events besides the default rise, set and maximum. These events are
determined by where in the pass a satellite (if at all) becomes illuminated, unobscured, or visible. A
satellite that enters or exits the Earth's shadow during a pass will have a 'first illuminated' or a
'last illuminated' event respectively. A satellite pass that occurs during sunrise or sunset will have
a 'last unobscured' or a 'first unobscured' event respectively. Events like these will only be displayed
when the illuminated or unobscured value changes at a single point during the pass. For example a satellite
pass that takes place in the middle of the day will have no first or last unobscured event, as it's unobscured
for the entire pass.

There are several other events that may provide utility in searching for specific pass types. If the event is
applicable, its value will be a `PositionInfo` instance, otherwise it will be `None`. Some events may be
the same `PositionInfo` instance. For example, in a pass that is visible during its entirety, the `maxInfo`,
`maxVisibleInfo`, `maxIlluminatedInfo`, and `maxUnobscuredInfo` attributes will be the same `PositionInfo`
instance. Similarly, all first info's will be the same as `riseInfo`, and all last info's will be the same
as `setInfo`. All possible event attribute names are: `riseInfo`, `setInfo`, `maxInfo`, `firstIlluminatedInfo`,
`lastIlluminatedInfo`, `maxIlluminatedInfo`, `firstUnobscuredInfo`, `lastUnobscuredInfo`, `maxUnobscuredInfo`,
`firstVisibleInfo`, `lastVisibleInfo`, and `maxVisibleInfo`.

License
-------
[MIT][mit license]

[tle wiki]: https://en.wikipedia.org/wiki/Two-line_element_set
[Pyevspace]: https://pypi.org/project/pyevspace/
[sattrack web app]: https://sattrack.qbizzle.com/
[PyPi]: https://pypi.org/
[github]: https://github.com/qbizzle68/sattrack/
[JD wiki]: https://en.wikipedia.org/wiki/Julian_day
[ut wiki]: https://en.wikipedia.org/wiki/Universal_Time
[gc wiki]: https://en.wikipedia.org/wiki/Gregorian_calendar
[topos wiki]: https://en.wikipedia.org/wiki/Horizontal_coordinate_system
[mit license]: https://choosealicense.com/licenses/mit/
