__all__ = []

# import subpackages
from .spacetime import *
__all__ += spacetime.__all__
from .util import *
__all__ += util.__all__

# import modules
from .coordinates import *
__all__ += coordinates.__all__
from .eclipse import *
__all__ += eclipse.__all__
from .exceptions import *
__all__ += exceptions.__all__
from .moon import *
__all__ += moon.__all__
from .orbit import *
__all__ += orbit.__all__
from .sampa import *
__all__ += sampa.__all__
from .sun import *
__all__ += sun.__all__
from .tle import *
__all__ += tle.__all__
from .topocentric import *
__all__ += topocentric.__all__

# if in debug mode, import internal modules as well
if __debug__ is True:
    from ._coordinates import *
    __all__ += _coordinates.__all__
    from ._orbit import *
    __all__ += _orbit.__all__
    from ._sampa import *
    __all__ += _sampa.__all__
    from ._topocentric import *
    __all__ += _topocentric.__all__

"""Compute positioning and angle data for satellites in earth orbit.

This package can create raw position and velocity data for earth satellites, and interpret the results depending on what
data you wish to find. The most simple are various position angles like right-ascension and declination, hour-angle and
azimuth and altitude. The more complex computations include if a geo-position is 'under' an orbital path or the next
pass above a horizon. Satellites can be created by using raw orbital parameters (raan, inclination, eccentricity etc.)
or by using TLE data published by various online sources, the latter of course being much more accurate for earth based
satellites.

Usage
_____

The easiest way to create a satellite is to create a `Satellite` object using a TLE from the `getTle()` method which
uses [Celestrak](http://www.celestrak.com) via the `requests` module.

>>> from sattrack.tle import getTle
>>> from sattrack.spacetime.juliandate import now
>>> from sattrack.orbit import Satellite
>>>
>>> tle = getTle('zarya') # get the most recent TLE for the ISS
>>> jd = now() # get a JulianDate object set to the current time
>>> iss = Satellite(tle) # create a Satellite object for the ISS
>>>
>>> state = iss.getState(jd) # get a tuple of vectors with the position and velocity vectors in TEME reference frame
>>> print(state[0])
[ -4502.126853, -1763.788424, -4787.006412 ]
>>> print(state[1])
[ 4.875115, -5.263727, -2.646547 ]

An orbit can also be created from raw orbital elements. The fromDegrees class method is used here for simplicity because
the default constructor takes values in radians.

>>> from sattrack.orbit import Elements, Orbit
>>>
>>> elements = Elements.fromDegrees(329.765, 51.6438, 59.0052, 0.000677, 6795.067, 184.78, jd)
>>> iss = Orbit(elements)
>>>
>>> state = iss.getState(jd) # get a tuple of vectors with the position and velocity vectors
>>> print(state[0])
[ -4501.715073, -1757.589387, -4783.376009 ]
>>> print(state[1])
[ 4.875803, -5.270390, -2.651550 ]

The `Orbit` and `Satellite` objects are derived from a common interface, which allows them to be used interchangabley by
other methods in the package. The `Satellite` object is the most accurate when handling satellites in Earth orbit, but
the `Orbit` object allows for handling any orbit around any orbitable body. The package also provides a Body object for
describing the Body around which an Orbit object orbits.

SatellitePass
-------------

The most sophisticated thing the package can do is compute the next time a satellite passes over a given location, as
well as detailed information about the pass. An overhead pass is possible if the satellite is above a given horizon, if
the sun is below the horizon, and the satellite is in sunlight (not in Earth's shadow). This does not account for other
light that can obscure the satellite (i.e. light pollution from highly populated areas).

A pass is considered the time between a satellite rising above the horizon, and setting back below it. There are at
least three major events in every pass: rise, set and maximum altitude. The information for these events are stored in
an object called `PositionInfo` in the `topocentric` module. The `PositionInfo` object holds the time, altitude, azimuth
and other visibility boolean values described next.

If the satellite is in the sunlight it is refered to as illuminated. If the sun is below the horizon while the
satellite is above the horizon, it is refered to as unobscured. A satellite is therefore refered to as visible if it is
illuminated and unobscured. During a pass, a satellite can become illuminated or unobscured, this happens when the
satellite leaves Earth's shadow or sunset occurs during the pass. There are therefore up to four additional events that
can occur during as pass, first-illuminated, last-illuminated, first-unobscured and last-unobscured. Each of these
events are stored in a `PositionInfo` within the `SatellitePass` object. Each `PositionInfo` object has an illumnated,
unobscured and visible boolean property that relates to those states for each event.

There is also an overall illumnated, unobscured and visible state for the `SatellitePass` object. Any of these values
are true if they are true for any event at any point during the pass. Keep in mind that a pass can be illuminated and
unobscured at different times, which therefore means it was never visible, so the pass visibility values will read
illuminated: True, unobscured: True, visible: False. This can occur if the satellite rises illuminated and obscured,
then becomes unilluminated by passing into the Earth's shadow, then becomes unobscured due to sunset. At no point in the
pass was the satellite illuminated and unobscured at the same time, therefore it is never visible. This is a contrived
example but is intended to show how to interpret a `SatellitePass` object.

A further note on passes: if an event is valid at rise or set time, it is not set in the SatellitePass object. For
example if a satellite rises illuminated, there will not be a first-illuminated `PositionInfo` set. Likewise, if a
satellite is unobscured when it rises and sets, there will not be a last-unobscured `PositionInfo` set. Only alternate
events that occur **during** the pass will be considered.
"""