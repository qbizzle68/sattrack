# Coordinate module documentation
The `coordinates` module supplies classes for representing coordinates for Earth and celestial reference frames.

## Types
*class* `coordinates.Coordinates`
- An abstract base class for creating coordinate objects.

*class* `coordinates.GeoPosition`
- A class to hold geo-position information.

*class* `coordinates.CelestialCoordinates`
- A class to hold position coordinates in the celestial coordinate system.

## Methods
The `coordinates` module imports three methods for getting Earth values at certain latitudes. Due to the oblate spheroid
shape of the earth, there are multiple types of latitudes namely 
[geodetic and geocentric](https://en.wikipedia.org/wiki/Latitude#Geodetic_and_geocentric_latitudes). There are three 
methods to computes values dependent on a form of latitude. A fourth module method computes the sub-point of a satellite.

*method* `coordinates.`**geocentricToGeodetic(geocentricLatitude)**
- Converts a geocentric latitude into a geodetic latitude. Both the geocentricLatitude parameter and return value are in
degrees.

*method* `coordinates.`**geodeticToGeocentric(geodeticLatitude)**
- Converts a geodetic latitude into a geocentric latitude. Both the geodeticLatitude parameter and return value are in
degrees.

*method* `coordinates.`**radiusAtLatitude(latitude)**
- Computes the radius of Earth at a given geodetic latitude, which varies between the equator and the poles. latitude is
the geodetic latitude in degrees and the return value is in kilometers.

*method* `coordinates.`**getSubPoint(position, jd)**
- Return the geo-position directly underneath a satellites position at a given time.
- Only the direction of the position vector is significant, not its magnitude.
- The elevation of the returned GeoPosition object is always 0.

## Coordinates objects
An abstract base class for deriving coordinate like objects for the module. 

### Constructor:
*class* `coordinates.`**Coordinates(latitude, longitude)**
- Latitude must be a numeric value between -90 and 90 degrees, otherwise a TypeError or ValueError is raised.
- Longitude must be a numeric value otherwise a TypeError is raised. 
- The longitude value is modulated into its valid ranges automatically.

### Instance attributes:
`Coordinates.`**latitude**
- Returns and sets the latitude of the coordinate.

`Coordinates.`**longitude**
- Returns and sets the longitude of the coordinate.

### Instance Methods:

`Coordinates.`**\_\_repr__()**
- Returns a JSON string representation of the coordinate. The \_\_iter__() method is not implemented, so for this
method to work for subclasses, the \_\_iter__() method must be overwritten.
  
## GeoPosition objects
- The `GeoPosition` class is a subclass of the `Coordinates` class, whose valid ranges of values are -90 to 90 degrees 
for latitude and -180 to 180 for longitude. Longitude values will automatically be modulated to values in their valid
range.
- Due to the oblate spheroid shape of earth, each point's latitude angle can be described as geodetic or geocentric.
Unless otherwise stated, when describing latitude, geodetic latitude is to be assumed.
- This class also contains an elevation value, which is used to derive more precise position vectors for a position if
the elevation is known.

### Constructor:
*class* `coordinates.`**GeoPosition(latitude, longitude, elevation)**
- All parameters must be numeric types otherwise raises a TypeError.
- latitude must be in the range of -90 to 90 otherwise a ValueError is raised, and longitude is modulated into a range 
of -180 to 180 if outside it.
- elevation is in kilometers.

### Instance attributes:
`GeoPosition.`**latitude**
- Returns and sets the latitude of the coordinate. When setting, the value must be a numeric type and in the range of
-90 to 90 degrees.

`GeoPosition.`**longitude**
- Returns and sets the longitude of the coordinate. When setting, the value must be a numeric type and will be modulated
if not in the range of -180 to 180 degrees.

`GeoPosition.`**elevation**
- Returns and sets the elevation in kilometers.

### Instance Methods:

`GeoPosition.`**\_\_iter__()**
- Yields a dictionary with the latitude, longitude and elevation values.

`GeoPosition.`**\__\_str__()**
- Returns a simple string of the form 'latitude: {lat}, longitude: {lng}, elevation: {elevation}'.

`GeoPosition.`**\_\_repr__()**
- Returns a JSON string representation of the coordinate.

`GeoPosition.`**\_\reduce__()**
- Provides support for pickle.

`GeoPosition.`**getRadius()**
- Returns the Earth radius at the surface point accounting for the geodetic latitude and the elevation. Return value in
kilometers.

`GeoPosition.`**getGeocentricLatitude()**
- Converts the geodetic latitude into the geocentric latitude in degrees. Guaranteed to be between -90 and 90.

`GeoPosition.`**getPositionVector(jd=None)**
- Computes the position vector of the geo-position.
- If jd is not None, the Earth offset angle at the given time is also accounted for.
- The magnitude of the position vector is equal to the radius at latitude plus the elevation.

`GeoPosition.`**getZenithVector(jd=None)**
- Computes the vector normal to the horizon plane of the position.
- If jd is not None, the Earth offset angle at the given time is also accounted for.
- The magnitude of the returned vector is always 1.0.

## CelestialCoordinates objects
The `CelestialCoordinates` object defines a position in the celestial coordinate system, where longitude angles are 
referred to as right-ascension, and latitude angles are referred to as declination.

### Constructor:
*class* `coordinates.`**CelestialCoordinates(rightAscension, declination)**
- declination must be a numeric type and between the values -90 to 90 degrees otherwise a TypeError or ValueError is
raised respectively.
- rightAscension must be a numeric type otherwise a TypeError is raised. The value will always be modulated between 
0 and 24 hours, where 15 degrees equals 1 hour.
- The parameters are in reverse order of the `Coordinates` constructor, because the idiom of referencing celestial
coordinates is 'right-ascension, declination' as opposed to 'latitude, longitude'.

### Instance attributes:
`CelestialCoordinates.`**latitude**
- Used to access the declination, return value is always between -90 and 90 degrees.
- Setting this value must be a numeric type between -90 and 90 degrees or a TypeError or ValueError is raised.

`CelestialCoordinates.`**longitude**
- Used to access the right-ascension, return value is always between 0 and 24 hours.
- Setting this value must be a numeric type or a TypeError is raised. Value will be modulated between 0 and 24.

### Instance Methods:

`CelestialCoordinates.`**\_\_iter__()**
- Yields a dictionary of the right-ascension and declination.

`CelestialCoordinates.`**\_\_str__()**
- Returns a simple string of the coordinates with the format 'right-ascension: {ra}, declination: {dec}'

`CelestialCoordinates.`**\_\_repr__()**
- Returns a JSON string representation of the coordinate.

`CelestialCoordinates.`**\_\_reduce__()**
- Provides support for pickle.
