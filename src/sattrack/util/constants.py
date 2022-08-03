from math import pi

# Mass constants
NEWTON_G = 6.67408e-11
EARTH_MU = 3.986004418e14 * 1e-9  # converted to km^3s^-2
SUN_MU = 1.32712440018e20 * 1e-9  # converted to km^3s^-2

#   Earth constants
EARTH_EQUITORIAL_RADIUS = 6378.135
EARTH_FLATTENING = 1.0 / 298.26
EARTH_POLAR_RADIUS = EARTH_EQUITORIAL_RADIUS * (1 - EARTH_FLATTENING)
CJ2 = -2.064734896e14

#   Sun constants
SUN_RADIUS = 6.957e5    # km

#   Distance constants
AU = 1.495978707e8  # km

#   Math constants
TWOPI = 2.0 * pi