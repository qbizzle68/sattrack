from sattrack.core.exceptions import SattrackException


class RootCountChange(SattrackException):
    """Raised to signal when the quantity of roots found in the zero-function changes."""
    pass


class PassedOrbitPathRange(SattrackException):
    """Raised while the time, during a search for a satellite pass, moves outside a valid range."""
    pass


class SatelliteAlwaysAbove(SattrackException):
    """Raised to signal that no satellite rise or set times exist due to the satellite position always being above
    the horizon. """
    pass


class NoPassException(SattrackException):
    """Raised when an overhead satellite pass does not occur at a given time and location."""
    pass


class NoTLEFound(SattrackException):
    """Raised when no TLE is found on Celestrak."""
    pass
