from sattrack.core.exceptions import SattrackException


class NoSatelliteEclipseException(SattrackException):
    """Exception raised when a satellite is always in sunlight, and does not enter Earth's shadow."""
    pass


class NoFunctionRootFound(SattrackException):
    """Exception raised when a root isn't found looking for eclipse times."""
    pass
