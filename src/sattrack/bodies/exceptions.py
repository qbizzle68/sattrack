from sattrack.core.exceptions import SattrackException


class SunRiseSetException(SattrackException):
    """Raised when the sun doesn't rise or set at a given time and location."""
    pass
