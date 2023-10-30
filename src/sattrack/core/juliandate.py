import json
import time
import datetime

from sattrack.util.constants import SECONDS_PER_DAY


def _jdToGregorian(number: int, fraction: float, timezone: float) -> (int, int, int, int, int, float):
    """Convert a Julian date to Gregorian calendar components."""

    # incorporate timezone to offset Julian date value
    extraDay, F = divmod(fraction + (timezone / 24.0) + 0.5, 1.0)
    # algorithm taken from the 'Solar Position Algorithm for Solar Radiation Applications' paper in appendix A.3
    Z = int(number + extraDay)
    if Z < 2299161:
        A = Z
    else:
        B = int((Z - 1867216.25) / 36524.25)
        A = Z + 1 + B - int(B / 4)
    C = A + 1524
    D = int((C - 122.1) / 365.25)
    G = int(365.25 * D)
    I = int((C - G) / 30.6001)
    # do not add the fractional part to the day
    d = C - G - int(30.6001 * I)
    m = I - 1 if I < 14 else I - 13
    y = D - 4716 if m > 2 else D - 4715

    # convert day fraction (measured from 0 hour) to time components
    s = F * 86400.0
    h = int(s / 3600.0)
    s -= h * 3600.0
    mi = int(s / 60.0)
    s -= mi * 60.0

    return y, m, d, h, mi, s


def _dateToJd(year: int, month: int, day: int, hour: int, minute: int, second: float, timezone: float = 0):
    """Logic to compute the Julian day number and fraction from Gregorian date components."""

    if month == 1 or month == 2:
        year = year - 1
        month = month + 12

    # We want to separate integer number and float fraction to maintain precision.
    # This is a modified version of the JD conversion, the 0.5 if moved 'up' to the float part.
    D = day + (hour / 24.0) + (minute / 1440.0) + (second / 86400.0) - (timezone / 24.0) - 0.5
    dayInt, dayFrac = divmod(D, 1.0)

    # We only add in integer part of the day here.
    dayNumber = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + dayInt - 1524

    if dayNumber > 2299160:
        A = int(year / 100)
        B = 2 - A + int(A / 4)
        dayNumber += B

    # We can take every precaution above to ensure dayNumber stays an integer in out computations, but the simplest
    # and most guaranteed way is to just convert it here.
    return int(dayNumber), dayFrac


class JulianDate:
    """Julian Date object which represents a moment in time as an integer day number and a fraction of the day after
    12 noon. An object can be set with Gregorian calendar components via the class constructor, or with a known day
    number and fraction via the fromNumber() class method. See now() for getting a JulianDate object with the current
    time."""

    __slots__ = '_dayNumber', '_dayFraction', '_timezone'

    def __init__(self, year: int, month: int, day: int, hour: int, minute: int, second: float, timezone: int = 0):
        number, fraction = _dateToJd(year, month, day, hour, minute, second, timezone)

        self._dayNumber = number
        self._dayFraction = fraction
        self._timezone = timezone

    @classmethod
    def fromNumber(cls, number: float, timezone: float = 0.0) -> 'JulianDate':
        """Creates a new JulianDate object directly from a Julian day number with optional timezone offset."""
        rtn = object.__new__(cls)

        rtn._dayNumber = int(number)
        rtn._dayFraction = number - rtn._dayNumber
        rtn._timezone = timezone

        return rtn

    @classmethod
    def fromDatetime(cls, date: datetime.datetime) -> 'JulianDate':
        """Creates a new JulianDate object from a Python datetime.datetime instance."""

        # Check if date is aware or naive.
        if date.tzinfo is None:
            tz = 0
        else:
            tz = date.tzinfo.utcoffset(None) / datetime.timedelta(hours=1)

        seconds = date.second + date.microsecond / 1e6
        return JulianDate(date.year, date.month, date.day, date.hour, date.minute, seconds, tz)

    def __str__(self) -> str:
        """Creates a string representation of the JulianDate."""

        return str(round(self.value, 6)) + ' --- ' + self.date(self._timezone)

    def __repr__(self) -> str:
        """Creates a string representation of the JulianDate."""

        y, m, d, h, mi, s = _jdToGregorian(self._dayNumber, self._dayFraction, self._timezone)
        return f'JulianDate({y}, {m}, {d}, {h}, {mi}, {s}, {self._timezone})'

    def toDict(self) -> dict:
        """Returns a dictionary of the JulianDate to create json formats of other types containing a GeoPosition."""

        return {"dayNumber": self._dayNumber, "dayFraction": self._dayFraction}

    def toJson(self) -> str:
        """Returns a string of the coordinate in json format."""

        return json.dumps(self, default=lambda o: o.toDict())

    # operators
    def __add__(self, other: datetime.timedelta) -> 'JulianDate':
        if isinstance(other, datetime.timedelta):
            deltaSolarDays = other / datetime.timedelta(days=1)
            return JulianDate.fromNumber(self.value + deltaSolarDays, self._timezone)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other: 'JulianDate') -> float:
        if isinstance(other, JulianDate):
            return self.value - other.value
        return NotImplemented

    def __eq__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value == other.value
        return NotImplemented

    def __ne__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value != other.value
        return NotImplemented

    def __lt__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value < other.value
        return NotImplemented

    def __le__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value <= other.value
        return NotImplemented

    def __gt__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value > other.value
        return NotImplemented

    def __ge__(self, other: 'JulianDate') -> bool:
        if isinstance(other, JulianDate):
            return self.value >= other.value
        return NotImplemented

    def __hash__(self):
        return hash((self._dayNumber, self._dayFraction, self._timezone))

    def __reduce__(self):
        return self.__class__.fromNumber, (self._dayNumber + self._dayFraction, self._timezone)

    # read-only properties
    @property
    def value(self) -> float:
        return self._dayNumber + self._dayFraction

    @property
    def number(self) -> int:
        return self._dayNumber

    @property
    def fraction(self) -> float:
        return self._dayFraction

    @property
    def timezone(self) -> float:
        return self._timezone

    def future(self, days: int | float) -> 'JulianDate':
        """Create a new JulianDate object in the future or past relative to the calling instance in solar days. A
        positive value moves forward in time, negative is backward."""

        return JulianDate.fromNumber(self._dayNumber + self._dayFraction + days, self._timezone)

    def date(self, timezone: float = None, n: int = 3) -> str:
        """Returns a string with the Julian date represented as a string as mm/dd/year hh:mm:ss +/- tz UTC."""
        if timezone is None:
            timezone = self._timezone

        y, m, d, h, mi, s = _jdToGregorian(self._dayNumber, self._dayFraction, timezone)
        secondRound = round(s, n)
        # It's possible for the seconds place to round to 60, so we need to bump everything.
        if secondRound == 60:
            y, m, d, h, mi, _ = _jdToGregorian(self._dayNumber, self._dayFraction + (1 / SECONDS_PER_DAY), timezone)
            secondRound = 0 if n is None else 0.0
        # force a leading zero in values < 0
        monthString = str(m) if m >= 10 else '0' + str(m)
        dayString = str(d) if d >= 10 else '0' + str(d)
        hourString = str(h) if h >= 10 else '0' + str(h)
        minuteString = str(mi) if mi >= 10 else '0' + str(mi)
        secondString = str(secondRound) if secondRound >= 10 else '0' + str(secondRound)
        # force a leading + on positive timezone offsets
        timezoneString = str(timezone) if timezone < 0 else '+' + str(timezone)

        return f'{y}/{monthString}/{dayString} {hourString}:{minuteString}:{secondString} {timezoneString} UTC'

    def day(self, timezone: float = None) -> str:
        """Return the day portion of the date represented as a string as mm/dd/year."""

        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = None, n: int = 3) -> str:
        """Return the time portion of the date represented as a string as hh:mm:ss +/- tz UTC."""

        return self.date(timezone, n).split(' ')[1]

    def toDatetime(self) -> datetime.datetime:
        """Converts the JulianDate to a Python datetime.datetime object."""

        year, month, day, hour, minutes, secondsFloat = _jdToGregorian(self._dayNumber, self._dayFraction,
                                                                       self._timezone)

        # Split seconds into seconds and microseconds and convert to integers.
        secondsWhole, microSeconds = divmod(secondsFloat, 1.0)
        microSeconds = int(microSeconds * 1000000)

        # Get timezone from datetime module.
        timezone = datetime.timezone(datetime.timedelta(hours=self._timezone))

        return datetime.datetime(year, month, day, hour, minutes, int(secondsWhole), microSeconds, timezone)

    def dayOfYear(self):
        y, m, d, _, _, _ = _jdToGregorian(self._dayNumber, self._dayFraction, self._timezone)

        n1 = int(275 * m / 9)
        n2 = int((m + 9) / 12)
        n3 = (1 + int((y - 4 * int(y / 4) + 2) / 3))

        return n1 - (n2 * n3) + d - 30


def now(timezone: float = None) -> JulianDate:
    """Create a JulianDate object with the current time. If timezone is None, try and get local timezone from Python
    time module."""

    # Use the time module to compute the timezone offset of the local machine.
    if timezone is None:
        gmtoff = time.localtime().tm_gmtoff or 0
        timezone = gmtoff / 3600.0

    delta = datetime.timedelta(hours=timezone)
    tz = datetime.timezone(delta)
    tm = datetime.datetime.now(tz)

    return JulianDate.fromDatetime(tm)


# JulianDate instance of the J2000 epoch.
J2000 = JulianDate(2000, 1, 1, 12, 0, 0, timezone=0)
