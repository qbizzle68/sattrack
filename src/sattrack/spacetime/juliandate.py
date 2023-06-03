import json
import time as _time
import datetime as _datetime
from operator import index as _index

if __debug__ is True:
    debug__all__ = ['_DAY_IN_MONTH', '_isLeap', '_getMonthDays', '_checkArgs', '_checkTimezone',
                    '_jdToGregorian']
else:
    debug__all__ = []

__all__ = ["JulianDate", "now", "J2000"] + debug__all__

# if __debug__ is True:
#     debug__all__ = ['_DAY_IN_MONTH', '_isLeap', '_getMonthDays', '_checkArgs', '_checkTimezone',
#                     '_jdToGregorian']
#     __all__ += debug__all__

# -1 is placeholder for indexing
_DAY_IN_MONTH = [-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _isLeap(year):
    """True if year is a leap year, false if not. A leap year is every 4 years, but not every 100 years, except
    for every 400 years."""
    return (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)


def _getMonthDays(month, year):
    """Returns the number of days in a month, using year to determine the number of days in February."""
    if month == 2 and _isLeap(year):
        return 29
    return _DAY_IN_MONTH[month]


def _checkArgs(month, day, year, hour, minute, second):
    """Raises exceptions if the argument types are not valid, returns arguments on success."""
    month = _index(month)
    day = _index(day)
    year = _index(year)
    hour = _index(hour)
    minute = _index(minute)
    if not isinstance(second, (int, float)):
        raise TypeError('second parameter must be an int or float type,'
                        ' not (%s)' % type(second).__name__)
    if not 1 <= month <= 12:
        raise ValueError('month must be between 1-12', month)
    monthDays = _getMonthDays(month, year)
    if not 1 <= day <= monthDays:
        raise ValueError('day must be between 1-%d' % monthDays, day)
    if not 0 <= hour <= 23:
        raise ValueError('hour must be between 0-23', hour)
    if not 0 <= minute <= 59:
        raise ValueError('minute must be between 0-59', minute)
    if not 0 <= second < 60:
        raise ValueError('second must be between 0-59.999', second)
    return month, day, year, hour, minute, second


def _checkTimezone(timezone):
    """Raises an exception if timezone is not a valid type for a timezone."""
    if not isinstance(timezone, (int, float)):
        raise TypeError('timezone parameter must be an int or float type', type(timezone))


def _jdToGregorian(value, timezone):
    """Convert a Julian date to Gregorian calendar components."""
    # incorporate timezone to offset Julian date value
    value += (timezone / 24.0) + 0.5
    # algorithm taken from the 'Solar Position Algorithm for Solar Radiation Applications' paper in appendix A.3
    Z = int(value)
    F = value - Z
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
    s = round(F * 86400.0, 3)
    h = int(s / 3600.0)
    s -= h * 3600.0
    mi = int(s / 60.0)
    s -= mi * 60.0

    return m, d, y, h, mi, s


class JulianDate:
    """Julian Date object which represents a moment in time as an integer day number and a fraction of the day after
    12 noon. An object can be set with Gregorian calendar components via the class constructor, or with a known day
    number and fraction via the fromNumber() class method. See now() for getting a JulianDate object with the current
    time."""

    __slots__ = '_dayNumber', '_dayFraction', '_timezone', '_hashcode'

    def __new__(cls, month: int, day: int, year: int, hour: int, minute: int, second: float, timezone: int = 0):
        """Creates a new JulianDate object from the Gregorian date components and optional timezone offset (negative
        being west of Greenwich, England, positive east)."""
        self = object.__new__(cls)
        self._hashcode = -1
        month, day, year, hour, minute, second = _checkArgs(month, day, year, hour, minute, second)
        _checkTimezone(timezone)
        self._dateToJd(month, day, year, hour, minute, second, timezone)
        return self

    @classmethod
    def fromNumber(cls, number: float, timezone: float = 0.0) -> 'JulianDate':
        """Creates a new JulianDate object directly from a Julian day number with optional timezone offset."""
        if isinstance(number, (int, float)):
            _checkTimezone(timezone)
            # create an 'empty' object with values in valid ranges
            rtn = cls(1, 1, 0, 0, 1, 1, timezone)
            rtn._dayNumber = int(number)
            rtn._dayFraction = number - rtn._dayNumber
            rtn._timezone = timezone
            return rtn
        raise TypeError('number parameter must be an int or float type')

    @classmethod
    def fromDatetime(cls, date: _datetime.datetime):
        """Creates a new JulianDate object from a Python datetime.datetime instance."""
        if not isinstance(date, _datetime.datetime):
            raise TypeError('fromDatetime() argument must be a datetime instance')
        # check if date is aware or naive
        if date.tzinfo is None:
            tz = 0
        else:
            tz = date.tzinfo.utcoffset(None) / _datetime.timedelta(hours=1)
        seconds = date.second + date.microsecond / 1e6
        return JulianDate(date.month, date.day, date.year, date.hour, date.minute, seconds, tz)

    def __str__(self) -> str:
        """Creates a string representation of the JulianDate."""
        return str(round(self.value, 6)) + ' --- ' + self.date(self._timezone)

    def __repr__(self) -> str:
        """Creates a string representation of the JulianDate."""
        m, d, y, h, mi, s = _jdToGregorian(self.value, self._timezone)
        return f'JulianDate({m}, {d}, {y}, {h}, {mi}, {s}, {self._timezone})'

    def toJson(self):
        """Returns a string of the coordinate in json format."""
        return json.dumps(self, default=lambda o: o.toDict())

    def toDict(self):
        """Returns a dictionary of the JulianDate to create json formats of other types containing a GeoPosition."""
        return {"dayNumber": self._dayNumber, "dayFraction": self._dayFraction}

    # operators
    def __add__(self, other: _datetime.timedelta):
        if isinstance(other, _datetime.timedelta):
            deltaSolarDays = other / _datetime.timedelta(days=1)
            return JulianDate.fromNumber(self.value + deltaSolarDays, self._timezone)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other: 'JulianDate') -> float:
        if isinstance(other, JulianDate):
            return self.value - other.value
        return NotImplemented

    def __eq__(self, other: 'JulianDate'):
        if isinstance(other, JulianDate):
            return self.value == other.value
        return NotImplemented

    def __ne__(self, other: 'JulianDate'):
        if isinstance(other, JulianDate):
            return self.value != other.value
        return NotImplemented

    def __lt__(self, other: 'JulianDate'):
        if isinstance(other, JulianDate):
            return self.value < other.value
        return NotImplemented

    def __le__(self, other: 'JulianDate'):
        if isinstance(other, JulianDate):
            return self.value <= other.value
        return NotImplemented

    def __gt__(self, other: 'JulianDate'):
        if isinstance(other, JulianDate):
            return self.value > other.value
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, JulianDate):
            return self.value >= other.value
        return NotImplemented

    def __hash__(self):
        if self._hashcode == -1:
            self._hashcode = hash((self._dayNumber, self._dayFraction))
        return self._hashcode

    def __reduce__(self):
        return self.__class__.fromNumber, (self._dayNumber + self._dayFraction, self._timezone)

    # read-only properties
    @property
    def value(self):
        return self._dayNumber + self._dayFraction

    @property
    def number(self):
        return self._dayNumber

    @property
    def fraction(self):
        return self._dayFraction

    @property
    def timezone(self):
        return self._timezone

    def future(self, days: int | float):
        """Create a new JulianDate object in the future or past relative to the calling instance in solar days. A
        positive value moves forward in time, negative is backward."""
        if isinstance(days, (int, float)):
            return JulianDate.fromNumber(self.value + days, self._timezone)
        raise TypeError('days parameter must be an int or float type')

    def difference(self, jd: int | float) -> float:
        """Compute the difference between a JulianDate and a Julian date number. This is a shorthand for
        (self.value - jd.value)."""
        if isinstance(jd, (int, float)):
            return self.value - jd
        raise TypeError('jd parameter must be an int or float type')

    def date(self, timezone: float = None) -> str:
        """Returns a string with the Julian date represented as a string as mm/dd/year hh:mm:ss +/- tz UTC."""
        if timezone is None:
            timezone = self._timezone
        else:
            _checkTimezone(timezone)

        m, d, y, h, mi, s = _jdToGregorian(self.value, timezone)
        secondRound = round(s, 3)
        # force a leading zero in values < 0
        dayString = str(d) if d >= 10 else '0' + str(d)
        hourString = str(h) if h >= 10 else '0' + str(h)
        minuteString = str(mi) if mi >= 10 else '0' + str(mi)
        secondString = str(secondRound) if secondRound >= 10 else '0' + str(secondRound)
        # force a leading + on positive timezone offsets
        timezoneString = str(timezone) if timezone < 0 else '+' + str(timezone)

        return f'{m}/{dayString}/{y} {hourString}:{minuteString}:{secondString} {timezoneString} UTC'

    def day(self, timezone: float = None) -> str:
        """Return the day portion of the date represented as a string as mm/dd/year."""
        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = None) -> str:
        """Return the time portion of the date represented as a string as hh:mm:ss +/- tz UTC."""
        return self.date(timezone).split(' ')[1]

    def toDatetime(self):
        """Converts the JulianDate to a Python datetime.datetime object."""
        month, day, year, hour, minutes, secondsFloat = _jdToGregorian(self.value, self._timezone)
        # split seconds into seconds and microseconds and convert to integers
        secondsInt = int(secondsFloat)
        microSeconds = round((secondsFloat - secondsInt) * 1000000)
        timezone = _datetime.timezone(_datetime.timedelta(hours=self._timezone))
        return _datetime.datetime(year, month, day, hour, minutes, secondsInt, microSeconds, timezone)

    def _dateToJd(self, month: int, day: int, year: int, hour: int, minute: int, second: float, timezone: int = 0):
        """Logic to compute the Julian day number and fraction from Gregorian date components."""

        # use methodology from the Julian date wiki to convert from Gregorian date
        term0 = int((1461 * (year + 4800 + int((month - 14) / 12))) / 4)
        term1 = int((367 * (month - 2 - (12 * int((month - 14) / 12)))) / 12)
        term2 = int((3 * int((year + 4900 + int((month - 14) / 12)) / 100)) / 4)
        self._dayNumber = term0 + term1 - term2 + day - 32075
        self._dayFraction = ((hour - 12) / 24.0) + (minute / 1440.0) + (second / 86400.0) - (timezone / 24.0)
        if self._dayFraction >= 1.0:
            self._dayNumber += 1
            self._dayFraction -= 1.0
        elif self._dayFraction < 0:
            self._dayNumber -= 1
            self._dayFraction += 1.0
        self._timezone = timezone


def now(timezone=None) -> JulianDate:
    """Create a JulianDate object with the current time. If timezone is None, try and get local timezone from Python
    time module."""
    if timezone is None:
        tz = _time.localtime().tm_gmtoff / 3600.0
    else:
        _checkTimezone(timezone)
        tz = timezone
    tm = _datetime.datetime.now(_datetime.timezone(_datetime.timedelta(hours=tz)))
    return JulianDate.fromDatetime(tm)


'''J2000 epoch'''
J2000 = JulianDate(1, 1, 2000, 12, 0, 0, timezone=0)
