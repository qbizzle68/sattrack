import json as _json
import time as _time
import datetime as _datetime
from operator import index as _index

__all__ = ("JulianDate", "now", "J2000")

# -1 is placeholder for indexing
_DAY_IN_MONTH = [-1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _is_leap(year):
    # leap year every 4 years, but not every 100, except for every 400
    return (year % 4 == 0) and (year % 100 != 0 or year % 400 == 0)


def _get_month_days(month, year):
    if month == 2 and _is_leap(year):
        return 29
    return _DAY_IN_MONTH[month]


def _check_args(month, day, year, hour, minute, second):
    month = _index(month)
    day = _index(day)
    year = _index(year)
    hour = _index(hour)
    minute = _index(minute)
    if not isinstance(second, (int, float)):
        raise TypeError('second parameter must be an int or float type,'
                        ' not %s' % type(second).__name__)
    if not 1 <= month <= 12:
        raise ValueError('month must be between 1-12', month)
    monthDays = _get_month_days(month, year)
    if not 1 <= day <= monthDays:
        raise ValueError('day must be between 1-%d' % monthDays, day)
    if not 0 <= hour <= 23:
        raise ValueError('hour must be between 0-23', hour)
    if not 0 <= minute <= 59:
        raise ValueError('minute must be between 0-59', minute)
    if not 0 <= second < 60:
        raise ValueError('second must be between 0-59.999', second)
    return month, day, year, hour, minute, second


def _check_timezone(timezone):
    if not isinstance(timezone, (int, float)):
        raise TypeError('timezone parameter must be an int or float type')


def _jd_to_gregorian(value, timezone):
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
    """
    Class implementation of a Julian Date, which represents a moment in time as an integer day number, and a fraction of
    said day after 12 noon. The JulianDate can be set from the components of a Gregorian calendar date via the __init__
    method, or to the current time using the now() module function. The future() and difference() method allow the
    ability to perform math between different dates, or to compute future or past dates. The solar day is the unit of
    time for these methods, which is equal to 24 hours, 1440 minutes or 86400 seconds.

    Constructors:

    __new__()
    fromNumber()
    fromDatetime()

    Operators:
    __repr__, __str__, __iter__, __reduce__, __int__
    __eq__, __ne__, __lt__, __le__, __gt__, __ge__, __hash__
    __add__, __radd__, __sub__ (add/radd only with timedelta arg, sub only with juliandate)

    Methods:

    setTime()
    future()
    difference()
    date()
    day()
    time()
    toDatetime()

    Properties (readonly):
    value, number, fraction, timezone
    """

    __slots__ = '_dayNumber', '_dayFraction', '_timezone', '_hashcode'

    def __new__(cls, month: int, day: int, year: int, hour: int, minute: int, second: float, timezone: int = 0):
        # constructor
        self = object.__new__(cls)
        self._hashcode = -1
        month, day, year, hour, minute, second = _check_args(month, day, year, hour, minute, second)
        _check_timezone(timezone)
        self._date_to_jd(month, day, year, hour, minute, second, timezone)
        return self

    # other constructors (class methods)

    @classmethod
    def fromNumber(cls, number: float, timezone: float = 0.0) -> 'JulianDate':
        if isinstance(number, (int, float)):
            _check_timezone(timezone)
            # create an 'empty' object with values in valid ranges
            rtn = cls(1, 1, 0, 0, 1, 1, timezone)
            rtn._dayNumber = int(number)
            rtn._dayFraction = number - rtn._dayNumber
            rtn._timezone = timezone
            return rtn
        raise TypeError('number parameter must be an int or float type')

    @classmethod
    def fromDatetime(cls, date: _datetime.datetime):
        if not isinstance(date, _datetime.datetime):
            raise TypeError('fromDatetime() argument must be a datetime instance')
        # check if date is aware or naive
        if date.tzinfo is None:
            tz = 0
        else:
            tz = date.tzinfo.utcoffset(None) / _datetime.timedelta(hours=1)
        seconds = date.second + date.microsecond / 1e6
        return JulianDate(date.month, date.day, date.year, date.hour, date.minute, seconds, tz)

    def __iter__(self):
        yield from {
            'dayNumber': self._dayNumber,
            'dayFraction': self._dayFraction
        }.items()

    def __str__(self) -> str:
        return str(round(self.value, 6)) + ' --- ' + self.date(self._timezone)

    def __repr__(self) -> str:
        # represent as a JSON
        return _json.dumps(dict(self))

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

    # numeric operators
    def future(self, days: int | float):
        # positive is forward in time, negative is backward in time
        if isinstance(days, (int, float)):
            return JulianDate.fromNumber(self.value + days, self._timezone)
        raise TypeError('days parameter must be an int or float type')

    def difference(self, jd: int | float) -> float:
        if isinstance(jd, (int, float)):
            # if self > jd return will be positive, if self < jd return will be negative
            return self.value - jd
        raise TypeError('jd parameter must be an int or float type')

    def date(self, timezone: float = None) -> str:
        if timezone is None:
            timezone = self._timezone
        else:
            _check_timezone(timezone)

        m, d, y, h, mi, s = _jd_to_gregorian(self.value, timezone)
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
        # day of the date is first chunk returned by date
        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = None) -> str:
        # time of the date is second chunk returned by date
        return self.date(timezone).split(' ')[1]

    def toDatetime(self):
        month, day, year, hour, minutes, secondsFloat = _jd_to_gregorian(self.value, self._timezone)
        # split seconds into seconds and microseconds and convert to integers
        secondsInt = int(secondsFloat)
        microSeconds = round((secondsFloat - secondsInt) * 1000000)
        timezone = _datetime.timezone(_datetime.timedelta(hours=self._timezone))
        return _datetime.datetime(year, month, day, hour, minutes, secondsInt, microSeconds, timezone)

    def _date_to_jd(self, month: int, day: int, year: int, hour: int, minute: int, second: float, timezone: int = 0):
        # month, day, year, hour, minute, second = _check_args(month, day, year, hour, minute, second)
        # if not isinstance(timezone, (int, float)):
        #     raise TypeError('timezone parameter must be an int or float type')

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
    if timezone is None:
        tz = _time.localtime().tm_gmtoff / 3600.0
    # elif isinstance(timezone, (int, float)):
        # raise TypeError('timezone parameter must be int or float type')
    else:
        _check_timezone(timezone)
        tz = timezone
    tm = _datetime.datetime.now(_datetime.timezone(_datetime.timedelta(hours=tz)))
    return JulianDate.fromDatetime(tm)


"""J2000 epoch"""
J2000 = JulianDate(1, 1, 2000, 12, 0, 0, timezone=0)
