import datetime
import json
import re
import time
import warnings
from collections import namedtuple
from functools import cached_property
from io import StringIO
from typing import Union

DateComponents = namedtuple('DateComponents', 'year month day hour minute seconds')
_HOUR = datetime.timedelta(hours=1)


def _isLeap(year: int) -> bool:
    """Returns True if the year is a leap year, else false. All negative multiples of
    4 returns True."""

    if year >= 0:
        return year % 4 == 0 and year % 100 != 0 or year % 400 == 0
    # All multiples of 4 are leap years when year < 0 according to the Julian Date
    return year % 4 == 0


def _julianDateToComponents(number: int, fraction: float, timezone: float) -> DateComponents:
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

    return DateComponents(y, m, d, h, mi, s)


def _componentsToJulianDate(year: int, month: int, day: int, hour: int, minute: int, second: float,
                            timezone: float = 0) -> (int, float):
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

    return int(dayNumber), dayFrac


def _findDayOfWeek(julianDate: float, timezone: float, standard: str = 'us') -> int:
    """Returns the day of the week as an integer. If standard is 'us' the value is in [0-6] where 0 is
    Sunday and 6 is Saturday. If standard is 'iso' the value is in [1-7] where Sunday is instead equal to 7.
    Any other value for standard raises a ValueError."""

    julianValue = julianDate + (timezone / 24)
    number, fraction = divmod(julianValue, 1.0)
    number = int(number)
    if fraction >= 0.5:
        number += 1

    if standard == 'us':
        return (number + 1) % 7
    elif standard == 'iso':  # pragma: no cover
        return (number % 7) + 1
    else:  # pragma: no cover
        raise ValueError('standard must be \'us\' or \'iso\'')


class DateFormatter:

    def __init__(self, components: DateComponents, timezone: datetime.tzinfo):
        self._components = components
        self._proxyDatetime = components.year <= 0
        self._timezone = timezone

    @property
    def components(self) -> DateComponents:  # pragma: no cover
        return self._components

    @property
    def proxy(self) -> bool:
        return self._proxyDatetime

    @property
    def timezone(self):  # pragma: no cover
        return self._timezone

    @cached_property
    def second(self) -> int:
        return int(self._components.seconds)

    @cached_property
    def microsecond(self) -> int:
        fraction = (self._components.seconds - self.second) * 1e6
        return round(fraction)

    @cached_property
    def datetime(self) -> datetime.datetime:
        components = self._components
        year = components.year
        if year > 0:
            return datetime.datetime(*components[:-1], self.second, self.microsecond, self._timezone)
        else:
            isLeap = _isLeap(year)
            _year = 4 if isLeap else 1
            _datetime = datetime.datetime(_year, *components[1:-1], self.second, self.microsecond,
                                          self._timezone)
            equivalentYear = self._getEquivalentYear(year, True)

            return _datetime.replace(year=equivalentYear)

    @staticmethod
    def _getEquivalentYear(year: int, identifiable: bool = False) -> int:  # pragma: no cover

        julianDate = sum(_componentsToJulianDate(year, 1, 1, 0, 0, 0, 0))
        dayOfWeek = _findDayOfWeek(julianDate, 0, 'us')
        isLeap = _isLeap(year)

        # todo: do we even use the years in the 2000's? should we just use the 6000's?
        match dayOfWeek:
            case 0:
                return (6688 if isLeap else 6665) if identifiable else (2012 if isLeap else 2023)
            case 1:
                return (6672 if isLeap else 6666) if identifiable else (2024 if isLeap else 2018)
            case 2:
                return (6684 if isLeap else 6661) if identifiable else (2008 if isLeap else 2019)
            case 3:
                return (6668 if isLeap else 6662) if identifiable else (2020 if isLeap else 2014)
            case 4:
                return (6680 if isLeap else 6663) if identifiable else (2004 if isLeap else 2015)
            case 5:
                return (6664 if isLeap else 6669) if identifiable else (2016 if isLeap else 2021)
            case 6:
                return (6676 if isLeap else 6670) if identifiable else (2000 if isLeap else 2022)
            case _:
                raise ValueError(f'US day of week is not in [0, 6]')

    def _adjustYearDependentSpec(self, char: str) -> str:
        """Any formatting done here is presumed to be on a date where the year is less than 0,
        i.e. a datetime.datetime object can not be instantiated to handle the formatting, so we
        need to intervene to ensure validity. char should be whatever the value of a single
        format spec after the '%' character."""

        components = self._components
        if char == 'y':
            result = components.year % 100
            result = result - 100 if result > 0 else result
            return f'-{abs(result):0>2}'
        elif char == 'Y':
            return f'-{abs(components.year):0>4}'
        elif char == 'z' or char == 'Z' or char == ':z':
            _datetime = self.datetime.replace(year=4)
            return _datetime.strftime(f'%{char}')
        elif char == 'c':
            # The issue here is if the year includes century or not. If we carefully choose a year, we can
            # ensure the number that appears is unique (not another unit of time) and make the replacement.

            # This choice of year is only valid for the next 4500 years or so. lol
            replacementYear = self._getEquivalentYear(components.year, True)
            _datetime = self.datetime.replace(year=replacementYear)
            result = _datetime.strftime('%c')
            # First, assume the century is included:
            match = re.search(rf'\D*{replacementYear}\D*', result)
            if match:
                year = abs(components.year)
                return result.replace(str(replacementYear), f'-{year:0>4}')
            # fixme: This seems very unlikely to be used, and we might not need it.
            # Otherwise, check for century omitted:
            replacementYear %= 100
            match = re.search(rf'\D*{replacementYear}\D*', result)
            if match:
                year = abs(components.year % 100 - 100)
                if year == 100:
                    year = 0
                return result.replace(str(replacementYear), f'-{year:0>2}')

            # If we end up here, no year is found in the output which should mean there was an error
            raise ValueError('unable to make year substitution in datetime result')  # pragma: no cover
        elif char == 'x':
            # Same approach as %c format spec
            replacementYear = self._getEquivalentYear(components.year, True)
            _datetime = self.datetime.replace(year=replacementYear)
            result = _datetime.strftime('%x')
            match = re.search(rf'\D*{replacementYear}\D*', result)
            if match:
                year = abs(components.year)
                return result.replace(str(replacementYear), f'-{year:0>4}')
            replacementYear %= 100
            match = re.search(rf'\D*{replacementYear}\D*', result)
            if match:
                year = abs(components.year % 100 - 100)
                if year == 100:
                    year = 0
                return result.replace(str(replacementYear), f'-{year:0>2}')

            raise ValueError('unable to make year substitution in datetime result')  # pragma: no cover
        elif char == 'G':
            _datetime = self.datetime
            year = int(_datetime.strftime('%G'))
            difference = _datetime.year - year
            result = components.year - difference
            if result == 0:
                return '0000'
            else:
                return f'-{abs(result):0>4}'

        return self.datetime.strftime(f'%{char}')

    def format(self, formatSpec: str) -> str:
        if not self._proxyDatetime:
            return self.datetime.strftime(formatSpec)
        else:
            buffer = StringIO()
            stringIterator = iter(formatSpec)
            for char in stringIterator:
                if char == '%':
                    try:
                        char = next(stringIterator)
                    except StopIteration:
                        raise ValueError('Invalid format string')

                    buffer.write(self._adjustYearDependentSpec(char))
                else:
                    buffer.write(char)

        return buffer.getvalue()


_DEFAULT_TIMEZONE = datetime.timezone.utc
Timezone = Union[float, int, datetime.tzinfo]


class JulianDate:
    __slots__ = '_number', '_fraction', '_timezone', '_utcOffset', '_components', '_formatter'

    def __init__(self, year: int, month: int, day: int, hour: int, minute: int, seconds: float,
                 timezone: Timezone = _DEFAULT_TIMEZONE):

        self._components = DateComponents(year, month, day, hour, minute, seconds)
        if isinstance(timezone, (int, float)):
            delta = datetime.timedelta(hours=timezone)
            timezone = datetime.timezone(delta)
        self._formatter = DateFormatter(self._components, timezone)
        self._timezone = timezone
        _datetime = self._formatter.datetime
        if year < 1:
            _datetime = _datetime.replace(year=4)
        self._utcOffset = timezone.utcoffset(_datetime) / _HOUR
        self._number, self._fraction = _componentsToJulianDate(year, month, day, hour, minute, seconds, self._utcOffset)

    @classmethod
    def fromNumber(cls, number: int, fraction: float, timezone: Timezone = _DEFAULT_TIMEZONE):
        # we want to allow fraction to be larger than one. as an implementation detail of
        # _julianDateToComponents, this is allowed; if this changes, we need to add more logic here
        # to ensure this assumption stays true

        if isinstance(timezone, (int, float)):
            delta = datetime.timedelta(hours=timezone)
            timezone = datetime.timezone(delta)

        utcComponents = _julianDateToComponents(number, fraction, 0.0)
        second, microsecond = divmod(utcComponents.seconds, 1.0)
        # We need to use a positive leap year to ensure the datetime is valid. Always using
        # a leap year is okay because we just need the utc offset.
        year = utcComponents[0] if utcComponents[0] > 0 else 4
        utcDatetime = datetime.datetime(year, *utcComponents[1:-1], int(second), int(microsecond * 1e6),
                                        datetime.timezone.utc)
        utcOffset = timezone.utcoffset(utcDatetime) / _HOUR
        components = _julianDateToComponents(number, fraction, utcOffset)

        # need to check if a value is extremely close to being rounded up, e.g. a second that is
        # 49.99999999999 instead of 50.0
        # need to check if rounding the seconds increases the integer part
        if int(round(components.seconds, 6)) > int(components.seconds):
            # right now we don't need to adjust number if fraction 'runs over' because it'll get
            # added in _julianDateToComponents
            oneSecond = 1 / 86400
            adjustedComponents = _julianDateToComponents(number, fraction + oneSecond, utcOffset)
            adjustedSeconds = round(components.seconds, 6) % 60
            components = DateComponents(*adjustedComponents[:-1], adjustedSeconds)

        self = object.__new__(cls)
        self.__init__(*components, timezone)
        return self

    @classmethod
    def fromDatetime(cls, _datetime: datetime.datetime):
        self = object.__new__(cls)
        seconds = _datetime.second + (_datetime.microsecond / 1e6)
        self.__init__(_datetime.year, _datetime.month, _datetime.day, _datetime.hour,
                      _datetime.minute, seconds, _datetime.tzinfo)
        return self

    @classmethod
    def now(cls, timezone: Timezone = None):
        if timezone is None:
            gmtoff = time.localtime().tm_gmtoff or 0
            timezone = gmtoff / 3600.0

        if isinstance(timezone, (int, float)):
            delta = datetime.timedelta(hours=timezone)
            tz = datetime.timezone(delta)
        elif isinstance(timezone, datetime.tzinfo):
            tz = timezone
        else:
            raise ValueError(f'timezone must be a Timezone type, not {type(timezone)}')

        _time = datetime.datetime.now(tz)
        return cls.fromDatetime(_time)

    @classmethod
    def fromisoformat(cls, dateString: str):
        _datetime = datetime.datetime.fromisoformat(dateString)
        return cls.fromDatetime(_datetime)

    @classmethod
    def strptime(cls, dateString: str, formatSpec: str):
        _datetime = datetime.datetime.strptime(dateString, formatSpec)
        return cls.fromDatetime(_datetime)

    def asTimezone(self, timezone: Timezone) -> 'JulianDate':
        return JulianDate.fromNumber(self._number, self._fraction, timezone)

    def __format__(self, format_spec: str) -> str:
        return self._formatter.format(format_spec)

    def strftime(self, format_spec: str) -> str:
        return self.__format__(format_spec)

    def __str__(self) -> str:
        _format = f'{self._number + self._fraction} --- %Y-%m-%d %H:%M:%S.%f %z UTC'
        return self.__format__(_format)

    def __repr__(self) -> str:
        componentString = '{}, {}, {}, {}, {}, {}'.format(*self._components)
        return f'{self.__class__.__name__}({componentString}, {self._timezone!r})'

    def __round__(self, n=None) -> 'JulianDate':
        microseconds = float(self.__format__('%f')) / 1e6
        rounded = round(microseconds, n)
        difference = (microseconds - rounded) / 86400

        fraction = self._fraction - difference
        numberAdjust, fraction = divmod(fraction, 1.0)

        return self.fromNumber(self._number + numberAdjust, fraction, self._timezone)

    def date(self, timezone: Timezone = None, n: int = 3) -> str:
        if timezone is not None:
            julianDate = self.asTimezone(timezone)
        else:
            julianDate = self

        if n == 0:
            return julianDate.__round__(0).__format__('%Y-%m-%d %H:%M:%S %z UTC')
        elif n == 6:
            return julianDate.__format__('%Y-%m-%d %H:%M:%S.%f %z UTC')
        elif n < 0:
            raise ValueError(f'n must be a non-negative integer, not {n}')
        # fixme: what if n > 6? should error check better

        formattedString = julianDate.__round__(n).__format__('%Y-%m-%d %H:%M:%S.%f %z UTC')
        microseconds = formattedString.split('.')[1].split(' ')[0]

        return formattedString.replace(microseconds, microseconds[:n])

    def day(self, timezone: float = None) -> str:
        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = None, n: int = 3) -> str:
        return self.date(timezone, n).split(' ')[1]

    @property
    def components(self) -> DateComponents:
        return self._components

    @property
    def number(self) -> int:
        return self._number

    @property
    def fraction(self) -> float:
        return self._fraction

    @property
    def value(self) -> float:
        return self._number + self._fraction

    @property
    def timezone(self) -> datetime.tzinfo:
        return self._timezone

    @property
    def utcOffset(self) -> float:
        return self._utcOffset

    def toDict(self) -> dict:
        return {"number": self._number, "fraction": self._fraction, "utcOffset": self._utcOffset}

    def toJson(self) -> str:
        return json.dumps(self, default=lambda o: o.toDict())

    def __add__(self, other: datetime.timedelta | float | int) -> 'JulianDate':
        if isinstance(other, datetime.timedelta):
            deltaSolarDays = other / datetime.timedelta(days=1)
            numberIncrease, fractionIncrease = divmod(deltaSolarDays, 1.0)
            return JulianDate.fromNumber(self._number + numberIncrease,
                                         self._fraction + fractionIncrease,
                                         self._timezone)
        elif isinstance(other, (int, float)):
            numberIncrease, fractionIncrease = divmod(other, 1.0)
            numberAdjust, fraction = divmod(self._fraction + fractionIncrease, 1.0)
            return JulianDate.fromNumber(self._number + numberIncrease + numberAdjust,
                                         fraction, self._timezone)
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other: 'JulianDate | float | int') -> float:
        if isinstance(other, JulianDate):
            return (self._number - other._number) + (self._fraction - other._fraction)
        elif isinstance(other, (int, float)):
            number, fraction = divmod(other, 1.0)
            return JulianDate.fromNumber(self._number - number, self._fraction - fraction, self._timezone)
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

    def __hash__(self) -> int:
        return hash((self._number, self._fraction, self._timezone))

    def __reduce__(self):
        return self.__class__, (*self._components, self._timezone)

    def toDatetime(self) -> datetime.datetime:
        if self._components.year <= 0:
            raise ValueError('unable to instantiate datetime.datetime object with non-positive year')

        return self._formatter.datetime

    def dayOfYear(self) -> int:
        month = self._components.month
        n1 = (275 * month) // 9
        n2 = (month + 9) // 12
        year = self._components.year
        n3 = 1 + (year - 4 * (year // 4) + 2) // 3

        return n1 - (n2 * n3) + self._components.day - 30

    def future(self, days: float | int) -> 'JulianDate':
        warnings.warn('JulianDate.future() will be removed in future versions, please use the + operator',
                      DeprecationWarning)

        # who cares about efficiency or preserving precision, this will be deprecated soon
        delta = datetime.timedelta(days=days)
        return self + delta

        numberAdjust, fraction = divmod(days + self._fraction, 1.0)
        if numberAdjust < 0:
            numberAdjust += 1
            fraction = 1 - fraction
        return JulianDate.fromNumber(self._number + numberAdjust, fraction, self._timezone)


J2000 = JulianDate(2000, 1, 1, 12, 0, 0, datetime.timezone.utc)


# ========================= OLD CODE BEGINS HERE =====================
# These functions aren't needed now, but we spent all that time figuring out some of them, I don't want to delete yet.


def _findDayOfYear(year: int, month: int, day: int) -> int:  # pragma: no cover
    """Find the day number of the year."""

    n1 = int(275 * month / 9)
    n2 = int((month + 9) / 12)
    n3 = (1 + int((year - 4 * int(year / 4) + 2) / 3))

    return n1 - (n2 * n3) + day - 30


def _findWeekOfYear(year: int, month: int, day: int, utcOffset: float, value: float = None,
                    sundayStart: bool = True) -> int:  # pragma: no cover
    """Computes the week of the year. If sundayStart is True, Sunday is the first day of the week,
    otherwise Monday is the first day of the week. Any days at the beginning of a year before the
    start of the week are considered to be in week 0."""

    if value is None:
        value = sum(_componentsToJulianDate(year, 1, 1, 0, 0, 0, utcOffset))

    # standard = 'us' if sundayStart is True else 'iso'
    yearStartWeekday = _findDayOfWeek(value, utcOffset, 'us')
    dayOfYear = _findDayOfYear(year, month, day)

    # Determine how many days occur during week 0.
    startWeekdayNumber = 0 if sundayStart else 1
    daysInWeekZero = (startWeekdayNumber - yearStartWeekday) % 7
    weekOneDayNumber = daysInWeekZero + 1

    daysSinceWeekOne = dayOfYear - weekOneDayNumber
    return (daysSinceWeekOne // 7) + 1


def _findWeekOfYearNonZero(year: int, month: int, day: int, utcOffset: float,
                           value: float = None) -> int:  # pragma: no cover
    dayOfYear = _findDayOfYear(year, month, day)
    jan1DayOfWeek = _findDayOfWeek(value, utcOffset, 'iso')
    if jan1DayOfWeek <= 4:
        if month == 12 and day >= 29:
            dec31DayOfWeek = _findDayOfWeek(value, utcOffset, 'iso')
            if (dec31DayOfWeek == 1 and day == 31) or \
               (dec31DayOfWeek == 2 and day >= 30) or \
               (dec31DayOfWeek == 3 and day >= 29):
                return 1

        decWeek1Count = jan1DayOfWeek - 1
        weekNumber = ((dayOfYear + decWeek1Count - 1) // 7) + 1
    else:
        janLastYearCount = 8 - jan1DayOfWeek
        if month == 1 and day <= janLastYearCount:
            return _findWeekOfYearNonZero(year - 1, 12, 31, utcOffset, value)

        weekNumber = ((dayOfYear - janLastYearCount - 1) // 7) + 1

    return weekNumber
