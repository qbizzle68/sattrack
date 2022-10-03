import json
import time
from datetime import datetime


class JulianDate:
    """
    Class implementation of a Julian Date, which represents a moment in time as an integer day number, and a fraction of
    said day after 12 noon. The JulianDate can be set from the components of a Gregorian calendar date via the __init__
    method, or to the current time using the now() module function. The future() and difference() method allow the
    ability to perform math between different dates, or to compute future or past dates. The solar day is the unit of
    time for these methods, which is equal to 24 hours, 1440 minutes or 86400 seconds.
    """

    def __init__(self, month: int, day: int, year: int, hour: int, minute: int, second: float, timeZone: int = 0):
        """
        Initialize an instance from the components of a Gregorian date.

        Args:
            month: Month of the year
            day: Day of the month
            year: Year of the date
            hour: Hours of the day
            minute: Minutes of the hour
            second: Second of the minute
            timeZone: Timezone offset of the time (default = 0)
        """
        self._dayNumber = 0
        self._dayFraction = 0
        self._timezone = timeZone
        self.setTime(month, day, year, hour, minute, second, timeZone)

    @classmethod
    def fromNumber(cls, number: float, timeZone: float = 0.0):
        """Creates a JulianDate instance directly from the Julian date number."""
        rtn = cls(0, 0, 0, 0, 0, 0, timeZone)
        rtn._dayNumber = int(number)
        rtn._dayFraction = number - rtn._dayNumber
        return rtn

        # rtn._day_number = int(day)
        # rtn._timezone = timezone
        # rtn._day_fraction = fraction
        # return rtn

    def __iter__(self):
        """Returns a generator object to iterate through the Julian date values."""
        yield from {
            'dayNumber': self._dayNumber,
            'dayFraction': self._dayFraction
        }.items()

    def __str__(self) -> str:
        """Returns the value and Gregorian date of the JulianDate as a string."""
        return str(round(self.value(), 6)) + ' --- ' + self.date(self._timezone)

    def __repr__(self) -> str:
        """Returns a JSON like representation of the attribute values."""
        return json.dumps(dict(self))

    def __int__(self):
        """Returns the day number of the Julian date."""
        return self._dayNumber

    def __add__(self, rhs: int | float):
        """
        Returns a future or past JulianDate relative to this instance.

        Args:
            rhs: The number of solar days in the future (positive) or past (negative).

        Returns:
            A new JulianDate object with the rhs argument added to it.
        """
        if type(rhs) is int or type(rhs) is float:
            val = self._day_number + self._day_fraction + rhs
            rtn = JulianDate(0, 0, 0, 0, 0, 0, self._timezone)
            rtn._day_number = int(val)
            rtn._day_fraction = val - rtn._day_number
            return rtn
        else:
            raise NotImplemented

    def __sub__(self, other) -> float:
        if type(other) is JulianDate:
            return self._dayNumber + self._dayFraction - other.value()
        elif type(other) is int or type(other) is float:
            return self.__add__(-other)
        else:
            raise NotImplemented

    def __lt__(self, other):
        if type(other) is int or type(other) is float:
            val = other
        elif type(other) is JulianDate:
            val = other.value()
        else:
            raise NotImplemented
        return (self._dayNumber + self._dayFraction) < val

    def __le__(self, other):
        if type(other) is int or type(other) is float:
            val = other
        elif type(other) is JulianDate:
            val = other.value()
        else:
            raise NotImplemented
        return (self._dayNumber + self._dayFraction) <= val

    def __gt__(self, other):
        if type(other) is int or type(other) is float:
            val = other
        elif type(other) is JulianDate:
            val = other.value()
        else:
            raise NotImplemented
        return (self._dayNumber + self._dayFraction) > val

    def __ge__(self, other):
        if type(other) is int or type(other) is float:
            val = other
        elif type(other) is JulianDate:
            val = other.value()
        else:
            raise NotImplemented
        return (self._dayNumber + self._dayFraction) >= val

    def setTime(self, month: int, day: int, year: int, hour: int, minute: int, sec: float, timeZone: int = 0):
        """
        Sets the value of the JulianDate from the components of a Gregorian date.

        Args:
            month: Month of the year
            day: Day of the month
            year: Year of the date
            hour: Hours of the day
            minute: Minutes of the hour
            sec: Second of the minute
            timeZone: Timezone offset of the time (default = 0)
        """

        t0 = int((1461 * (year + 4800 + int((month - 14) / 12))) / 4)
        t1 = int((367 * (month - 2 - (12 * int((month - 14) / 12)))) / 12)
        t2 = int((3 * int((year + 4900 + int((month - 14) / 12)) / 100)) / 4)
        self._dayNumber = t0 + t1 - t2 + day - 32075
        self._dayFraction = ((hour - 12) / 24.0) + (minute / 1440.0) + (sec / 86400.0) - (timeZone / 24.0)
        if self._dayFraction >= 1.0:
            self._dayNumber += 1
            self._dayFraction -= 1.0
        elif self._dayFraction < 0:
            self._dayNumber -= 1
            self._dayFraction += 1.0
        self._timezone = timeZone

    def value(self):
        """Returns the value of the JulianDate."""
        return self._dayNumber + self._dayFraction

    def getNumber(self):
        """Returns the Julian day number."""
        return self._dayNumber

    def getFraction(self):
        """Returns the fraction part of the JulianDate."""
        return self._dayFraction

    def future(self, days: float):
        """
        Returns a future or past date, relative to the calling date. A positive argument computes a future date, a
        negative value computes a past date.

        Args:
            days: The offset between dates in solar days, where a positive value indicates forward, negative backwards.

        Returns:
            The future or past JulianDate.
        """

        rtn = JulianDate(0, 0, 0, 0, 0, 0)
        int_day = int(days)
        rtn._dayNumber = self._dayNumber + int_day
        rtn._dayFraction = self._dayFraction + (days - int_day)
        if rtn._dayFraction >= 1.0:
            rtn._dayNumber += 1
            rtn._dayFraction -= 1.0
        rtn._timezone = self._timezone
        return rtn

    def difference(self, jd) -> float:
        """
        Returns the difference in time between two JulianDates. The value is relative to the calling object, so if the
        argument value is before the calling object, the returned value is positive. The opposite returns a negative
        value.

        Args:
            jd: Relative JulianDate to compute the difference between.

        Returns:
            The difference between the calling date and the argument date in solar days.
        """

        return self.value() - jd.value()

    def date(self, timeZone: float = None) -> str:
        """
        Generates a string with the Gregorian calendar date equivalent of the calling date. The timezone offset adjusts
        to the correct time. Printing this string allows the printed string to accommodate a timezone offset, something
        the __str__ method can't do.

        Args:
            timeZone: Timezone offset if the desired time is not UTC. (Default = 0)

        Returns:
            A string representing the date as a Gregorian calendar date.
        """

        if timeZone is None:
            timeZone = self._timezone

        val = self._dayNumber + self._dayFraction + 0.5 + (timeZone / 24.0)
        Z = int(val)
        F = val - Z
        if Z < 2299161:
            A = Z
        else:
            B = int((Z - 1867216.25) / 36524.25)
            A = Z + 1 + B - int(B / 4)
        C = A + 1524
        D = int((C - 122.1) / 365.25)
        G = int(365.25 * D)
        I = int((C - G) / 30.6001)
        d = C - G - int(30.6001 * I) + F
        m = I - 1 if I < 14 else I - 13
        y = D - 4716 if m > 2 else D - 4715

        frac = val - int(val)
        s = round(frac * 86400.0, 3)
        h = int(s / 3600.0)
        s -= h * 3600.0
        mi = int(s / 60.0)
        s -= mi * 60.0

        return ('{}/{}/{} {}:{}:{} {} UTC'
                .format(m, int(d) if d >= 10 else '0' + str(int(d)),
                        y, h if h >= 10 else '0' + str(h),
                        mi if mi >= 10 else '0' + str(mi),
                        str(round(s, 3)) if s >= 10 else '0' + str(round(s, 3)),
                        str(timeZone) if timeZone < 0 else '+' + str(timeZone)
                        ))

    def day(self, timezone: float = None) -> str:
        """
        Returns a string with the day portion of the Gregorian calendar date equivalent.

        Args:
            timezone: Timezone offset to adjust the day number if applicable. (Default = 0.0)

        Returns:
            A string of the day portion of the Gregorian calendar date equivalent.
        """

        if timezone is None:
            timezone = self._timezone

        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = None) -> str:
        """
        Returns a string with the time portion of the Gregorian calendar date equivalent.

        Args:
            timezone: Timezone offset to adjust the time. (Default = 0.0)

        Returns:

        """

        if timezone is None:
            timezone = self._timezone

        return self.date(timezone).split(' ')[1]

    def getTimeZone(self):
        return self._timezone


def now() -> JulianDate:
    """Returns a new JulianDate object equal to the current UTC time."""
    # tm = datetime.utcnow()
    tm = datetime.now()
    return JulianDate(tm.month, tm.day, tm.year, tm.hour, tm.minute, tm.second + (tm.microsecond / 1000000),
                      time.localtime().tm_gmtoff / 3600.0)


J2000 = JulianDate(1, 1, 2000, 12, 0, 0)
