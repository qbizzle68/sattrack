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
        self._day_number = 0
        self._day_fraction = 0
        self.setTime(month, day, year, hour, minute, second, timeZone)

    def __str__(self) -> str:
        """Returns the value and Gregorian date of the JulianDate as a string."""
        return str(round(self.value(), 6)) + ' --- ' + self.date()

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
        self._day_number = t0 + t1 - t2 + day - 32075
        self._day_fraction = ((hour - 12) / 24.0) + (minute / 1440.0) + (sec / 86400.0) - (timeZone / 24.0)
        if self._day_fraction >= 1.0:
            self._day_number += 1
            self._day_fraction -= 1.0
        elif self._day_fraction < 0:
            self._day_number -= 1
            self._day_fraction += 1.0

    def value(self):
        """Returns the value of the JulianDate."""
        return self._day_number + self._day_fraction

    def number(self):
        """Returns the Julian day number."""
        return self._day_number

    def fraction(self):
        """Returns the fraction part of the JulianDate."""
        return self._day_fraction

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
        rtn._day_number = self._day_number + int_day
        rtn._day_fraction = self._day_fraction + (days - int_day)
        if rtn._day_fraction >= 1.0:
            rtn._day_number += 1
            rtn._day_fraction -= 1.0
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

    def date(self, timezone: float = 0.0) -> str:
        """
        Generates a string with the Gregorian calendar date equivalent of the calling date. The timezone offset adjusts
        to the correct time. Printing this string allows the printed string to accommodate a timezone offset, something
        the __str__ method can't do.

        Args:
            timezone: Timezone offset if the desired time is not UTC. (Default = 0)

        Returns:
            A string representing the date as a Gregorian calendar date.
        """

        val = self.value()
        Z = int(val + 0.5 + (timezone / 24.0))
        F = (val + 0.5 + (timezone / 24.0)) - Z
        if Z < 2299161:
            A = Z
        else:
            B = int((Z - 1867216.25) / 36524.25)
            A = Z + 1 + B - (B / 4)
        C = A + 1524
        D = int((C - 122.1) / 365.25)
        G = int(365.25 * D)
        _I = int((C - G) / 30.6001)
        d = C - G - int(30.6001 * _I) + F
        m = _I - 1 if _I < 14 else _I - 13
        y = D - 4716 if m > 2 else D - 4715

        dayFrac = d - int(d)
        h = int(dayFrac * 24)
        mi = int((dayFrac - (h / 24.0)) * 1440.0)
        s = (dayFrac - (h / 24.0) - (mi / 1440.0)) * 86400.0

        return ('{}/{}/{} {}:{}:{} {} UTC'
                .format(m, int(d) if d >= 10 else '0' + str(int(d)),
                        y, h if h >= 10 else '0' + str(h),
                        mi if mi >= 10 else '0' + str(mi),
                        str(round(s, 3)) if s >= 10 else '0' + str(round(s, 3)),
                        str(timezone) if timezone < 0 else '+' + str(timezone)
                        ))

    def day(self, timezone: float = 0.0) -> str:
        """
        Returns a string with the day portion of the Gregorian calendar date equivalent.

        Args:
            timezone: Timezone offset to adjust the day number if applicable. (Default = 0.0)

        Returns:
            A string of the day portion of the Gregorian calendar date equivalent.
        """

        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = 0.0) -> str:
        """
        Returns a string with the time portion of the Gregorian calendar date equivalent.

        Args:
            timezone: Timezone offset to adjust the time. (Default = 0.0)

        Returns:

        """

        return self.date(timezone).split(' ')[1]


def now() -> JulianDate:
    """Returns a new JulianDate object equal to the current UTC time."""
    tm = datetime.utcnow()
    return JulianDate(tm.month, tm.day, tm.year, tm.hour, tm.minute, tm.second + (tm.microsecond / 1000000))


J2000 = JulianDate(1, 1, 2000, 12, 0, 0)
