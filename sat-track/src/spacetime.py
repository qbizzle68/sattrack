from datetime import datetime


class JulianDate:

    def __init__(self, month: int, day: int, year: int, hour: int, minute: int, second: float, timeZone: int = 0):
        self._date = 0
        self.setTime(month, day, year, hour, minute, second, timeZone)

    def __str__(self) -> str:
        """Returns the value and Gregorian date of the JulianDate."""
        return str(round(self._date, 6)) + ' --- ' + self.date()

    def setTime(self, m: int, d: int, y: int, hr: int, mi: int, sec: float, timeZone: int = 0):
        """Sets the value of a JulianDate object based on a Gregorian calendar date.
        Parameters:
        m:        Gregorian date month.
        d:        Gregorian date day.
        y:        Gregorian date year.
        hr:       Gregorian date hour.
        mi:       Gregorian date minute.
        sec:      Gregorian date second, includes fractional portions.
        timeZone: Timezone of the Gregorian date time, defaults to 0."""
        t0 = int((1461 * (y + 4800 + int((m - 14) / 12))) / 4)
        t1 = int((367 * (m - 2 - (12 * int((m - 14) / 12)))) / 12)
        t2 = int((3 * int((y + 4900 + int((m - 14) / 12)) / 100)) / 4)
        self._date = t0 + t1 - t2 + d - 32075 + (
                    ((hr - 12) / 24.0) + (mi / 1440.0) + (sec / 86400.0) - (timeZone / 24.0))

    def value(self):
        """Returns the value of the JulianDate."""
        return self._date

    def number(self):
        """Returns the Julian day number of this JulianDate."""
        return int(self._date)

    def fraction(self):
        """Returns the fractional part of this JulianDate."""
        return self._date - (int(self._date))

    def future(self, days: float):
        """Computes a future JulianDate from this date.
        Parameters:
        days:   Number of solar days in the future of this date. A negative value
                indicates a computing a past date."""
        rtn = JulianDate(0, 0, 0, 0, 0, 0)
        rtn._date = self._date + days
        return rtn

    def difference(self, jd) -> float:
        """Returns the difference between two JulianDates, measured in solar days.
        Parameters:
        jd: A JulianDate object to find the difference from this object."""
        return self._date - jd._date

    def date(self, timezone: float = 0.0) -> str:
        """Generates a string with the Gregorian date equivalent of this JulianDate.
        Parameters:
        timezone: The timezone offset to create the Gregorian date."""
        Z = int(self._date + 0.5 + (timezone / 24.0))
        F = (self._date + 0.5 + (timezone / 24.0)) - Z
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
        """Returns a string with the day portion of the Gregorian date equivalent
        of this JulianDate.
        Parameters:
        timezone: The timezone used to generate the Gregorian date, defaults to 0."""
        return self.date(timezone).split(' ')[0]

    def time(self, timezone: float = 0.0) -> str:
        """Returns a string with the time portion of the Gregorian date equivalent
        of this JulianDate.
        Parameters:
        timezone: The timezone used to generate the Gregorian date, defaults to 0."""
        return self.date(timezone).split(' ')[1]


def now() -> JulianDate:
    """Returns a new JulianDate object equal to the current UTC time."""
    tm = datetime.utcnow()
    return JulianDate(tm.month, tm.day, tm.year, tm.hour, tm.minute, tm.second + (tm.microsecond / 1000000))


J2000 = JulianDate(1, 1, 2000, 12, 0, 0)
