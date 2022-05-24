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


def siderealTime(jd: JulianDate) -> float:
    """Converts a certain solar time into sidereal time, most useful for computing Earth's
    geographic reference frame offset from the geocentric equitorial celestial reference frame.
    I.e. this is the offset angle of GMT to the first point of Aries measured in hours.
    Parameters:
    jd: The time to convert.
    return: The local sidereal time of Greenwich, England."""
    dt = jd.difference(J2000)
    tmp = (18.697_374_558 + 24.065_709_824_419_08 * dt) % 24.0
    return tmp + 24.0 if tmp < 0 else tmp


def localSiderealTime(jd: JulianDate, lng: float) -> float:
    """Converts a certain solar time into local sidereal time, useful for computing your
    longitudinal offset from the geocentric equitorial celestial reference frame.
    Parameters:
    jd:     The time to convert.
    lng:    The local longitude.
    Returns the local sidereal time in hours."""
    if lng < 0:
        lng += 360.0
    lng2RA = lng / 360.0 * 24.0
    lst = siderealTime(jd) + lng2RA
    return lst % 24.0


def earthOffsetAngle(jd: JulianDate) -> float:
    """Computes the earth offset angle from the celestial reference frame for a
    given solar time.
    Parameters:
    jd: The time to compute.
    Returns the offset angle in degrees."""
    return siderealTime(jd) / 24.0 * 360.0
