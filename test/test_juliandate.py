import datetime
import unittest
from unittest import mock
from copy import copy

from sattrack.core.juliandate import now, J2000, JulianDate

values = ((2000, 1, 1, 12, 0, 0),
          (1999, 1, 1, 0, 0, 0),
          (1988, 6, 19, 12, 0, 0),
          (1988, 1, 27, 0, 0, 0),
          (1987, 6, 19, 12, 0, 0),
          (1987, 1, 27, 0, 0, 0),
          (1900, 1, 1, 0, 0, 0),
          (1600, 12, 31, 0, 0, 0),
          (1600, 1, 1, 0, 0, 0),
          (837, 4, 10, 7, 12, 0),
          (-122, 1, 1, 0, 0, 0),
          (-123, 12, 31, 0, 0, 0),
          (-1000, 7, 12, 12, 0, 0),
          (-1000, 2, 29, 0, 0, 0),
          (-1001, 8, 17, 21, 36, 0),
          (-4712, 1, 1, 12, 0, 0))

answers = (2451545, 2451179.5, 2447332, 2447187.5, 2446966, 2446822.5, 2415020.5, 2305812.5, 2305447.5,
           2026871.8, 1676497.5, 1676496.5, 1356001, 1355866.5, 1355671.4, 0.0)


class TestJuliandate(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        jdList = []
        for args in values:
            jd = JulianDate(*args)
            jdList.append(jd)

        cls.jds = jdList

    def testValues(self):
        for jd, answer in zip(self.jds, answers):
            with self.subTest(jd=jd, number=answer):
                self.assertEqual(jd.value, answer)

    @unittest.skip('Work in progress')
    def testValueWithTimezone(self):
        pass

    def testString(self):

        results = ('2000/01/01 12:00:00.0', '1999/01/01 00:00:00.0', '1988/06/19 12:00:00.0',
                   '1988/01/27 00:00:00.0', '1987/06/19 12:00:00.0', '1987/01/27 00:00:00.0',
                   '1900/01/01 00:00:00.0', '1600/12/31 00:00:00.0', '1600/01/01 00:00:00.0',
                   '837/04/10 07:12:00.0', '-122/01/01 00:00:00.0', '-123/12/31 00:00:00.0',
                   '-1000/07/12 12:00:00.0', '-1000/02/29 00:00:00.0', '-1001/08/17 21:36:00.0',
                   '-4712/01/01 12:00:00.0')

        for jd, result in zip(self.jds, results):
            with self.subTest(jd=jd.value):
                answer = f'{jd.value} --- {result} +0 UTC'
                self.assertEqual(str(jd), answer)

    def testJson(self):
        answer = '{"dayNumber": 2451545, "dayFraction": 0.0}'
        self.assertEqual(self.jds[0].toJson(), answer)

        answer = '{"dayNumber": 2451179, "dayFraction": 0.5}'
        self.assertEqual(self.jds[1].toJson(), answer)

    def testEquality(self):
        jd1 = JulianDate(2000, 1, 1, 12, 0, 0)
        jd2 = JulianDate(2000, 1, 1, 12, 0, 0)
        self.assertEqual(jd1, jd2)

        # Equality with timezone offset.
        jd3 = JulianDate(2000, 1, 1, 7, 0, 0, -5)
        self.assertEqual(jd1, jd3)

        # Not equal
        self.assertNotEqual(self.jds[0], self.jds[1])
        # With timezone offset.
        self.assertFalse(jd1 != jd3)

    def testComparators(self):
        for i, jd in enumerate(self.jds[:-1]):
            with self.subTest(i=i, jd=jd.value):
                self.assertLess(self.jds[i + 1], jd)
                self.assertLessEqual(self.jds[i + 1], jd)

        for i, jd in enumerate(self.jds[:-1]):
            with self.subTest(i=i, jd=jd.value):
                self.assertGreater(jd, self.jds[i + 1])
                self.assertGreaterEqual(jd, self.jds[i + 1])

    def testAddTimedelta(self):
        delta = datetime.timedelta(days=1, hours=9)
        jd = self.jds[0] + delta
        answer = JulianDate(2000, 1, 2, 21, 0, 0)
        self.assertEqual(jd, answer)

        # test __radd__
        jd = delta + self.jds[0]
        self.assertEqual(jd, answer)

    def testSubtract(self):
        answer = answers[1] - answers[0]
        diff = self.jds[1] - self.jds[0]
        self.assertEqual(diff, answer)

        # Negative difference
        answer = answers[4] - answers[5]
        diff = self.jds[4] - self.jds[5]
        self.assertEqual(diff, answer)

    def testCopyingAndHashing(self):
        cpy = copy(self.jds[0])
        self.assertEqual(cpy, self.jds[0])
        self.assertIsNot(cpy, self.jds[0])

        self.assertEqual(hash(cpy), hash(self.jds[0]))

    def testProperties(self):
        for jd, number in zip(self.jds, answers):
            with self.subTest(jd=jd, number=number, msg='value'):
                self.assertEqual(jd.value, number)
            with self.subTest(jd=jd, number=number, msg='number'):
                self.assertEqual(jd.number, int(number))
            with self.subTest(jd=jd, number=number, msg='fraction'):
                self.assertAlmostEqual(jd.fraction, number % 1)

        self.assertEqual(self.jds[0].timezone, 0)
        jd = JulianDate(2023, 10, 6, 0, 0, 0, -5)
        self.assertEqual(jd.timezone, -5)
        jd = JulianDate(2023, 10, 6, 0, 0, 0, -5.5)
        self.assertEqual(jd.timezone, -5.5)

    def testFuture(self):
        jd = JulianDate(2000, 1, 1, 12, 0, 0)
        future = jd.future(1.5)
        answer = JulianDate(2000, 1, 3, 0, 0, 0)
        self.assertEqual(future, answer)

        future = jd.future(-1.5)
        answer = JulianDate(1999, 12, 31, 0, 0, 0)
        self.assertEqual(future, answer)

    def testStringMethods(self):
        # Most of these should be tested from the __str__ method, but need to test the timezone and n (rounding) args.
        # Testing timezones.
        jd = JulianDate(2000, 1, 1, 12, 0, 0)
        self.assertEqual(jd.date(), '2000/01/01 12:00:00.0 +0 UTC')
        self.assertEqual(jd.date(-5), '2000/01/01 07:00:00.0 -5 UTC')
        self.assertEqual(jd.date(5), '2000/01/01 17:00:00.0 +5 UTC')

        self.assertEqual(jd.day(), '2000/01/01')
        self.assertEqual(jd.day(-13), '1999/12/31')
        self.assertEqual(jd.day(13), '2000/01/02')

        self.assertEqual(jd.time(), '12:00:00.0')
        self.assertEqual(jd.time(-5), '07:00:00.0')
        self.assertEqual(jd.time(5), '17:00:00.0')

        # Testing rounding.
        jd = JulianDate(2000, 1, 1, 12, 0, 12.3456789)
        self.assertEqual(jd.date(n=5), '2000/01/01 12:00:12.34568 +0 UTC')
        self.assertEqual(jd.date(n=2), '2000/01/01 12:00:12.35 +0 UTC')
        self.assertEqual(jd.date(n=0), '2000/01/01 12:00:12.0 +0 UTC')
        self.assertEqual(jd.date(n=None), '2000/01/01 12:00:12 +0 UTC')

        self.assertEqual(jd.time(n=5), '12:00:12.34568')
        self.assertEqual(jd.time(n=2), '12:00:12.35')
        self.assertEqual(jd.time(n=0), '12:00:12.0')
        self.assertEqual(jd.time(n=None), '12:00:12')

    def testDatetimeConversion(self):
        datetimeList = [datetime.datetime(2000, 1, 1, 12, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1999, 1, 1, 0, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1988, 6, 19, 12, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1988, 1, 27, 0, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1987, 6, 19, 12, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1987, 1, 27, 0, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1900, 1, 1, 0, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1600, 12, 31, 0, 0, 0, 0, datetime.timezone.utc),
                        datetime.datetime(1600, 1, 1, 0, 0, 0, 0, datetime.timezone.utc),
                        # This one just rounds differently, no need to focus on making it perfect.
                        datetime.datetime(837, 4, 10, 7, 11, 59, 999999, datetime.timezone.utc)]

        for jd, dtime in zip(self.jds[:10], datetimeList):
            with self.subTest(jd=jd):
                self.assertEqual(jd.toDatetime(), dtime)

    def testDayOfYear(self):
        jd = JulianDate(2023, 9, 24, 12, 0, 0, 0)
        self.assertEqual(jd.dayOfYear(), 267)
        jd = JulianDate(2024, 9, 24, 12, 0, 0, 0)
        self.assertEqual(jd.dayOfYear(), 268)

    @staticmethod
    def _datetimeSideEffect(*args, **kwargs):
        return datetime.datetime(*args, **kwargs)

    def testNow(self):
        # todo: this doesn't test the timezone mechanism inside now()
        # With tm_gmtoff as None in time.localtime()
        args = (2023, 10, 26, 23, 50, 29)
        tmp = datetime.datetime(*args, 0)
        with mock.patch('sattrack.core.juliandate.datetime.datetime') as nowMock:
            nowMock.now.return_value = tmp
            nowMock.side_effect = self._datetimeSideEffect
            self.assertEqual(now(), JulianDate(*args[:6], 0))

        # With timezone set
        args = (2023, 10, 26, 23, 50, 29)
        tz = datetime.timezone(datetime.timedelta(hours=-5))
        tmp = datetime.datetime(*args, 0, tz)
        with mock.patch('sattrack.core.juliandate.datetime.datetime') as nowMock:
            nowMock.now.return_value = tmp
            nowMock.side_effect = self._datetimeSideEffect
            self.assertEqual(now(), JulianDate(*args[:6], -5))

    def testJ2000(self):
        self.assertEqual(J2000.number, 2451545)
        self.assertEqual(J2000.fraction, 0.0)
        self.assertEqual(J2000.timezone, 0.0)


if __name__ == '__main__':
    unittest.main()
