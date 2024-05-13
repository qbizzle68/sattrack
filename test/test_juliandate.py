import datetime
import pickle
import unittest

from sattrack.core.juliandate import J2000, JulianDate, DateComponents

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


class TestJulianDateNew(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls) -> None:
    #     jds = []
    #     for args in values:

    def testInit(self):
        for comps, answer in zip(values, answers):
            with self.subTest(comps=comps, number=answer):
                jd = JulianDate(*comps)
                self.assertEqual(jd.value, answer)

    def testClassInit(self):
        # Test the date May 9th, 2024, 23:38:00 -0500 UTC.
        answer = JulianDate(2024, 5, 9, 23, 38, 0, -5)
        jd = JulianDate.fromNumber(2460440, 0.6930555555555564, -5)
        self.assertEqual(jd, answer, 'testing 2024-05-09 23:38:00 -0500 from JulianDate.fromNumber')

        timezone = datetime.timezone(datetime.timedelta(hours=-5))
        _datetime = datetime.datetime(2024, 5, 9, 23, 38, 0, tzinfo=timezone)
        jd = JulianDate.fromDatetime(_datetime)
        self.assertEqual(jd, answer, 'testing 2024-05-09 23:38:00 -0500 from JulianDate.fromDatetime')

        jd = JulianDate.fromisoformat('2024-05-09 23:38:00 -0500')
        self.assertEqual(jd, answer, 'testing 2024-05-09 23:38:00 -0500 from JulianDate.fromisoformat')

        # this screws up other type checking
        # with mock.patch('sattrack.core.juliandate.time') as structMock:
        #     struct = time.struct_time((2024, 5, 9, 23, 38, 0, 0, 0, 1, None, -5*3600))
        #     structMock.localtime.return_value = struct
        #     with mock.patch('sattrack.core.juliandate.datetime') as nowMock:
        #         nowMock.datetime.now.return_value = _datetime
        #         jd = JulianDate.now()

        # Since we're completely dependent on the datetime module here, just test the simple case
        jd = JulianDate.strptime('2024-05-09 23:38:00 -0500', '%Y-%m-%d %H:%M:%S %z')
        self.assertEqual(jd, answer, 'testing 2024-05-09 23:38:00 -0500 from JulianDate.strptime')

    def testClassNegativeInit(self):
        # Test the date February 29th, -1000 12:34:56 -0200 UTC and March 1st, -1001 12:34:56 -0200 UTC.
        answer1 = JulianDate(-1000, 2, 29, 12, 34, 56, -2)
        answer2 = JulianDate(-1001, 3, 1, 12, 34, 56, -2)

        jd = JulianDate.fromNumber(1355867, 0.1075925925, -2)
        self.assertEqual(jd, answer1, 'testing -1000-02-29 12:34:56 -0200 from JulianDate.fromNumber')
        jd = JulianDate.fromNumber(1355502, 0.1075925925, -2)
        self.assertEqual(jd, answer2, 'testing -1001-03-01 12:34:56 -0200 from JulianDate.fromNumber')

        # We can't handle negative years with this right now.
        # jd = JulianDate.fromisoformat('-1000-02-29 12:34:56 -0200')
        # self.assertEqual(jd, answer1, 'testing -1000-02-29 12:34:56 -0200 from JulianDate.fromisoformat')
        # jd = JulianDate.fromisoformat('-1001-03-01 12:34:56 -0200')
        # self.assertEqual(jd, answer2, 'testing -1000-02-29 12:34:56 -0200 from JulianDate.fromisoformat')

        # We can't handle negative years with this right now.
        # jd = JulianDate.strptime('-1000-02-29 12:34:56 -0200', '%Y-%m-%d %H:%M:%S %z')
        # self.assertEqual(jd, answer1, 'testing -1000-02-29 12:34:56 -0200 from JulianDate.strptime')
        # jd = JulianDate.strptime('-1001-03-01 12:34:56 -0200', '%Y-%m-%d %H:%M:%S %z')
        # self.assertEqual(jd, answer2, 'testing -1000-02-29 12:34:56 -0200 from JulianDate.strptime')

    def testAsTimezone(self):
        jd = JulianDate(2024, 5, 10, 4, 38, 0, 0)
        answer = JulianDate(2024, 5, 9, 23, 38, 0, -5)

        jd = jd.asTimezone(-5)
        self.assertEqual(jd.value, answer.value, 'adjusting a timezone with asTimezone()')

    def testStrings(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.8888888, -5)

        self.assertEqual(str(jd), '2460440.6930658435 --- 2024-05-09 23:38:00.888889 -0500 UTC',
                         'testing __str__')
        self.assertEqual(repr(jd), 'JulianDate(2024, 5, 9, 23, 38, 0.8888888, '
                                   'datetime.timezone(datetime.timedelta(days=-1, seconds=68400)))',
                         'testing __repr__')
        self.assertEqual(jd.date(), '2024-05-09 23:38:00.889 -0500 UTC', 'testing .date()')
        self.assertEqual(jd.date(0), '2024-05-10 04:38:00.889 +0000 UTC', 'testing .date() with timezone=0')
        self.assertEqual(jd.date(n=6), '2024-05-09 23:38:00.888889 -0500 UTC', 'testing .date() with n=6')
        self.assertEqual(jd.date(n=0), '2024-05-09 23:38:01 -0500 UTC', 'testing .date() with n=0')

        self.assertEqual(jd.day(), '2024-05-09', 'testing .day()')
        self.assertEqual(jd.time(), '23:38:00.889', 'testing .time()')

    def testFormat(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.8888888, -5)

        specs = ('a', 'A', 'w', 'd', 'b', 'B', 'm', 'y', 'Y', 'H', 'I', 'p', 'M', 'S', 'f', 'z', 'Z',
                 'j', 'U', 'W', 'c', 'x', 'X', 'G', 'u', 'V')
        results = ('Thu', 'Thursday', '4', '09', 'May', 'May', '05', '24', '2024', '23', '11', 'PM',
                   '38', '00', '888889', '-0500', 'UTC-05:00', '130', '18', '19',
                   'Thu May  9 23:38:00 2024', '05/09/24', '23:38:00', '2024', '4', '19')
        for spec, result in zip(specs, results):
            with self.subTest('testing __format__', spec=spec):
                self.assertEqual(format(jd, f'%{spec}'), result, f'testing __format__ with spec {spec}')

        for spec, result in zip(specs, results):
            with self.subTest('testing strftime', spec=spec):
                self.assertEqual(jd.strftime(f'%{spec}'), result, f'testing .strftime with spec {spec}')

        self.assertEqual(format(jd, 'day is %d'), 'day is 09', 'testing year formatting with text')

        with self.assertRaises(ValueError, msg='testing invalid format spec exception'):
            format(jd, '%P')

    def testRound(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.8888888, -5)

        self.assertAlmostEqual(round(jd, 3).value, 2460440.693065845, 5, 'testing __round__ with n=3')
        self.assertAlmostEqual(round(jd).value, 2460440.6930671297, 5, 'testing __round__ with n=None')

    def testRoundFormat(self):
        jd = JulianDate(2024, 5, 9, 23, 59, 59.9999)

        self.assertEqual(jd.time(n=4), '23:59:59.9999', 'testing not rounding up with n=4')
        self.assertEqual(jd.time(n=2), '00:00:00.00', 'testing rounding up with n=2')
        self.assertEqual(jd.date(n=2), '2024-05-10 00:00:00.00 +0000 UTC', 'testing rounding up with n=2')

        self.assertEqual(jd.time(n=0), '00:00:00', 'testing rounding up with n=0')
        self.assertEqual(jd.time(n=6), '23:59:59.999900', 'testing not rounding up with n=6')
        with self.assertRaises(ValueError, msg='testing negative n exception'):
            jd.date(n=-1)

    def testProperties(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)

        self.assertEqual(type(jd.components), DateComponents, 'testing components property is correct type')
        self.assertEqual(jd.number, 2460440, 'testing number property equality')
        self.assertAlmostEqual(jd.fraction, 0.6930569844444445, msg='testing fraction property equality')
        self.assertAlmostEqual(jd.value, 2460440.6930569843, msg='testing value property equality')
        self.assertEqual(jd.timezone, datetime.timezone(datetime.timedelta(hours=-5)),
                         'testing timezone property equality')
        self.assertEqual(jd.utcOffset, -5, 'testing utcOffset property equality')

    def testSerialization(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)

        self.assertEqual(jd.toDict(), {'number': 2460440, 'fraction': 0.6930569844444445, 'utcOffset': -5.0},
                         'testing toDict equality')
        self.assertEqual(jd.toJson(), '{"number": 2460440, "fraction": 0.6930569844444445, "utcOffset": -5.0}',
                         'testing toJson equality')
        reduceResults = (JulianDate, (*jd.components, jd.timezone))
        self.assertEqual(jd.__reduce__(), reduceResults, 'testing __reduce__ return values')
        pickled = pickle.dumps(jd)
        pickleAnswer = b'\x80\x04\x95\x80\x00\x00\x00\x00\x00\x00\x00\x8c\x18sattrack.core.juliandate' \
                       b'\x94\x8c\nJulianDate\x94\x93\x94(M\xe8\x07K\x05K\tK\x17K&G?\xbf\x9a\xcf\xfa~\xb6' \
                       b'\xbf\x8c\x08datetime\x94\x8c\x08timezone\x94\x93\x94h\x03\x8c\ttimedelta\x94\x93' \
                       b'\x94J\xff\xff\xff\xffJ0\x0b\x01\x00K\x00\x87\x94R\x94\x85\x94R\x94t\x94R\x94.'
        self.assertEqual(pickled, pickleAnswer, 'testing pickle.dumps works with jd')

    def testMath(self):
        jd1 = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)
        jd2 = JulianDate(2024, 5, 11, 5, 38, 0.123456, -5)

        self.assertAlmostEqual(jd2-jd1, 1.25, msg='testing JulianDate subtraction')
        self.assertEqual(jd1 + 1.25, jd2, 'testing JulianDate float addition')
        self.assertEqual(1.25 + jd1, jd2, 'testing JulianDate __radd__')
        self.assertEqual(jd2 - 1.25, jd1, 'testing JulianDate float subtraction')

        delta = datetime.timedelta(hours=5, minutes=3, seconds=2.1)
        additionAnswer = JulianDate(2024, 5, 10, 4, 41, 2.223456, -5)
        self.assertEqual(jd1 + delta, additionAnswer, 'testing JulianDate datetime.timedelta addition')
        self.assertEqual(delta + jd1, additionAnswer, 'testing JulianDate datetime.timedelta __radd__')

    def testComparison(self):
        jd1 = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)
        jd2 = JulianDate(2024, 5, 11, 5, 38, 0.123456, -5)
        jd3 = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)

        self.assertEqual(jd1, jd3, 'testing JulianDate equality')
        self.assertNotEqual(jd1, jd2, 'testing JulianDate equality is false')
        self.assertLess(jd1, jd2, 'testing JulianDate less than')
        self.assertLessEqual(jd1, jd3, 'testing JulianDate less than or equal')
        self.assertGreater(jd2, jd1, 'testing JulianDate greater than')
        self.assertGreaterEqual(jd1, jd3, 'testing JulianDate greater than or equal')

        self.assertEqual(hash(jd1), -5442201358416014943, 'testing JulianDate hash')

    def testMiscellaneous(self):
        jd = JulianDate(2024, 5, 9, 23, 38, 0.123456, -5)

        timezone = datetime.timezone(datetime.timedelta(hours=-5))
        _datetime = datetime.datetime(2024, 5, 9, 23, 38, 0, 123456, timezone)
        self.assertEqual(jd.toDatetime(), _datetime, 'testing .toDatetime() method')
        jdNegativeYear = JulianDate(-1, 1, 1, 1, 1, 1)
        with self.assertRaises(ValueError, msg='testing negative year .toDatetime exception'):
            jdNegativeYear.toDatetime()

        self.assertEqual(jd.dayOfYear(), 130, 'testing .dayOfYear() method')

    def testJ2000(self):
        self.assertEqual(J2000.value, 2451545.0)

    def testNegativeYearFormatting(self):
        jd = JulianDate(-1000, 2, 2, 2, 2, 2, -5)

        self.assertTrue(jd._formatter.proxy, 'testing formatter proxy property value')

        self.assertEqual(format(jd, '%y'), '-00', 'testing negative year formatting %y')
        self.assertEqual(format(jd, '%Y'), '-1000', 'testing negative year formatting %Y')
        # this test requires some brute force since locale.nl_langinfo is not always defined
        _datetime = jd._formatter.datetime
        sampleFormat = _datetime.strftime('%c').replace('6688', '-1000').replace('88', '00')
        self.assertEqual(format(jd, '%c'), sampleFormat, 'testing negative year formatting %c')
        self.assertEqual(format(jd, '%z'), '-0500', 'testing negative year formatting %z')
        _datetime = jd._formatter.datetime
        sampleFormat = _datetime.strftime('%x').replace('6688', '-1000').replace('88', '-00')
        self.assertEqual(format(jd, '%x'), sampleFormat, 'testing negative year formatting %x')
        self.assertEqual(format(jd, '%G'), '-1000', 'testing negative year formatting %G')
        self.assertEqual(format(JulianDate(0, 1, 1, 0, 0, 0), '%G'), '0000', 'testing year zero formatting %G')
        self.assertEqual(format(jd, '%d'), '02', 'testing negative year other formatting')

    def testTimezones(self):
        import zoneinfo

        timezone = zoneinfo.ZoneInfo('US/Central')
        jd = JulianDate(2024, 5, 9, 23, 38, 0.123456, timezone)
        self.assertEqual(str(jd), '2460440.6930569843 --- 2024-05-09 23:38:00.123456 -0500 UTC',
                         'testing zoneinfo timezone string')
        self.assertEqual(repr(jd), "JulianDate(2024, 5, 9, 23, 38, 0.123456, zoneinfo.ZoneInfo(key='US/Central'))",
                         'testing zoneinfo timezone repr')
        self.assertEqual(jd.utcOffset, -5, 'testing timezone in daylight savings time')
        jd = JulianDate(2024, 1, 1, 0, 0, 0, timezone)
        self.assertEqual(jd.utcOffset, -6, 'testing timezone in standard time')


if __name__ == '__main__':
    unittest.main()
