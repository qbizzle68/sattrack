from sattrack.satpass import ShadowController
from sattrack.spacetime.juliandate import now
from sattrack.structures.tle import getTle

starlink_numbers = [
    1007, 1011, 1022, 1029,
    1073, 1071, 1110, 1085,
    1132, 1141, 1173, 1138,
    1201, 1234, 1279, 1301,
    1306, 1294, 1346, 1332,
    1441, 1393, 1461, 1475,
    1391, 1523, 1524, 1515,
    1585, 1588, 1654, 1686,
    1550, 1644, 1687, 1774
]

starlink_tles = {}
for n in starlink_numbers:
    print(f'fetching tle for starlink {n}')
    tle = getTle(f'starlink-{n}')
    starlink_tles[n] = tle
    print('received')


if __name__ == '__main__':
    jd = now()
    print('jd:', jd.date(-5))
    for k, v in starlink_tles.items():
        print('starlink', k)
        sc = ShadowController(v)
        sc.computeValues(jd)
        times = sc.getTimes()
        phis = sc.getAnomalies()
        print('entry time:', times[0].date(-5))
        print('exit time:', times[1].date(-5))
        print('entry anomaly:', phis[0])
        print('exit anomaly:', phis[1])
        if times[0].value() < times[1].value():
            print(f'starlink {k} ------------------------- PASSED')
        else:
            print(f'starlink {k} ------------------------- FAILED')
