from sattrack.api import *

tle = getTle('zarya')
iss = Satellite(tle)

jd = now()
geo = GeoPosition(38.0608, -97.9298)
np = getNextPass(iss, geo, jd)
plist = getPassList(iss, geo, jd, 2)
