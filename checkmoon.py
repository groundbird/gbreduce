import astropy as ap
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
# from astropy import Time

telescope = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)
time = ap.time.Time('2019-11-12T01:00:00')
print(time)
print(ap.coordinates.get_moon(time,location=telescope))
