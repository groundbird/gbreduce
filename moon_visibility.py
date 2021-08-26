import astropy as ap
import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon, get_body
import matplotlib.pyplot as plt
import datetime
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)


tenerife = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)

starttime = [Time('2020-07-01 0:00:00'),Time('2020-08-01 0:00:00'),Time('2020-09-01 0:00:00'),Time('2020-10-01 0:00:00'),Time('2020-11-01 0:00:00'),Time('2020-12-01 0:00:00'),Time('2021-01-01 0:00:00'),Time('2021-02-01 0:00:00'),Time('2021-03-01 0:00:00')]
endtime = [Time('2020-08-01 00:00:00'),Time('2020-09-01 0:00:00'),Time('2020-10-01 0:00:00'),Time('2020-11-01 0:00:00'),Time('2020-12-01 0:00:00'),Time('2021-01-01 0:00:00'),Time('2021-02-01 0:00:00'),Time('2021-03-01 0:00:00'),Time('2021-04-01 0:00:00')]
sampling = 5000
formatter = DateFormatter('%d')
for k in range(0,len(starttime)):
	print('\n')
	print(starttime[k])
	fig, ax = plt.subplots(figsize=(12, 6))
	time_delta = (endtime[k].mjd - starttime[k].mjd)/sampling
	times_mjd = np.arange(starttime[k].mjd,endtime[k].mjd,time_delta)
	times = Time(times_mjd,format='mjd')

	sun = get_sun(times)
	altaz = sun.transform_to(AltAz(obstime=times,location=tenerife))
	plt.plot_date(times.plot_date,altaz.alt.deg,label='Sun',fmt='-')
	moon = get_moon(times,tenerife)
	altaz = moon.transform_to(AltAz(obstime=times,location=tenerife))
	plt.plot_date(times.plot_date,altaz.alt.deg,label='Moon',fmt='-')
	venus = get_body('venus',times,tenerife)
	altaz = venus.transform_to(AltAz(obstime=times,location=tenerife))
	plt.plot_date(times.plot_date,altaz.alt.deg,label='Venus',fmt='-')

	plt.legend()
	plt.ylim([60,90])
	plt.xlabel('Time (UTC)')
	plt.ylabel('Elevation (deg)')
	plt.title(starttime[k].to_value('iso', subfmt='date'))
	ax.xaxis.set_major_formatter(formatter)
	plt.savefig('moon/visibility_'+str(starttime[k])+'.png')
	plt.clf()
