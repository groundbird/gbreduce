#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Make a map of the sky showing the time during the day to observe it with GB
# 
# Version history:
#
# 11-Nov-2019  M. Peel       Started

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from astropy.time import Time
from astropy.coordinates import Angle
import astropy_speedups
from scipy import optimize
import os
from astroutils import *
from gbreduce import gbreduce
import datetime
import time

now = datetime.datetime.now()
el = 70.0

run = gbreduce(nside=512)
prefix=''

times = []
az = []
for i in range(0,24):
	for j in range(0,60):
		for k in np.arange(0,360,0.1):
			times.append(str(now.year)+'-'+str('{0:02d}'.format(now.month))+'-'+str('{0:02d}'.format(now.day))+'T'+str('{0:02d}'.format(i))+':'+str('{0:02d}'.format(j))+':00.0')
			az.append(k)
print(now.day)
# print(k)
# print(times)
print(len(times))
# exit()
t = Time(times, format='isot', scale='utc')

skypos = run.calc_positions(az, el, t.jd)
healpix_pixel, centralpos = run.calc_healpix_pixels(skypos)

skymap = np.zeros(hp.nside2npix(512), dtype=np.float)
for i in range(0,len(healpix_pixel)):
		skymap[healpix_pixel[i]] = times[i][11:13]

run.write_healpix_map(skymap,prefix,'timemap.fits')

hp.mollview(skymap,title='timemap')
plt.savefig('timemap.png')
plt.clf()
