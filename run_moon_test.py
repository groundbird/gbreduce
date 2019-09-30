import matplotlib.pyplot as plt
import numpy as np
import math

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.modeling import models, fitting
import scipy.fftpack
from scipy import signal, optimize
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from astropy.time import Time
from astropy.coordinates import Angle
import astropy_speedups
from scipy.optimize import curve_fit
import os
import datetime
import time
import ephem
import multiprocessing as mp

from scipy.signal import find_peaks

coordsys = 1
nside = 512
npix = 12*nside*nside

def combine_sky_maps(skymaps,hitmaps,prefix,outputname,centralpos=(0,0),plotlimit=0.0):

	skymap = np.zeros(npix, dtype=np.float)
	hitmap = np.zeros(npix, dtype=np.float)

	nummaps = len(skymaps)
	for i in range(0,nummaps):
		try:
			print(skymaps[i])
			inputmap = hp.read_map(skymaps[i])
			inputhitmap = hp.read_map(hitmaps[i])
		except:
			continue
		for j in range(0,npix):
			if inputhitmap[j] > 0:
				skymap[j] = skymap[j] + inputmap[j]*inputhitmap[j]
				hitmap[j] = hitmap[j] + inputhitmap[j]

	# We now have a combined map, time to normalise it
	for i in range(0,npix):
		if hitmap[i] >= 1:
			skymap[i] = skymap[i]/hitmap[i]
		else:
			skymap[i] = hp.pixelfunc.UNSEEN

	write_healpix_map(skymap,prefix,outputname+'_skymap.fits')
	write_healpix_map(hitmap,prefix,outputname+'_hitmap.fits')

	hp.mollview(skymap)
	plt.savefig(outputname+'_skymap.png')
	plt.close()
	plt.clf()
	hp.mollview(hitmap)
	plt.savefig(outputname+'_hitmap.png')
	if plotlimit != 0.0:
		hp.gnomview(skymap,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
	else:
		hp.gnomview(skymap,rot=centralpos,reso=5)
	plt.savefig(outputname+'_zoom.png')
	plt.close()
	plt.clf()

	return



def calc_healpix_pixels(skypos):
	if coordsys == 0:
		healpix_pixel = hp.ang2pix(nside, (np.pi/2)-Angle(skypos.b).radian, Angle(skypos.l).radian)
		pos = (np.median(Angle(skypos.l).degree),np.median((Angle(skypos.b).degree)))
	else:
		print(skypos)
		healpix_pixel = hp.ang2pix(nside, (np.pi/2)-Angle(skypos.dec).radian, Angle(skypos.ra).radian)
		pos = (np.median(Angle(skypos.ra).degree),np.median((Angle(skypos.dec).degree)))
	return healpix_pixel, pos

def write_healpix_map(data, prefix, outputname,headerinfo=[]):
	extra_header = []
	extra_header.append(("instrum",("tfgi")))
	# extra_header.append(("tgfi_v",(str(self.ver))))
	extra_header.append(("obsname",(prefix)))
	now = datetime.datetime.now()
	extra_header.append(("writedat",(now.isoformat())))
	# for args in extra_header:
	# 	print(args[0].upper())
	# 	print(args[1:])
	hp.write_map(outputname,data,overwrite=True,extra_header=extra_header)
	return 


indir = "/Users/mpeel/Desktop/GB/"
indir = "/Users/mpeel/Documents/git/groundbird/rhea_comm/"
outdir = indir+'plots/'
datafile = 'moon2.dat'
onsky = [0, 1, 4, 5, 8, 9, 12, 13]

# maps = []
# hitmaps = []
# for val in onsky:
# 	maps.append(outdir+'skymap_'+str(val+1)+'.fits')
# 	hitmaps.append(outdir+'hitmap_'+str(val+1)+'.fits')
# combine_sky_maps(maps, hitmaps, 'moontest', outdir+'combined', centralpos=[323.43,20.73])
# exit()

data = np.loadtxt(indir+datafile)
datarate = 1000
rpm = 2.045
deg_per_sec = rpm*360.0/60.0
numsamples = len(data[:])
numrotations = rpm*numsamples/(60.0*datarate)
print(numsamples/(60.0*datarate))
print(numrotations)
print(rpm)
print(deg_per_sec)
print(deg_per_sec*60.0)
az = np.arange(numsamples)*deg_per_sec/datarate
az = az % 360
el = 70.0

timearr = Time(2458747.5+np.arange(numsamples)/(24.0*60.0*60.0*datarate), format='jd')
telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)


position = AltAz(az=az*u.deg,alt=el*np.ones(numsamples)*u.deg,location=telescope,obstime=timearr)
skypos = position.transform_to(ICRS)
healpix_pixel, centralpos = calc_healpix_pixels(skypos)


plt.plot(az)
plt.savefig(outdir+'test.png')
plt.clf()

# for chan in range(2):
for chan in range(16):
	offset = np.median(data[:,chan+1])
	data[:,chan+1] = data[:,chan+1]-offset
	scale = np.max(data[:,chan+1])
	scale2 = np.min(data[:,chan+1])
	if np.abs(scale) > np.abs(scale2):
		data[:,chan+1] = data[:,chan+1]/scale
	else:
		data[:,chan+1] = data[:,chan+1]/scale2


	try:
		peaks, _ = find_peaks(data[:,chan+1], height=0.6)

		plt.plot(data[:,0]/datarate,data[:,chan+1])
		plt.plot(data[peaks,0]/datarate, data[peaks,chan+1], "x")
		plt.savefig(outdir+'ch'+str(chan+1)+'.png')
		plt.clf()
		plt.plot(az,data[:,chan+1])
		plt.savefig(outdir+'ch'+str(chan+1)+'_vs_az.png')
		plt.clf()

		plt.plot(data[range(peaks[0]-100,peaks[0]+400),0]/datarate,data[range(peaks[0]-100,peaks[0]+400),chan+1])
		plt.savefig(outdir+'ch'+str(chan+1)+'_zoom.png')
		plt.clf()
		
		print(peaks)
		print(az[peaks])
		print(data[peaks,0])
		print(data[peaks,chan+1])
	except:
		print('huh')


	skymap = np.zeros(npix, dtype=np.float)
	hitmap = np.zeros(npix, dtype=np.float)
	for i in range(0,len(healpix_pixel)):
		skymap[healpix_pixel[i]] = skymap[healpix_pixel[i]] + data[i,chan+1]
		hitmap[healpix_pixel[i]] = hitmap[healpix_pixel[i]] + 1
	for i in range(0,len(skymap)):
		if hitmap[i] >= 1:
			skymap[i] = skymap[i]/hitmap[i]
		else:
			skymap[i] = hp.pixelfunc.UNSEEN

	write_healpix_map(skymap,'moontest',outdir+'/skymap_'+str(chan+1)+'.fits')
	write_healpix_map(hitmap,'moontest',outdir+'/hitmap_'+str(chan+1)+'.fits')
	hp.mollview(skymap,title='Channel ' + str(chan+1))
	plt.savefig(outdir+'/skymap_'+str(chan+1)+'.png')
	plt.clf()
	hp.mollview(hitmap,title='Channel ' + str(chan+1))
	plt.savefig(outdir+'/hitmap_'+str(chan+1)+'.png')
	plt.clf()

	hp.gnomview(skymap,rot=[323.43,20.73],reso=1.5,title='Channel ' + str(chan+1))
	plt.savefig(outdir+'/skymap_'+str(chan+1)+'zoom.png')
	plt.clf()

	# if chan in onsky:
	# 	plt.plot(data[:,0]/datarate,np.abs(data[:,chan+1])-np.abs(data[:,chan+3]))
	# 	plt.savefig(outdir+'ch'+str(chan+1)+'-ch'+str(chan+3)+'.png')
	# 	plt.clf()

