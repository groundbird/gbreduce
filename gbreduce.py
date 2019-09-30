#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# The main class for gbreduce
# 
# Version history:
#
# 30-Sep-2019  M. Peel       Started

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
from scipy.optimize import curve_fit
import os
from astroutils import *
import datetime
import time

from gbreduce_functions import *
from gbreduce_read import *

class gbreduce:
	def __init__(self,datadir="/net/nas/proyectos/quijote2/tod/",outdir='',nside=512,datarate=1000):
		# Directories
		self.datadir = datadir
		self.outdir = outdir

		# Basic telescope information
		self.telescope = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)
		self.datarate = datarate

		# Map-making information
		self.nside = nside
		self.npix = hp.nside2npix(self.nside)
		self.coordsys = 1 # 0 = Galactic, 1 = ICRS.

		# Constants
		self.k = 1.380e-23 # Boltzman constant m2 kg s-2 K-1
		self.c = 2.9979e8

		# Stable version control
		self.ver = 0.0

	def calc_JytoK(self,beam,freq):
		return 1e-26*((self.c/freq)**2)/(2.0*self.k*calc_beam_area(beam))

	def calc_farfield(self,diameter,frequency=0.0,wavelength=0.0):
		if wavelength != 0.0:
			return 2.0 * diameter**2 / wavelength
		else:
			return 2.0 * diameter**2 * frequency / self.c


	def calc_positions(self, az, el, jd):

		timearr = Time(jd, format='jd')
		# Apply the pointing model
		# print('Applying pointing model')
		# az,el = self.pointing_model(az[0],el[0])

		position = AltAz(az=az*u.deg,alt=el*u.deg,location=self.telescope,obstime=timearr)
		if self.coordsys == 0:
			skypos = position.transform_to(Galactic)
			pa = []
		else:
			skypos = position.transform_to(ICRS)
			# print('Calculating position angles')
			# pa = self.apobserver.parallactic_angle(timearr,skypos)
			# print('Done with positions')
		return skypos#,pa

	def write_fits_file(self, dataset,prefix,outputname):
		hdu = fits.PrimaryHDU(dataset)
		hdul = fits.HDUList([hdu])
		hdul.writeto(outputname)

	def write_healpix_map(self, data, prefix, outputname,headerinfo=[]):
		extra_header = []
		extra_header.append(("instrum",("tfgi")))
		extra_header.append(("tgfi_v",(str(self.ver))))
		extra_header.append(("obsname",(prefix)))
		now = datetime.datetime.now()
		extra_header.append(("writedat",(now.isoformat())))
		# for args in extra_header:
		# 	print(args[0].upper())
		# 	print(args[1:])
		hp.write_map(outputname,data,overwrite=True,extra_header=extra_header)
		return 

	def calc_healpix_pixels(self, skypos):
		if self.coordsys == 0:
			healpix_pixel = hp.ang2pix(self.nside, (np.pi/2)-Angle(skypos.b).radian, Angle(skypos.l).radian)
			pos = (np.median(Angle(skypos.l).degree),np.median((Angle(skypos.b).degree)))
		else:
			print(skypos)
			healpix_pixel = hp.ang2pix(self.nside, (np.pi/2)-Angle(skypos.dec).radian, Angle(skypos.ra).radian)
			pos = (np.median(Angle(skypos.ra).degree),np.median((Angle(skypos.dec).degree)))
		return healpix_pixel, pos

	def analyse_tod(self, prefix, filename, el=70,rpm=2,datarate=0,starttime=2458747.5,pixelrange=range(0,16), plotlimit=0.0, quiet=False, dofft=False, plottods=True, plotmap=True, dopol=False, plotcombination=True, numfiles=50):

		# Sort out the output directories
		if prefix[-1] != '/':
			prefix = prefix + "/"
		plotext = 'plots/'
		print(self.outdir+prefix)
		ensure_dir(self.outdir+prefix)
		ensure_dir(self.outdir+prefix+plotext)

		# Read in the data
		dataset = np.asarray(read_rhea_data(filename),dtype=np.float64)
		# Trim the first and last entries, as they can be bad.
		dataset = dataset[1:-1]

		# Use the default datarate if not otherwise set.
		if datarate == 0:
			datarate = self.datarate

		# Make up the az+el information
		deg_per_sec = rpm*360.0/60.0
		numsamples = len(dataset[:])
		numrotations = rpm*numsamples/(60.0*datarate)
		print(numsamples/(60.0*datarate))
		print(numrotations)
		print(rpm)
		print(deg_per_sec)
		print(deg_per_sec*60.0)
		az = np.arange(numsamples)*deg_per_sec/datarate
		az = az % 360

		jd = starttime+np.arange(numsamples)/(24.0*60.0*60.0*datarate)

		skypos = self.calc_positions(az, el, jd)
		healpix_pixel, centralpos = self.calc_healpix_pixels(skypos)

		plot_tod(az,self.outdir+prefix+plotext+'/plot_az.png')
		plot_tod(skypos.ra,self.outdir+prefix+plotext+'/plot_ra.png')
		plot_tod(skypos.dec,self.outdir+prefix+plotext+'/plot_dec.png')
		# plot_tod(pa,self.outdir+prefix+plotext+'/plot_pa.png')

		for pix in pixelrange:

			plot_tod(dataset[:,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'.png')

			# Correct for any offset
			offset = np.median(dataset[:,pix+1])
			print(offset)
			dataset[:,pix+1] -= offset

			# Also rescale to 1.0 for now
			scale = np.max(dataset[:,pix+1])
			scale2 = np.min(dataset[:,pix+1])
			if np.abs(scale) < np.abs(scale2):
				scale = scale2
			dataset[:,pix+1] /= scale

			plot_tod(dataset[:,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'_rescale.png')

			skymap = np.zeros(self.npix, dtype=np.float)
			hitmap = np.zeros(self.npix, dtype=np.float)
			for i in range(0,len(healpix_pixel)):
				skymap[healpix_pixel[i]] = skymap[healpix_pixel[i]] + dataset[i,pix+1]
				hitmap[healpix_pixel[i]] = hitmap[healpix_pixel[i]] + 1
			for i in range(0,len(skymap)):
				if hitmap[i] >= 1:
					skymap[i] = skymap[i]/hitmap[i]
				else:
					skymap[i] = hp.pixelfunc.UNSEEN

			self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_'+str(pix+1)+'.fits')
			self.write_healpix_map(hitmap,prefix,self.outdir+prefix+'/hitmap_'+str(pix+1)+'.fits')

			hp.mollview(skymap,title='Channel ' + str(pix+1))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix+1)+'.png')
			plt.clf()
			hp.mollview(hitmap,title='Channel ' + str(pix+1))
			plt.savefig(self.outdir+prefix+'/hitmap_'+str(pix+1)+'.png')
			plt.clf()

			hp.gnomview(skymap,rot=centralpos,reso=20,title='Channel ' + str(pix+1))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix+1)+'view.png')
			plt.clf()

			hp.gnomview(skymap,rot=[323.43,20.73],reso=1.5,title='Channel ' + str(pix+1))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix+1)+'zoom.png')
			plt.clf()

		# Write out the TOD to disk
		ra_col = fits.Column(name='ra',format='E',array=np.array(Angle(skypos.ra).degree))
		dec_col = fits.Column(name='dec',format='E',array=np.array(Angle(skypos.dec).degree))
		jd_col = fits.Column(name='jd',format='E',array=jd)
		col1 = fits.Column(name='ch1',format='E',array=dataset[:,1])
		col2 = fits.Column(name='ch2',format='E',array=dataset[:,2])
		col3 = fits.Column(name='ch3',format='E',array=dataset[:,3])
		col4 = fits.Column(name='ch4',format='E',array=dataset[:,4])
		col5 = fits.Column(name='ch5',format='E',array=dataset[:,5])
		col6 = fits.Column(name='ch6',format='E',array=dataset[:,6])
		col7 = fits.Column(name='ch7',format='E',array=dataset[:,7])
		col8 = fits.Column(name='ch8',format='E',array=dataset[:,8])
		col9 = fits.Column(name='ch9',format='E',array=dataset[:,9])
		col10 = fits.Column(name='ch10',format='E',array=dataset[:,10])
		col11 = fits.Column(name='ch11',format='E',array=dataset[:,11])
		col12 = fits.Column(name='ch12',format='E',array=dataset[:,12])
		col13 = fits.Column(name='ch13',format='E',array=dataset[:,13])
		col14 = fits.Column(name='ch14',format='E',array=dataset[:,14])
		col15 = fits.Column(name='ch15',format='E',array=dataset[:,15])
		col16 = fits.Column(name='ch16',format='E',array=dataset[:,16])
		cols = fits.ColDefs([ra_col, dec_col, jd_col, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16])
		hdu = fits.BinTableHDU.from_columns(cols)

		hdr = fits.Header()
		hdr['INSTRUM'] = 'GroundBIRD'
		hdr['VERSION'] = str(self.ver)
		hdr['OBSNAME'] = prefix
		primary_hdu = fits.PrimaryHDU(header=hdr)
		hdul = fits.HDUList([primary_hdu, hdu])
		hdul.writeto(self.outdir+prefix+prefix[:-1]+'_tod.fits',overwrite=True)
		return
