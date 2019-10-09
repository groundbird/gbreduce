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
from scipy import optimize
import os
from astroutils import *
import datetime
import time
import emcee
import corner

from gbreduce_functions import *
from gbreduce_read import *
from rhea_comm.lib_read_rhea import *
import mkid_pylibs as klib

class gbreduce:
	def __init__(self,datadir="/net/nas/proyectos/quijote2/tod/",outdir='',nside=512,datarate=1000,use_mkidpylibs=True):
		# Directories
		self.datadir = datadir
		self.outdir = outdir

		# Basic telescope information
		self.telescope = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)
		self.datarate = datarate

		# Analysis options
		self.use_mkidpylibs = use_mkidpylibs

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
		now = datetime.now()
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

	def analyse_swp(self, prefix, filename, datarate=0,freqrange=3.0,freqstep=0.01,centerfreqs=[], lo=4.99e9):
		print(freqrange)
		print(freqstep)
		print(centerfreqs)
		# Sort out the output directories
		if prefix[-1] != '/':
			prefix = prefix + "/"
		plotext = 'plots/'
		print(self.outdir+prefix)
		ensure_dir(self.outdir+prefix)
		ensure_dir(self.outdir+prefix+plotext)

		# Don't do this for now, just use the standard library.
		if not self.use_mkidpylibs:
			# Use the datarate from the file or the default datarate if not otherwise set.
			if datarate == 0:
				# try:
				# 	datarate = dataset[0]['rate']
				# except:
				datarate = self.datarate

			freqset, dataset = read_rhea_swp_data(filename)
			print(freqset)
			print(len(freqset))
			print(dataset[:,0])
			print(len(dataset[:,1]))
			# Approximate rescale
			dataset /= ((2**28) * 200.e6 / datarate)

			# Remove values that are close to zero
			# print(len(dataset[(np.abs(dataset[:,1])+np.abs(dataset[:,2]))<1e-5,1]))
			# dataset = dataset[(np.abs(dataset[:,1])+np.abs(dataset[:,2]))>1e-5,:]

			print(dataset)
			samplesize = np.shape(dataset)
			samplefreqs = np.arange(samplesize[0])*freqstep
			# print(samplesize)
			# print(samplefreqs)
			# print(np.shape(dataset))
			# print(len(dataset[np.abs(dataset[:,2])<1e-5,2]))
			# print(dataset[0:100,1])
			# print(centerfreqs[0])
			samplefreqs = (samplefreqs+centerfreqs[0]-1.5) * 1e6
			# print(samplefreqs)
			samplefreqs = samplefreqs + lo
			# print(samplefreqs)

		# Get the first fit from mkids_lib
		if self.use_mkidpylibs:
			kid = klib.kidana.KidAnalyzer()
			paramset = []
			for i in range(len(centerfreqs)):
				kid.swpfromfile('rhea',filename,lo,i)
				kid.swp.fitIQ(nfwhm=3, fitter='gaolinbg')
				# kid.swp.fitresult.report()
				mkid_fit_params = kid.swp.fitresult.values()
				print(mkid_fit_params)
				# print(mkid_fit_params['c'])
				# params = np.array([mkid_fit_params['absa'], mkid_fit_params['arga'], mkid_fit_params['tau']*1e7, mkid_fit_params['c']*1e8, mkid_fit_params['fr']*1e-9, mkid_fit_params['Qr']*1e-4, mkid_fit_params['Qc']*1e-4, mkid_fit_params['phi0']])
				# print(params)
				paramset.append(mkid_fit_params)
			# exit()
		else:
			# params = np.array([10.0, 1.0, 1.0, 10.0, (centerfreqs[0]+1.5)*1e6+4e9, 1.0, 1.0, 1.0])

			# samplefreqs = (samplefreqs) * 1e6
			# params = np.array([10.0, 1.0, 1.0, 10.0, (1.5)*1e6, 1.0, 1.0, 1.0])

			# params = np.array([10.0, 1.0, 1.0, 10.0, 1.5, 1.0, 1.0, 1.0])
			# params = np.array([0.12, 305, 6.5e-7, -2.2e-8, 4.9074e9, 25232, 85000, -0.28])
			params = np.array([0.12244757113931891, 305.97531201109734, 6.516879931537661, -2.2256369833132223, 4907366610.603888, 25232.341646851175, 85014.87075854921, -0.2812287453378864])
			#	a, w, t, c, fr, qr, qc, p0 = param
			#=absa, arga, tau, c, fr, Qr, Qc, phi0 in Honda's code

		# Focus on the cut-out
		# dataset = dataset[1300:1800,:]
		# samplefreqs = samplefreqs[1300:1800]

		if not self.use_mkidpylibs:
			err_estimate = np.ones(len(dataset[:,0]))*(np.std(dataset[1:30,0]) + 1j*np.std(dataset[1:30,1]))
			# print(err_estimate)

			plt.plot(dataset[:,0],dataset[:,1])
			testfit = fit_mkids(samplefreqs,params)
			plt.plot(testfit.real,testfit.imag)
			plt.savefig(self.outdir+prefix+plotext+'plot_before.png')
			plt.clf()

			print(len(params))
			print(len(samplefreqs))
			# # tofit = dataset[:,1]+dataset[:,2]*1j
			# # print(len(tofit))
			param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_mkids, params, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j)),full_output=True)
			print(params)
			print(param_est)
			print(mesg_result)
			plt.plot(dataset[:,0],dataset[:,1])
			testfit = fit_mkids(samplefreqs,param_est)
			plt.plot(testfit.real,testfit.imag)
			plt.savefig(self.outdir+prefix+plotext+'plotfit.png')
			plt.clf()

			nll = lambda *args: -lnlike(*args)

			ndim = len(params)
			bounds = np.zeros((ndim,2))
			for i in range(0,ndim):
				# bounds[i][0] = -10.0*params[i]
				# bounds[i][1] =  10.0*params[i]
				# if bounds[i][0] > bounds[i][1]:
				# 	bounds[i][0] =  10.0*params[i]
				# 	bounds[i][1] = -10.0*params[i]
				bounds[i][0] = None
				bounds[i][1] = None
				if i == 4: # frequency
					bounds[i][0] = params[i]-0.00005
					bounds[i][1] = params[i]+0.00005
					# bounds[i][0] = min(samplefreqs)*1e-9
					# bounds[i][1] = max(samplefreqs)*1e-9
			print(bounds)

			result = optimize.minimize(nll, param_est, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(err_estimate)),method='L-BFGS-B',bounds=bounds)#, options={'gtol': 1e-30, 'disp': False, 'maxiter': 4000})		result = optimize.minimize(nll, params, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(0.1*dataset[:,0]+0.1*dataset[:,1]*1j)),method='L-BFGS-B')#, options={'gtol': 1e-30, 'disp': False, 'maxiter': 4000})
			maxlikelihood = result["x"]
			print("Done")
			print(result)
			print(maxlikelihood)
			print(result['success'])
			print(result['message'])

			plt.plot(dataset[:,0],dataset[:,1])
			testfit = fit_mkids(samplefreqs,maxlikelihood)
			plt.plot(testfit.real,testfit.imag)
			plt.savefig(self.outdir+prefix+plotext+'plotfit_likelihood.png')
			plt.clf()
			# exit()
			# Let's do some MCMC fitting!
			nwalkers = 100

			# not_fixed *= 1e-4 # Use a random distribution only for values that aren't fixed.
			pos = [params + 1e-4*params*np.random.randn(ndim) for i in range(nwalkers)]
			sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=8,args=(bounds, samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(err_estimate)))
			sampler.run_mcmc(pos, 100)

			samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

			plt.plot(samples[:,0])
			plt.savefig(self.outdir+prefix+plotext+'mc_test.png')
			plt.close()
			plt.plot(samples[:,1])
			plt.savefig(self.outdir+prefix+plotext+'mc_test1.png')
			plt.close()
			plt.plot(samples[:,2])
			plt.savefig(self.outdir+prefix+plotext+'mc_test2.png')
			plt.close()
			plt.plot(samples[:,3])
			plt.savefig(self.outdir+prefix+plotext+'mc_test3.png')
			plt.close()
			fig = corner.corner(samples)
			fig.savefig(self.outdir+prefix+plotext+"mc_triangle.png")
			plt.close()

			p1,p2,p3,p4,p5,p6,p7,p8 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
	                             zip(*np.percentile(samples, [16, 50, 84],
	                                                axis=0)))
			print(p1)
			print(p2)
			print(p3)
			print(p4)
			print(p5)
			print(p6)
			print(p7)
			print(p8)
			mcparams = [p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0]]
			plt.plot(dataset[:,0],dataset[:,1])
			testfit = fit_mkids(samplefreqs,mcparams)
			plt.plot(testfit.real,testfit.imag)
			plt.savefig(self.outdir+prefix+plotext+'plotfit_mc.png')
			plt.clf()
			print(params)

			for i in range(0,int(samplesize[1]/2)):
				plot_val_tod(dataset[:,i*2],dataset[:,i*2+1],self.outdir+prefix+plotext+'IQ_'+str(i)+'.png',formatstr='b')

				# plot_val_tod(dataset[:,1]-dataset[:,3],dataset[:,2]-dataset[:,4],self.outdir+prefix+plotext+'test2.png')
				plot_tod(abs(dataset[:,i*2]),self.outdir+prefix+plotext+'I_'+str(i)+'.png',formatstr='b')
				plot_tod(abs(dataset[:,i*2+1]),self.outdir+prefix+plotext+'Q_'+str(i)+'.png',formatstr='b')
				plot_tod(abs(dataset[:,i*2]+dataset[:,i*2+1]*1j),self.outdir+prefix+plotext+'amp_'+str(i)+'.png',formatstr='b')
				plot_tod(angle(dataset[:,i*2]+dataset[:,i*2+1]*1j),self.outdir+prefix+plotext+'phase_'+str(i)+'.png',formatstr='b')

		return paramset

	def analyse_tod(self, prefix, filename, el=70,rpm=2,datarate=0,starttime=2458747.5,numpix=4, plotlimit=0.0, quiet=False, dofft=False, plottods=True, plotmap=True, dopol=False, plotcombination=True, numfiles=50,swp_params=[],centerfreqs=[],lo=4.99e9):

		# Sort out the output directories
		if prefix[-1] != '/':
			prefix = prefix + "/"
		plotext = 'plots/'
		print(self.outdir+prefix)
		ensure_dir(self.outdir+prefix)
		ensure_dir(self.outdir+prefix+plotext)

		print(centerfreqs)
		freqs = np.asarray(centerfreqs) * 1.0e6
		freqs = freqs +lo
		print(freqs)

		# Read in the data
		# rawdata = read_rhea_tod(filename)
		# if swpname != '':
		# 	swp = read_rhea_swp(swpname)

		# numchannels = len(rawdata)
		# print(np.shape(rawdata[0]['phase']))
		# dataset = np.asarray([rawdata[0]['phase']])
		# print(np.shape(dataset))
		# for i in range(1,numchannels):
		# 	dataset = np.append(dataset, [rawdata[i]['phase'][:]], axis=0)
		# print(np.shape(dataset))
		# dataset = np.transpose(dataset)
		# print(np.shape(dataset))
		# dataset = [rawdata[0]['phase'], rawdata[1]['phase']]

		# Use the datarate from the file or the default datarate if not otherwise set.
		if datarate == 0:
			# try:
			# 	datarate = dataset[0]['rate']
			# except:
			datarate = self.datarate

		dataset = read_rhea_data(filename)
		# dataset = read_rhea_data(swpname)
		print(dataset)
		# Approximate rescale
		dataset[:,1:] /= ((2**28) * 200.e6 / datarate)
		# Trim the first and last entries, as they can be bad.
		dataset = dataset[1:-1]
		print(dataset)

		plot_val_tod(dataset[:,1],dataset[:,2],self.outdir+prefix+plotext+'test.png')
		plot_val_tod(dataset[:,1]-dataset[:,4],dataset[:,2]-dataset[:,4],self.outdir+prefix+plotext+'test2.png')
		plot_tod(abs(dataset[:,1]+dataset[:,2]*1j),self.outdir+prefix+plotext+'test3.png')
		plot_tod(angle(dataset[:,1]+dataset[:,2]*1j),self.outdir+prefix+plotext+'test4.png')

		test = dataset[:,1] + dataset[:,2]*1j
		plot_val_tod(test.real, test.imag, self.outdir+prefix+plotext+'test5.png')
		# exit()

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
		# print(jd)
		# exit()
		skypos = self.calc_positions(az, el, jd)
		healpix_pixel, centralpos = self.calc_healpix_pixels(skypos)

		plot_tod(az,self.outdir+prefix+plotext+'/plot_az.png')
		plot_tod(skypos.ra,self.outdir+prefix+plotext+'/plot_ra.png')
		plot_tod(skypos.dec,self.outdir+prefix+plotext+'/plot_dec.png')
		# plot_tod(pa,self.outdir+prefix+plotext+'/plot_pa.png')

		for pix in range(0,numpix*2):
			rewound = gao_rewind(freqs[pix],dataset[:,pix*2+1]+1j*dataset[:,pix*2+2], swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
			# rewound_offpeak = gao_rewind(dataset[:,pix*4]+3, dataset[:,pix*4]+4, swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
			dataset[:,pix*2+1] = abs(rewound)
			dataset[:,pix*2+2] = angle(rewound)
			# dataset[:,pix*4+2] = abs(rewound_offpeak)
			# dataset[:,pix*4+3] = angle(rewound_offpeak)

		for pix in range(0,numpix*4):

			plot_tod(dataset[:,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'.png')

			# Correct for any offset
			offset = np.median(dataset[:,pix+1])
			print(offset)
			dataset[:,pix+1] -= offset

			# If we're looking at phase, unwind if needed
			dataset[dataset[:,pix+1]>np.pi,pix+1] -= 2.0*np.pi

			# If we're highly negative, then we probably need to invert.
			if np.max(dataset[:,pix+1]) < -np.min(dataset[:,pix+1])/2.0:
				dataset[:,pix+1] = -dataset[:,pix+1]

			# Also rescale to 1.0 for now
			# scale = np.max(dataset[:,pix+1])
			# scale2 = np.min(dataset[:,pix+1])
			# if np.abs(scale) < np.abs(scale2):
			# 	scale = scale2
			# dataset[:,pix+1] /= scale

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

		# # Write out the TOD to disk
		# ra_col = fits.Column(name='ra',format='E',array=np.array(Angle(skypos.ra).degree))
		# dec_col = fits.Column(name='dec',format='E',array=np.array(Angle(skypos.dec).degree))
		# jd_col = fits.Column(name='mjd',format='E',array=jd-2400000.5)
		# col1 = fits.Column(name='ch1',format='E',array=dataset[:,1])
		# col2 = fits.Column(name='ch2',format='E',array=dataset[:,2])
		# col3 = fits.Column(name='ch3',format='E',array=dataset[:,3])
		# col4 = fits.Column(name='ch4',format='E',array=dataset[:,4])
		# col5 = fits.Column(name='ch5',format='E',array=dataset[:,5])
		# col6 = fits.Column(name='ch6',format='E',array=dataset[:,6])
		# col7 = fits.Column(name='ch7',format='E',array=dataset[:,7])
		# col8 = fits.Column(name='ch8',format='E',array=dataset[:,8])
		# col9 = fits.Column(name='ch9',format='E',array=dataset[:,9])
		# col10 = fits.Column(name='ch10',format='E',array=dataset[:,10])
		# col11 = fits.Column(name='ch11',format='E',array=dataset[:,11])
		# col12 = fits.Column(name='ch12',format='E',array=dataset[:,12])
		# col13 = fits.Column(name='ch13',format='E',array=dataset[:,13])
		# col14 = fits.Column(name='ch14',format='E',array=dataset[:,14])
		# col15 = fits.Column(name='ch15',format='E',array=dataset[:,15])
		# col16 = fits.Column(name='ch16',format='E',array=dataset[:,16])
		# cols = fits.ColDefs([ra_col, dec_col, jd_col, col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16])
		# hdu = fits.BinTableHDU.from_columns(cols)

		# hdr = fits.Header()
		# hdr['INSTRUM'] = 'GroundBIRD'
		# hdr['VERSION'] = str(self.ver)
		# hdr['OBSNAME'] = prefix
		# primary_hdu = fits.PrimaryHDU(header=hdr)
		# hdul = fits.HDUList([primary_hdu, hdu])
		# hdul.writeto(self.outdir+prefix+prefix[:-1]+'_tod.fits',overwrite=True)
		return