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
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5, get_moon, get_sun
from astropy.time import Time
from astropy.coordinates import Angle
from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator
# import astropy_speedups
from scipy import optimize, signal
import os
from astroutils import *
import datetime
import time
import emcee
import corner
import gc
import timeit

from gbreduce_functions import *
from gbreduce_read import *
from rhea_comm.lib_read_rhea import *
import mkid_pylibs as klib
from matplotlib.gridspec import GridSpec
from merge_az import merge_az

class gbreduce:
	def __init__(self,datadir="/net/nas/proyectos/quijote2/tod/",outdir='',azdir='',eldir='',pixinfo='',elcompressed=False,domedir='',tempdir='',nside=512,datarate=1000,use_mkidpylibs=True,focallength=500.0, boresight=96.0, astrom_interp_rate = 300 * u.s):
		# Directories
		self.datadir = datadir
		self.outdir = outdir
		self.azdir = azdir
		self.eldir = eldir
		self.domedir = domedir
		self.tempdir = tempdir
		self.pixinfo = pixinfo

		# Basic telescope information
		self.telescope = EarthLocation(lat=28.300467*u.deg, lon=-16.510288*u.deg, height=2390*u.m)
		self.datarate = datarate
		self.focallength = focallength
		self.boresight = boresight

		# Analysis options
		self.use_mkidpylibs = use_mkidpylibs
		self.elcompressed = elcompressed

		# Map-making information
		self.nside = nside
		self.npix = hp.nside2npix(self.nside)
		self.coordsys = 1 # 0 = Galactic, 1 = ICRS.

		# Constants
		self.k = 1.380e-23 # Boltzman constant m2 kg s-2 K-1
		self.c = 2.9979e8

		self.astrom_interp_rate = astrom_interp_rate

		# Stable version control
		self.ver = 0.0

		# Preload some information
		if self.pixinfo != '':
			self.pixel_az, self.pixel_el, self.pixel_freq = get_pixinfo(self.pixinfo, self.focallength, self.boresight)


	def calc_JytoK(self,beam,freq):
		return 1e-26*((self.c/freq)**2)/(2.0*self.k*calc_beam_area(beam))

	def calc_farfield(self,diameter,frequency=0.0,wavelength=0.0):
		if wavelength != 0.0:
			return 2.0 * diameter**2 / wavelength
		else:
			return 2.0 * diameter**2 * frequency / self.c

	def get_second_brightest(self, timestream):
		maxpix = np.argmax(np.abs(timestream))
		timestream[maxpix] = np.median(timestream)
		return np.argmax(np.abs(timestream))

	def calc_positions(self, az, el, jd):

		timearr = Time(jd, format='jd')
		# Apply the pointing model
		# print('Applying pointing model')
		# az,el = self.pointing_model(az[0],el[0])

		# Check to see if the elevation is under 90°, and set it to 90° if not
		el[el>90.0]=90.0

		with erfa_astrom.set(ErfaAstromInterpolator(self.astrom_interp_rate)):
			position = AltAz(az=az*u.deg,alt=el*u.deg,location=self.telescope,obstime=timearr)
			# position = AltAz(az=u.Quantity(az, u.deg, copy=False),alt=u.Quantity(el, u.deg, copy=False),location=self.telescope,obstime=timearr)
			if self.coordsys == 0:
				skypos = position.transform_to(Galactic())
				pa = []
			else:
				skypos = position.transform_to(ICRS())

				# print('Calculating position angles')
				# pa = self.apobserver.parallactic_angle(timearr,skypos)
				# print('Done with positions')
		return skypos#,pa

	def write_fits_file(self, dataset,prefix,outputname):
		hdu = fits.PrimaryHDU(dataset)
		hdul = fits.HDUList([hdu])
		hdul.writeto(outputname)

	def write_fits(self, ra, dec, jd, A, P, hdr, outfile,moondist=[],sundist=[],flag=[]):
		colset = [fits.Column(name='ra',format='E',array=np.array(ra)), fits.Column(name='dec',format='E',array=np.array(dec)), fits.Column(name='mjd',format='D',array=jd), fits.Column(name='A',format='E',array=A), fits.Column(name='P',format='E',array=P)]
		if moondist != []:
			colset.append(fits.Column(name='moondist',format='E',array=moondist))
		if sundist != []:
			colset.append(fits.Column(name='sundist',format='E',array=sundist))
		if flag != []:
			colset.append(fits.Column(name='flag',format='I',array=flag))
		cols = fits.ColDefs(colset)
		hdu = fits.BinTableHDU.from_columns(cols)
		primary_hdu = fits.PrimaryHDU(header=hdr)
		hdul = fits.HDUList([primary_hdu, hdu])
		hdul.writeto(outfile,overwrite=True)
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

	def write_healpix_map(self, data, prefix, outputname,headerinfo=[]):
		extra_header = []
		extra_header.append(("instrum",("gb")))
		extra_header.append(("gbreduc",(str(self.ver))))
		extra_header.append(("obsname",(prefix)))
		now = datetime.now()
		extra_header.append(("writedat",(now.isoformat())))
		# for args in extra_header:
		# 	print(args[0].upper())
		# 	print(args[1:])
		hp.write_map(outputname,data,overwrite=True,extra_header=extra_header)
		return

	def get_skymap(self, healpix_pixel,dataset, flag=[], rethitmap=True):
		skymap = np.zeros(self.npix, dtype=np.float)
		hitmap = np.zeros(self.npix, dtype=np.float)
		if len(flag) == 0:
			flag = np.ones(len(healpix_pixel))
		for i in range(0,len(healpix_pixel)):
			if flag[i] == 1:
				skymap[healpix_pixel[i]] = skymap[healpix_pixel[i]] + dataset[i]
				hitmap[healpix_pixel[i]] = hitmap[healpix_pixel[i]] + 1
		for i in range(0,len(skymap)):
			if hitmap[i] >= 1:
				skymap[i] = skymap[i]/hitmap[i]
			else:
				skymap[i] = hp.pixelfunc.UNSEEN
		if rethitmap:
			return skymap, hitmap
		else:
			return skymap

	def combine_sky_maps(self,skymaps,hitmaps,prefix,outputname,centralpos=(0,0),plotlimit=0.0):

		skymap = np.zeros(self.npix, dtype=np.float)
		hitmap = np.zeros(self.npix, dtype=np.float)

		nummaps = len(skymaps)
		for i in range(0,nummaps):
			try:
				inputmap = hp.read_map(skymaps[i])
				inputhitmap = hp.read_map(hitmaps[i])
			except:
				continue
			for j in range(0,self.npix):
				if inputhitmap[j] > 0:
					skymap[j] = skymap[j] + inputmap[j]*inputhitmap[j]
					hitmap[j] = hitmap[j] + inputhitmap[j]

		# We now have a combined map, time to normalise it
		for i in range(0,self.npix):
			if hitmap[i] >= 1:
				skymap[i] = skymap[i]/hitmap[i]
			else:
				skymap[i] = hp.pixelfunc.UNSEEN
				hitmap[i] = hp.pixelfunc.UNSEEN

		self.write_healpix_map(skymap,prefix,outputname+'_skymap.fits')
		self.write_healpix_map(hitmap,prefix,outputname+'_hitmap.fits')

		hp.mollview(skymap)
		plt.savefig(outputname+'_skymap.png')
		plt.close()
		plt.clf()

		hp.mollview(skymap,norm='hist')
		plt.savefig(outputname+'_skymap_hist.png')
		plt.close()
		plt.clf()

		skymap[skymap != hp.pixelfunc.UNSEEN] -= np.median(skymap[skymap != hp.pixelfunc.UNSEEN])
		std = np.std(skymap[skymap != hp.pixelfunc.UNSEEN])
		hp.mollview(skymap,max=2.0*std,min=-2.0*std)
		plt.savefig(outputname+'_skymap_std.png')
		plt.close()
		plt.clf()

		hp.mollview(hitmap,xsize=1400,norm='log')
		plt.savefig(outputname+'_hitmap.png')
		hp.mollview(hitmap,xsize=1400,max=np.median(hitmap[hitmap != hp.pixelfunc.UNSEEN]))
		plt.savefig(outputname+'_hitmap_median.png')
		if plotlimit != 0.0:
			hp.gnomview(skymap,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
		else:
			hp.gnomview(skymap,rot=centralpos,reso=5)
		plt.savefig(outputname+'_zoom.png')
		plt.close()
		plt.clf()


		hitmap_128 = hp.ud_grade(hitmap, nside_out=128,power=-2)
		hp.mollview(hitmap_128,xsize=1400,norm='log')
		plt.savefig(outputname+'_hitmap_128.png')
		hp.mollview(hitmap,xsize=1400,max=np.median(hitmap_128[hitmap_128 != hp.pixelfunc.UNSEEN]))
		plt.savefig(outputname+'_hitmap_128_median.png')
		plt.clf()

		return


	def runset(self, subdir='',ext='',skipfirst=0,doswp=True):
		folderlist = os.listdir(self.datadir+subdir)
		todolist = []
		for folder in folderlist:
			if '.' not in folder:
				subfolderlist = os.listdir(self.datadir+subdir+folder)
				search = 'KSPS'
				test = [f for f in subfolderlist if search in f]
				print(test)
				if test != []:
					print('No data')
					continue
				test = [f for f in subfolderlist if 'tod' in f]
				if test == []:
					print('No data')
					continue
				subfolderlist = [f for f in subfolderlist if search not in f]
				if len(subfolderlist) > 2:
					todolist.append(subdir+folder)
		print(todolist)

		badrun = []
		if skipfirst == 0:
			trip = 1
		else:
			trip = 0
		for item in todolist:
			if trip < skipfirst:
				print('Skipping ' + item)
				trip += 1
				continue
			print('Running ' + item)
			folder = item+'/'
			starttime = datetime(int(folder[8:12]),int(folder[12:14]),int(folder[14:16]),int(folder[22:24]),int(folder[24:26]),int(folder[26:28]),tzinfo=pytz.timezone('UTC')).timestamp()
			name=folder[8:-1].replace('/','_')+ext

			date = datetime.utcfromtimestamp(starttime)
			contents = os.listdir(self.datadir+folder)
			todo = []
			for line in contents:
				if 'kids' in line and 'list' in line:
					todo.append(line.replace('kids','').replace('list',''))
			for line in todo:
				print(line)
				if 'GB01' in line:
					continue
				swpfile = 'swp'+line+'rawdata'
				todfile = 'tod'+line+'rawdata'
				kidparams = 'kids'+line+'list'
				# try:
				if doswp:
					swp_params = self.analyse_swp(name+line.replace('.',''),self.datadir+folder+swpfile,kidparams=self.datadir+folder+kidparams)
				else:
					swp_params = []
				tod_analysis = self.analyse_tod(name+line.replace('.',''),self.datadir+folder+todfile,kidparams=self.datadir+folder+kidparams,swp_params=swp_params,starttime=starttime)
				# except:
					# badrun.append(name+line)
		if len(badrun) > 0:
			print('There were some failed runs: ' + str(badrun))
		return


	def analyse_swp(self, prefix, filename, datarate=0,freqrange=0.0,freqstep=0.0,centerfreqs=[], lo=4.99e9,starttime=0.0,kidparams=[]):
		print(prefix)
		print(filename)
		if kidparams == []:
			# We have an old style file, read in the parameters from the filename
			fileinfo = parse_filename(filename)
			# [fileformat, freqs, samplespeed, date]
			if fileinfo[0] != 'swp':
				print('Wrong file format!')
				return 0
			else:
				if freqrange == 0.0:
					freqrange = fileinfo[2]
				if freqstep == 0.0:
					freqstep = fileinfo[3]
				if centerfreqs == []:
					centerfreqs = fileinfo[1]
				if starttime == 0.0:
					starttime = fileinfo[4]
		else:
			kids = read_kidslist(kidparams)
			print(kids)
			centerfreqs = kids['kids_freqs']
			lo = kids['sg_freq']

			# TODO Need to modify this to get the frequency info from somewhere other than the filename!

			print(freqrange)
			print(freqstep)
			print(centerfreqs)
			print(starttime)
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
			print('Dataset shape:')
			print(np.shape(dataset))
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
			samplefreqs = (samplefreqs+centerfreqs[0]-1.5)
			# print(samplefreqs)
			samplefreqs = samplefreqs + lo
			# print(samplefreqs)

		# Get the first fit from mkids_lib
		if self.use_mkidpylibs == True:
			# kid = klib.kidana.KidAnalyzer()
			paramset = []
			for i in range(len(centerfreqs)*2):
				try:
					swpset = klib.readfile_swp('rhea',filename,i,lo)

					# Special case for the RIKEN KID
					if centerfreqs[i] == -31900000:
						swpset.fitIQ(nfwhm=-1, frqrange=[4.732e9,4.734e9], fitter='gaolinbg')
					else:
						swpset.fitIQ(nfwhm=3, fitter='gaolinbg')
					mkid_fit_params = swpset.fitresult.values()

					fig,ax = klib.plotter.plotSwp(swpset,kind='rawdata',title='rawdata',color='blue',ls=':',lw=1,marker='.')
					klib.plotter.plotSwp(swpset,kind='fitdata.rawdata',ax=ax,color='red',marker='',lw=2,ls='-',)
					fig.savefig(self.outdir+prefix+plotext+'raw_kidfit_'+str(i)+'.png')
					fig.clf()
					fit.close()
					fig,ax = klib.plotter.plotSwp(swpset,kind='rwdata',title='rwdata',color='blue',ls=':',lw=1,marker='.')
					klib.plotter.plotSwp(swpset,kind='fitdata.rwdata',ax=ax,color='red',marker='',lw=2,ls='-',)
					fig.savefig(self.outdir+prefix+plotext+'rw_kidfit_'+str(i)+'.png')
					fig.clf()
					fit.close()
					del swpset
					gc.collect()

					print(mkid_fit_params)
					# print(mkid_fit_params['c'])
					# params = np.array([mkid_fit_params['absa'], mkid_fit_params['arga'], mkid_fit_params['tau']*1e7, mkid_fit_params['c']*1e8, mkid_fit_params['fr']*1e-9, mkid_fit_params['Qr']*1e-4, mkid_fit_params['Qc']*1e-4, mkid_fit_params['phi0']])
					# print(params)
					paramset.append(mkid_fit_params)
				except:
					# Sweep fitting has failed for some reason, just use the one for the previous channel for now, if there was one
					if len(paramset) >= 1:
						paramset.append(paramset[-1])
					else:
						paramset.append([])
			# exit()
		# else:
		# 	# dataset = dataset[1300:1800,:]
		# 	# samplefreqs = samplefreqs[1300:1800]
		# 	# print(samplefreqs)
		# 	# print(dataset)
		# 	print(len(samplefreqs))
		# 	# First do a fit to a skewed Lorentizian.
		# 	params = Fit_SkewedLorentizian(samplefreqs[1300:1800], dataset[1300:1800,0]**2+dataset[1300:1800,1]**2)
		# 	params1 = params
		# 	print(params)
		# 	params = Fit_7para(samplefreqs, dataset[:,0]+dataset[:,1]*1j, fitresult=params)
		# 	print(params)
		#
		# 	t = dataset[:,0]**2 + dataset[:,1]**2
		# 	fig = plt.figure(figsize=(10, 8))
		# 	gs = GridSpec(8, 2, figure=fig)
		# 	ax1 = fig.add_subplot(gs[0, 0])
		# 	ax1.plot(samplefreqs, dataset[:,0])
		# 	ax2 = fig.add_subplot(gs[0, 1])
		# 	ax2.plot(samplefreqs, dataset[:,1])
		# 	ax3 = fig.add_subplot(gs[1:3, :])
		# 	ax3.plot(dataset[:,0], dataset[:,1])
		# 	# if self.fitresult_sp is not None:
		# 	ax3.plot(Fit_7para_Func(params.params, samplefreqs), Fit_7para_Func(params.params, samplefreqs))
		# 	ax3.set_aspect('equal')
		# 	ax4 = fig.add_subplot(gs[4:6, :])
		# 	ax4.plot(samplefreqs, t)
		# 	ax4.plot(samplefreqs, Fit_SkewedLorentizian_Func(params1.params, samplefreqs))
		# 	ax5 = fig.add_subplot(gs[7:9, :])
		# 	ax5.plot(samplefreqs, Fit_SkewedLorentizian_Func(params1.params, samplefreqs, data=t))
		# 	fig.tight_layout()
		# 	plt.show()
		#
		# 	# samplefreqs = (samplefreqs) * 1e6
		# 	# params = np.array([10.0, 1.0, 1.0, 10.0, (1.5)*1e6, 1.0, 1.0, 1.0])
		#
		# 	# params = np.array([10.0, 1.0, 1.0, 10.0, 1.5, 1.0, 1.0, 1.0])
		# 	# params = np.array([0.12, 305, 6.5e-7, -2.2e-8, 4.9074e9, 25232, 85000, -0.28])
		# 	params = np.array([0.12244757113931891, 305.97531201109734, 6.516879931537661, -2.2256369833132223, 4907366610.603888, 25232.341646851175, 85014.87075854921, -0.2812287453378864])
		# 	#	a, w, t, c, fr, qr, qc, p0 = param
		# 	#=absa, arga, tau, c, fr, Qr, Qc, phi0 in Honda's code

		# Focus on the cut-out
		# dataset = dataset[1300:1800,:]
		# samplefreqs = samplefreqs[1300:1800]

		# if not self.use_mkidpylibs:
		# 	err_estimate = np.ones(len(dataset[:,0]))*(np.std(dataset[1:30,0]) + 1j*np.std(dataset[1:30,1]))
		# 	# print(err_estimate)
		#
		# 	plt.plot(dataset[:,0],dataset[:,1])
		# 	testfit = fit_mkids(samplefreqs,params)
		# 	plt.plot(testfit.real,testfit.imag)
		# 	plt.savefig(self.outdir+prefix+plotext+'plot_before.png')
		# 	plt.clf()
		#
		# 	print(len(params))
		# 	print(len(samplefreqs))
		# 	# # tofit = dataset[:,1]+dataset[:,2]*1j
		# 	# # print(len(tofit))
		# 	param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_mkids, params, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j)),full_output=True)
		# 	print(params)
		# 	print(param_est)
		# 	print(mesg_result)
		# 	plt.plot(dataset[:,0],dataset[:,1])
		# 	testfit = fit_mkids(samplefreqs,param_est)
		# 	plt.plot(testfit.real,testfit.imag)
		# 	plt.savefig(self.outdir+prefix+plotext+'plotfit.png')
		# 	plt.clf()
		#
		# 	nll = lambda *args: -lnlike(*args)
		#
		# 	ndim = len(params)
		# 	bounds = np.zeros((ndim,2))
		# 	for i in range(0,ndim):
		# 		# bounds[i][0] = -10.0*params[i]
		# 		# bounds[i][1] =  10.0*params[i]
		# 		# if bounds[i][0] > bounds[i][1]:
		# 		# 	bounds[i][0] =  10.0*params[i]
		# 		# 	bounds[i][1] = -10.0*params[i]
		# 		bounds[i][0] = None
		# 		bounds[i][1] = None
		# 		if i == 4: # frequency
		# 			bounds[i][0] = params[i]-0.00005
		# 			bounds[i][1] = params[i]+0.00005
		# 			# bounds[i][0] = min(samplefreqs)*1e-9
		# 			# bounds[i][1] = max(samplefreqs)*1e-9
		# 	print(bounds)
		#
		# 	result = optimize.minimize(nll, param_est, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(err_estimate)),method='L-BFGS-B',bounds=bounds)#, options={'gtol': 1e-30, 'disp': False, 'maxiter': 4000})		result = optimize.minimize(nll, params, args=(samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(0.1*dataset[:,0]+0.1*dataset[:,1]*1j)),method='L-BFGS-B')#, options={'gtol': 1e-30, 'disp': False, 'maxiter': 4000})
		# 	maxlikelihood = result["x"]
		# 	print("Done")
		# 	print(result)
		# 	print(maxlikelihood)
		# 	print(result['success'])
		# 	print(result['message'])
		#
		# 	plt.plot(dataset[:,0],dataset[:,1])
		# 	testfit = fit_mkids(samplefreqs,maxlikelihood)
		# 	plt.plot(testfit.real,testfit.imag)
		# 	plt.savefig(self.outdir+prefix+plotext+'plotfit_likelihood.png')
		# 	plt.clf()
		# 	# exit()
		# 	# Let's do some MCMC fitting!
		# 	nwalkers = 100
		#
		# 	# not_fixed *= 1e-4 # Use a random distribution only for values that aren't fixed.
		# 	pos = [params + 1e-4*params*np.random.randn(ndim) for i in range(nwalkers)]
		# 	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=8,args=(bounds, samplefreqs, complex_to_real(dataset[:,0]+dataset[:,1]*1j), complex_to_real(err_estimate)))
		# 	sampler.run_mcmc(pos, 100)
		#
		# 	samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
		#
		# 	plt.plot(samples[:,0])
		# 	plt.savefig(self.outdir+prefix+plotext+'mc_test.png')
		# 	plt.close()
		# 	plt.plot(samples[:,1])
		# 	plt.savefig(self.outdir+prefix+plotext+'mc_test1.png')
		# 	plt.close()
		# 	plt.plot(samples[:,2])
		# 	plt.savefig(self.outdir+prefix+plotext+'mc_test2.png')
		# 	plt.close()
		# 	plt.plot(samples[:,3])
		# 	plt.savefig(self.outdir+prefix+plotext+'mc_test3.png')
		# 	plt.close()
		# 	fig = corner.corner(samples)
		# 	fig.savefig(self.outdir+prefix+plotext+"mc_triangle.png")
		# 	plt.close()
		#
		# 	p1,p2,p3,p4,p5,p6,p7,p8 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
		# 						 zip(*np.percentile(samples, [16, 50, 84],
		# 											axis=0)))
		# 	print(p1)
		# 	print(p2)
		# 	print(p3)
		# 	print(p4)
		# 	print(p5)
		# 	print(p6)
		# 	print(p7)
		# 	print(p8)
		# 	mcparams = [p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0]]
		# 	plt.plot(dataset[:,0],dataset[:,1])
		# 	testfit = fit_mkids(samplefreqs,mcparams)
		# 	plt.plot(testfit.real,testfit.imag)
		# 	plt.savefig(self.outdir+prefix+plotext+'plotfit_mc.png')
		# 	plt.clf()
		# 	print(params)
		#
		# 	for i in range(0,int(samplesize[1]/2)):
		# 		plot_val_tod(dataset[:,i*2],dataset[:,i*2+1],self.outdir+prefix+plotext+'IQ_'+str(i)+'.png',formatstr='b')
		#
		# 		# plot_val_tod(dataset[:,1]-dataset[:,3],dataset[:,2]-dataset[:,4],self.outdir+prefix+plotext+'test2.png')
		# 		plot_tod(abs(dataset[:,i*2]),self.outdir+prefix+plotext+'I_'+str(i)+'.png',formatstr='b')
		# 		plot_tod(abs(dataset[:,i*2+1]),self.outdir+prefix+plotext+'Q_'+str(i)+'.png',formatstr='b')
		# 		plot_tod(abs(dataset[:,i*2]+dataset[:,i*2+1]*1j),self.outdir+prefix+plotext+'amp_'+str(i)+'.png',formatstr='b')
		# 		plot_tod(angle(dataset[:,i*2]+dataset[:,i*2+1]*1j),self.outdir+prefix+plotext+'phase_'+str(i)+'.png',formatstr='b')

		return paramset

	def analyse_tod(self, prefix, filename, el=[],rpm=2,datarate=1000.0,starttime=0.0,numpix=4, plotlimit=0.0, quiet=False, dofft=False, plottods=True, plotmap=True, dopol=False, plotcombination=True, numfiles=50,swp_params=[],centerfreqs=[],lo=4.99e9,kidparams=[]):

		# Sort out the output directories
		if prefix[-1] != '/':
			prefix = prefix + "/"
		plotext = 'plots/'
		print(self.outdir+prefix)
		ensure_dir(self.outdir+prefix)
		ensure_dir(self.outdir+prefix+plotext)

		logfile = open(self.outdir+prefix+"/_log.txt","w")
		logfile.write(str(swp_params))

		on_channels = [0,2,4,6]
		if kidparams == []:
			# We have an old style file, read in the parameters from the filename
			fileinfo = parse_filename(filename)
			# [fileformat, freqs, samplespeed, date]
			if fileinfo[0] != 'tod':
				print('Wrong file format!')
				return 0
			else:
				if datarate == 0.0:
					datarate = fileinfo[2]
				if centerfreqs == []:
					centerfreqs = fileinfo[1]
				if starttime == 0.0:
					starttime = fileinfo[3]
		else:
			kids = read_kidslist(kidparams)
			print(kids)
			logfile.write(str(kids)+'\n')
			centerfreqs = kids['kids_freqs']
			lo = kids['sg_freq']
			on_channels=kids['kids']
			print('Frequencies:')
			print(kids['kids_freqs'])
			print(kids['blinds_freqs'])
			numpix = len(kids['kids_freqs'])
		print('Number of pixels: ' + str(numpix))
		print('Data rate: ' + str(datarate))
		print('Centre frequencies: ' + str(centerfreqs))
		print('Start time: ' + str(starttime))

		freqs = np.asarray(centerfreqs)
		freqs = freqs +lo
		print('Frequencies: ' + str(freqs))
		logfile.write('Frequencies ' + str(freqs) + '\n\n')

		# Use the datarate from the file or the default datarate if not otherwise set.
		if datarate == 0:
			# try:
			# 	datarate = dataset[0]['rate']
			# except:
			datarate = self.datarate

		start = time.time()
		dataset = read_rhea_data_new(filename)
		print(time.time() - start)
		print(dataset[0:20])
		exit()

		start = time.time()
		dataset = read_rhea_data(filename)
		print(time.time() - start)
		print(dataset[0:20])

		startsync = -1
		i=0
		while startsync == -1:
			startsync = dataset[i,-2]
			i+=1
		precount = i
		print(startsync)
		endsync = dataset[-1,-2]
		print(endsync)

		try:
			startsynctime = find_az_synctime(self.azdir, starttime-100, startsync)
		except:
			logfile.write("Bad azimuth data, can't find start time, skipping")
			logfile.close()
			return

		print(startsynctime)
		print(float(precount/datarate))
		starttime = startsynctime - float(precount/datarate)
		print(starttime)

		# Approximate rescale
		dataset[:,1:-2] /= ((2**28) * 200.e6 / datarate)
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

		# Calculate the timing
		numsamples = len(dataset[:])
		timestream = starttime + np.arange(numsamples)/(datarate)
		print(min(timestream))
		print(max(timestream))



		# Either get azimuth data, or make some up...
		aztimes, azdata, azrot, azoffset = fetch_azdata(self.azdir, starttime-200, starttime+len(dataset[:])/datarate+100, compressed=False)
		# print(azdata)
		if(len(azdata) < 1):
			print('No azimuth data!')
			exit()
			# # Make up the az+el information
			# deg_per_sec = rpm*360.0/60.0
			# numrotations = rpm*numsamples/(60.0*datarate)
			# print(numsamples/(60.0*datarate))
			# print(numrotations)
			# print(rpm)
			# print(deg_per_sec)
			# print(deg_per_sec*60.0)
			# az = np.arange(numsamples)*deg_per_sec/datarate
			# az = az % 360
		else:

			# This is the new code
			print(np.shape(dataset))
			try:
				dataset, az = merge_az(dataset[:,0:-2], dataset[:,-2], dataset[:,-1], aztimes, azdata, azrot, azoffset)
				timestream = dataset[:,0]
			except:
				logfile.write('Bad azimuth data, skipping')
				logfile.close()
				return
			print(np.shape(dataset))
			# exit()
			# This is the old code, still in use for now.
			# az = np.zeros(len(timestream))
			# aztimes_pos = 0
			# for i in range(0,len(timestream)):
			# 	while(aztimes[aztimes_pos] < timestream[i] and aztimes_pos < len(aztimes)-1):
			# 		aztimes_pos += 1
			# 	az[i] = azdata[aztimes_pos-1]

		del aztimes
		del azdata
		del azrot
		del azoffset
		gc.collect()

		# # Fetch the dome data
		domet, domef = fetch_domedata(self.domedir, min(timestream), max(timestream))
		flag = np.ones(len(timestream))
		if len(domet) >= 1:
			# We can use the dome opening/closing to do an initial cut of the data
			dometimes_pos = 0
			for i in range(0,len(timestream)):
				while(domet[dometimes_pos] < timestream[i] and dometimes_pos < len(domet)-1):
					dometimes_pos += 1
				if dometimes_pos == len(domet):
					flag[i] = int(domef[dometimes_pos])
				else:
					flag[i] = int(domef[dometimes_pos-1])

		del domet
		del domef
		gc.collect()

		# Fetch the temperature data
		temperaturet, temperaturevals = fetch_tempdata(self.tempdir, min(timestream), max(timestream))
		if len(temperaturet) >= 1:
			# We can use the dome opening/closing to do an initial cut of the data
			temptimes_pos = 0
			for i in range(0,len(timestream)):
				while(temperaturet[temptimes_pos] < timestream[i] and temptimes_pos < len(temperaturet)-1):
					temptimes_pos += 1
				# if temperaturevals[temptimes_pos][1] > 0.28:
				if temperaturevals[temptimes_pos][1] > 0.4:
					flag[i] = 2

		del temperaturet
		del temperaturevals
		gc.collect()

		# If the elevation parameter isn't set, read in the encoder values
		if el == []:
			eltimes, eldata = fetch_eldata(self.eldir, min(timestream)-200, max(timestream)+100, compressed=self.elcompressed)
			el = np.zeros(len(timestream))
			eltimes_pos = 0
			print(len(eltimes))
			for i in range(0,len(timestream)):
				# print(eltimes_pos)
				# print(eltimes)
				while(eltimes[eltimes_pos] < timestream[i] and eltimes_pos < len(eltimes)-1):
					if eltimes_pos >= len(eltimes):
						break
					eltimes_pos += 1

				el[i] = (eldata[eltimes_pos-1])+90.0 # Convert from zenith angle to elevation

		del eltimes
		del eldata
		gc.collect()

		t = Time(timestream, format='unix', scale='utc')

		# Plot the azimuth, elevation, and flags
		plot_tod(az,self.outdir+prefix+plotext+'/plot_az.png')
		azdiff = np.diff(az)
		azdiff[azdiff < -300] += 360.0
		plot_tod(azdiff,self.outdir+prefix+plotext+'/plot_az_diff.png')
		plot_tod(el,self.outdir+prefix+plotext+'/plot_el.png')
		plot_tod(flag,self.outdir+prefix+plotext+'/plot_flag.png')

		del azdiff
		gc.collect()

		if swp_params != []:
			for pix in range(0,numpix):
				pix_I = pix*2+1
				pix_Q = pix*2+2
				newamp = subtract_blindtone_JS(dataset[:,pix_I]+1j*dataset[:,pix_Q], dataset[:,kids['blinds'][0]*2+3]+1j*dataset[:,kids['blinds'][0]*2+4])

				newamp = gao_rewind(freqs[pix],newamp, swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])

				dataset[:,pix_I] = abs(newamp)
				# Need to correct the phase of the angle to avoid negatives
				angletmp = -1.0*angle(newamp)
				angletmp = angletmp - np.sign(angletmp)*np.pi
				# And apply the correction for linearity
				dataset[:,pix_Q] = 2.0*np.tan(angletmp/2.0)
				del newamp
				del angletmp
				gc.collect()

				# rewound = gao_rewind(freqs[pix],dataset[:,pix*4+1]+1j*dataset[:,pix*4+2], swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
				# # rewound_offpeak = gao_rewind(dataset[:,pix*4]+3, dataset[:,pix*4]+4, swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
				# dataset[:,pix*2+1] = abs(rewound)
				# dataset[:,pix*2+2] = angle(rewound)
				# # dataset[:,pix*4+2] = abs(rewound_offpeak)
				# # dataset[:,pix*4+3] = angle(rewound_offpeak)

		# sos = signal.butter(10, 0.333, 'hp', fs=1000, output='sos')

		# Calculating moon positions
		print('Moon and sun calculations...')
		print(len(t))

		numinterp = 10000
		moonpos_az, moonpos_el, moonpos_ra, moonpos_dec = get_body_azel_radec(self.telescope, t, numinterp,body='moon')
		sunpos_az, sunpos_el, sunpos_ra, sunpos_dec = get_body_azel_radec(self.telescope, t, numinterp,body='sun')
		print(len(moonpos_az))
		print(len(moonpos_az))
		# moonpos = get_moon(t,self.telescope)
		# moonpos_azel = moonpos.transform_to(AltAz(location=self.telescope,obstime=t))
		print('Moon and sun calculations done.')
		for pix in range(0,numpix):
			pix_A = pix*2+1
			pix_P = pix*2+2
			print('Pixel number ' + str(pix))
			logfile.write('\n\nPixel '+str(pix)+'\n')
			start = 100
			end = -100

			matchfreq = int(freqs[pix]/1e6)
			match = -1
			print('Match frequency: ' + str(matchfreq))
			try:
				match = self.pixel_freq.index(float(matchfreq))
			except:
				try:
					match = self.pixel_freq.index(float(matchfreq)-1.0)
				except:
					try:
						match = self.pixel_freq.index(float(matchfreq)+1.0)
					except:
						match = -1

			logfile.write('Match number is ' + str(match)+'\n')
			if match >= 0:
				logfile.write('Using az offset ' + str(self.pixel_az[match])+'\n')
				logfile.write('Using el offset ' + str(self.pixel_el[match])+'\n')
				pixel_az = az + self.pixel_az[match]
				pixel_el = el + self.pixel_el[match]
			else:
				pixel_az = az
				pixel_el = el

			# Calculate the positions for this pixel
			print('Calculate skypos...')
			skypos = self.calc_positions(pixel_az, pixel_el, t.jd)
			print('Calculate healpix positions...')
			healpix_pixel, centralpos = self.calc_healpix_pixels(skypos)
			print('Plot positions...')
			plot_tod(pixel_az,self.outdir+prefix+plotext+'/plot_az_'+str(pix)+'.png')
			plot_tod(pixel_el,self.outdir+prefix+plotext+'/plot_el_'+str(pix)+'.png')
			plot_tod(skypos.ra,self.outdir+prefix+plotext+'/plot_ra_'+str(pix)+'.png')
			plot_tod(skypos.dec,self.outdir+prefix+plotext+'/plot_dec'+str(pix)+'.png')
			# plot_tod(pa,self.outdir+prefix+plotext+'/plot_pa.png')


			plot_tod(dataset[:,pix],self.outdir+prefix+plotext+'/plot_'+str(pix)+'.png')

			plot_ground(az[start:end],dataset[start:end,pix],self.outdir+prefix+plotext+'/plot_'+str(pix)+'_az_binned.png')


			# High-pass filter the data
			# dataset[:,pix+1] = signal.sosfilt(sos, dataset[:,pix+1])
			print('Subtract baselines...')
			num_to_average = 100
			if 'GB02' in prefix:
				num_to_average = 500
			dataset_sub = subtractbaseline(dataset[:,pix_A],option=0,navg=num_to_average)
			dataset_sub2 = subtractbaseline(dataset[:,pix_P],option=0,navg=num_to_average)

			# # Correct for any offset
			# offset = np.median(dataset[:,pix+1])
			# print(offset)
			# dataset[:,pix+1] -= offset

			# If we're looking at phase, unwind if needed
			# dataset[dataset[:,pix+1]>np.pi,pix+1] -= 2.0*np.pi

			# If we're highly negative, then we probably need to invert.
			# if np.max(dataset[:,pix+1]) < -np.min(dataset[:,pix+1])/2.0:
			# 	dataset[:,pix+1] = -dataset[:,pix+1]

			# Also rescale to 1.0 for now
			# scale = np.max(dataset[:,pix+1])
			# scale2 = np.min(dataset[:,pix+1])
			# if np.abs(scale) < np.abs(scale2):
			# 	scale = scale2
			# dataset[:,pix+1] /= scale

			plot_tod(dataset_sub,self.outdir+prefix+plotext+'/plot_'+str(pix)+'_I_rescale.png')
			plot_tod(dataset_sub2,self.outdir+prefix+plotext+'/plot_'+str(pix)+'_Q_rescale.png')
			plot_val_tod(az[start:end],dataset_sub[start:end],self.outdir+prefix+plotext+'/plot_'+str(pix)+'_az.png')
			plot_val_tod(az[start:end],dataset[start:end,pix_P],self.outdir+prefix+plotext+'/plot_'+str(pix)+'_az_raw.png')

			try:
				maxpix = self.get_second_brightest(dataset_sub)
				maxpix2 = self.get_second_brightest(dataset_sub2)
			except:
				logfile.write('Something odd happened with getting argmax of the dataset!')
				return 0
			logfile.write('At array index ' + str(maxpix)+'\n')
			logfile.write('Maximum value in A is ' + str(dataset_sub[maxpix])+'\n')
			logfile.write('At azimuth ' + str(pixel_az[maxpix])+'\n')
			logfile.write('At elevation ' + str(pixel_el[maxpix])+'\n')
			logfile.write('Maximum value in phase is ' + str(dataset_sub2[maxpix2])+'\n')
			logfile.write('At azimuth ' + str(pixel_az[maxpix2])+'\n')
			logfile.write('At elevation ' + str(pixel_el[maxpix2])+'\n')
			# moonpos = get_moon(t[maxpix],self.telescope)
			# logfile.write(str(moonpos)+'\n')
			logfile.write('Moon was at azimuth' + str(moonpos_az[maxpix])+'\n')
			logfile.write('Moon was at elevation ' + str(moonpos_el[maxpix])+'\n')
			logfile.write('At RA in A' + str(skypos[maxpix].ra.deg)+'\n')
			logfile.write('At Dec in A' + str(skypos[maxpix].dec.deg)+'\n')
			logfile.write('At RA in phase' + str(skypos[maxpix2].ra.deg)+'\n')
			logfile.write('At Dec in phase' + str(skypos[maxpix2].dec.deg)+'\n')
			logfile.write('Moon was at RA ' + str(moonpos_ra[maxpix])+'\n')
			logfile.write('Moon was at Dec ' + str(moonpos_dec[maxpix])+'\n')
			az_correction = pixel_az[maxpix]-moonpos_az[maxpix]#.az.deg
			el_correction = pixel_el[maxpix]-moonpos_el[maxpix]#.alt.deg
			logfile.write('Difference in A is ' + str(az_correction) + ' in azimuth and ' + str(el_correction) + ' in elevation.\n')
			az_correction2 = pixel_az[maxpix2]-moonpos_az[maxpix2]#.az.deg
			el_correction2 = pixel_el[maxpix2]-moonpos_el[maxpix2]#.alt.deg
			logfile.write('Difference in phase is ' + str(az_correction2) + ' in azimuth and ' + str(el_correction2) + ' in elevation.\n')
			# exit()

			print('Get maps...')
			skymap, hitmap = self.get_skymap(healpix_pixel,dataset_sub,flag)
			skymap2 = self.get_skymap(healpix_pixel,dataset_sub2,flag,rethitmap=False)
			self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_'+str(pix)+'_A.fits')
			self.write_healpix_map(hitmap,prefix,self.outdir+prefix+'/hitmap_'+str(pix)+'_P.fits')

			hp.mollview(skymap,title='Pixel ' + str(pix))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_A.png')
			plt.clf()
			plt.close()
			hp.mollview(skymap2,title='Pixel ' + str(pix))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_P.png')
			plt.clf()
			plt.close()
			hp.mollview(hitmap,title='Pixel ' + str(pix))
			plt.savefig(self.outdir+prefix+'/hitmap_'+str(pix)+'.png')
			plt.clf()
			plt.close()

			hp.gnomview(skymap,rot=centralpos,reso=23,title='Pixel ' + str(pix))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_A_view.png')
			plt.clf()
			plt.close()
			hp.gnomview(skymap2,rot=centralpos,reso=23,title='Pixel ' + str(pix))
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_P_view.png')
			plt.clf()
			plt.close()

			skymap[skymap != hp.pixelfunc.UNSEEN] -= np.median(skymap[skymap != hp.pixelfunc.UNSEEN])
			std = np.std(skymap[skymap != hp.pixelfunc.UNSEEN])
			hp.mollview(skymap,max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_A_std.png')
			plt.clf()
			plt.close()

			hp.gnomview(skymap,rot=centralpos,reso=23,title='Channel ' + str(pix),max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_A_view_std.png')
			plt.clf()
			plt.close()

			skymap2[skymap2 != hp.pixelfunc.UNSEEN] -= np.median(skymap2[skymap2 != hp.pixelfunc.UNSEEN])
			std = np.std(skymap2[skymap2 != hp.pixelfunc.UNSEEN])
			hp.mollview(skymap2,max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_P_std.png')
			plt.clf()
			plt.close()

			hp.gnomview(skymap2,rot=centralpos,reso=23,title='Channel ' + str(pix),max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix)+'_P_view_std.png')
			plt.clf()
			plt.close()

			# Calculate the positions for this pixel, including the offset calculated from the moon
			print('Moon and sun separation...')
			moondist = skypos.separation(SkyCoord(moonpos_ra*u.deg,moonpos_dec*u.deg,frame='icrs'))
			sundist = skypos.separation(SkyCoord(sunpos_ra*u.deg,sunpos_dec*u.deg,frame='icrs'))
			radist = skypos.ra.deg-moonpos_ra
			decdist = skypos.dec.deg-moonpos_dec
			usevals = np.where(np.sqrt(radist**2+decdist**2) < 10.0)
			if sum(usevals[0]) > 10:
				healpix_pixel2 = hp.ang2pix(self.nside, (np.pi/2)-decdist*np.pi/180.0, radist*np.pi/180.0)
				skymap, hitmap = self.get_skymap(healpix_pixel2,dataset_sub,flag)
				skymap2 = self.get_skymap(healpix_pixel2,dataset_sub2,flag,rethitmap=False)
				self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_moon_'+str(pix)+'_A.fits')
				hp.mollview(skymap,max=3.0*std,min=-3.0*std)
				plt.savefig(self.outdir+prefix+'/skymap_moon_'+str(pix)+'_A.png')
				plt.clf()
				plt.close()
				self.write_healpix_map(skymap2,prefix,self.outdir+prefix+'/skymap_moon_'+str(pix)+'_P.fits')
				hp.mollview(skymap2,max=3.0*std,min=-3.0*std)
				plt.savefig(self.outdir+prefix+'/skymap_moon_'+str(pix)+'_P.png')
				plt.clf()
				plt.close()
				self.write_healpix_map(hitmap,prefix,self.outdir+prefix+'/hitmap_moon_'+str(pix)+'.fits')
				hp.mollview(hitmap)
				plt.savefig(self.outdir+prefix+'/hitmap_moon_'+str(pix)+'_A.png')
				plt.clf()
				plt.close()

			del dataset_sub
			del dataset_sub2
			gc.collect()

			# Write out the TOD for this pixel to disk
			hdr = fits.Header()
			hdr['INSTRUM'] = 'GroundBIRD'
			hdr['REDCODE'] = 'gbreduce'
			hdr['VERSION'] = str(self.ver)
			hdr['OBSNAME'] = prefix[:-1]
			hdr['PIX'] = str(pix)
			hdr['CENFREQ'] = str(freqs[int(pix)])
			hdr['ARGA'] = str(swp_params[int(pix)]['arga'])
			hdr['ABSA'] = str(swp_params[int(pix)]['absa'])
			hdr['TAU'] = str(swp_params[int(pix)]['tau'])
			hdr['FR'] = str(swp_params[int(pix)]['fr'])
			hdr['QR'] = str(swp_params[int(pix)]['Qr'])
			hdr['QC'] = str(swp_params[int(pix)]['Qc'])
			hdr['PHI'] = str(swp_params[int(pix)]['phi0'])
			hdr['C'] = str(swp_params[int(pix)]['c'])

			if match >= 0:
				hdr['AZOFF'] = str(self.pixel_az[match])
				hdr['ELOFF'] = str(self.pixel_el[match])
			else:
				hdr['AZOFF'] = str(0.0)
				hdr['ELOFF'] = str(0.0)

			self.write_fits(Angle(skypos.ra).degree, Angle(skypos.dec).degree, t.mjd, dataset[:,pix_A], dataset[:,pix_P], hdr, self.outdir+prefix+prefix[:-1]+'_'+str(int(pix))+'_tod.fits',moondist=moondist.deg,sundist=sundist.deg,flag=flag)
			# ra_col = fits.Column(name='ra',format='E',array=np.array(Angle(skypos.ra).degree))
			# dec_col = fits.Column(name='dec',format='E',array=np.array(Angle(skypos.dec).degree))
			# jd_col = fits.Column(name='mjd',format='D',array=t.mjd)
			# col_list = [ra_col, dec_col, jd_col]
			# # for i in range(0,numpix*4):
			# # 	col_list.append(fits.Column(name='ch'+str(i),format='E',array=dataset[:,i]))
			# col_list.append(fits.Column(name='I',format='E',array=dataset[:,pix-2]))
			# col_list.append(fits.Column(name='Q',format='E',array=dataset[:,pix-1]))
			# col_list.append(fits.Column(name='A',format='E',array=dataset[:,pix]))
			# col_list.append(fits.Column(name='P',format='E',array=dataset[:,pix+1]))
			# cols = fits.ColDefs(col_list)
			# hdu = fits.BinTableHDU.from_columns(cols)
			#
			#
			# primary_hdu = fits.PrimaryHDU(header=hdr)
			# hdul = fits.HDUList([primary_hdu, hdu])
			# hdul.writeto(self.outdir+prefix+prefix[:-1]+'_'+str(int(pix/4))+'_tod.fits',overwrite=True)

		del dataset
		del skypos
		del moondist
		del sundist
		del radist
		del decdist
		del moonpos_az
		del moonpos_el
		del skymap
		del skymap2
		del hitmap
		gc.collect()

		skymap = np.zeros(self.npix, dtype=np.float)
		for i in range(0,len(healpix_pixel)):
			skymap[healpix_pixel[i]] = az[i]
		self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_az.fits')
		logfile.close()



		return
