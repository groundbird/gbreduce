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
import astropy_speedups
from scipy import optimize, signal
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
from matplotlib.gridspec import GridSpec
from merge_az import merge_az

class gbreduce:
	def __init__(self,datadir="/net/nas/proyectos/quijote2/tod/",outdir='',azdir='',eldir='',pixinfo='',elcompressed=False,domedir='',tempdir='',nside=512,datarate=1000,use_mkidpylibs=True):
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
		self.focallength = 600.0 # mm
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

		# Check to see if the elevation is under 90°, and set it to 90° if not
		el[el>90.0]=90.0

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


	def runset(self, subdir='',ext='',skipfirst=0):
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

		if skipfirst == 0:
			trip = 1
		else:
			trip = 0
		for item in todolist:
			if trip < skipfirst:
				trip += 1
				continue
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
				swpfile = 'swp'+line+'rawdata'
				todfile = 'tod'+line+'rawdata'
				kidparams = 'kids'+line+'list'
				swp_params = self.analyse_swp(name+line.replace('.',''),self.datadir+folder+swpfile,kidparams=self.datadir+folder+kidparams)
				print(swp_params)
				tod_analysis = self.analyse_tod(name+line.replace('.',''),self.datadir+folder+todfile,kidparams=self.datadir+folder+kidparams,swp_params=swp_params,starttime=starttime)
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
			kid = klib.kidana.KidAnalyzer()
			paramset = []
			for i in range(len(centerfreqs)*2):
				try:
					kid.swpfromfile('rhea',filename,lo,i)
					kid.swp.fitIQ(nfwhm=3, fitter='gaolinbg')
					# kid.swp.fitresult.report()
					mkid_fit_params = kid.swp.fitresult.values()

					fig,ax = klib.plotter.plotSwp(kid.swp,kind='rawdata',title='rawdata',color='blue',ls=':',lw=1,marker='.')
					klib.plotter.plotSwp(kid.swp,kind='fitdata.rawdata',ax=ax,color='red',marker='',lw=2,ls='-',)
					fig.savefig(self.outdir+prefix+plotext+'raw_kidfit_'+str(i)+'.png')
					fig.clf()
					fig,ax = klib.plotter.plotSwp(kid.swp,kind='rwdata',title='rwdata',color='blue',ls=':',lw=1,marker='.')
					klib.plotter.plotSwp(kid.swp,kind='fitdata.rwdata',ax=ax,color='red',marker='',lw=2,ls='-',)
					fig.savefig(self.outdir+prefix+plotext+'rw_kidfit_'+str(i)+'.png')

					print(mkid_fit_params)
					# print(mkid_fit_params['c'])
					# params = np.array([mkid_fit_params['absa'], mkid_fit_params['arga'], mkid_fit_params['tau']*1e7, mkid_fit_params['c']*1e8, mkid_fit_params['fr']*1e-9, mkid_fit_params['Qr']*1e-4, mkid_fit_params['Qc']*1e-4, mkid_fit_params['phi0']])
					# print(params)
					paramset.append(mkid_fit_params)
				except:
					# Sweep fitting has failed for some reason, just use the one for the previous channel for now.
					paramset.append(paramset[-1])
			# exit()
		else:
			# dataset = dataset[1300:1800,:]
			# samplefreqs = samplefreqs[1300:1800]
			# print(samplefreqs)
			# print(dataset)
			print(len(samplefreqs))
			# First do a fit to a skewed Lorentizian.
			params = Fit_SkewedLorentizian(samplefreqs[1300:1800], dataset[1300:1800,0]**2+dataset[1300:1800,1]**2)
			params1 = params
			print(params)
			params = Fit_7para(samplefreqs, dataset[:,0]+dataset[:,1]*1j, fitresult=params)
			print(params)

			t = dataset[:,0]**2 + dataset[:,1]**2            
			fig = plt.figure(figsize=(10, 8))
			gs = GridSpec(8, 2, figure=fig)
			ax1 = fig.add_subplot(gs[0, 0])
			ax1.plot(samplefreqs, dataset[:,0])
			ax2 = fig.add_subplot(gs[0, 1])
			ax2.plot(samplefreqs, dataset[:,1])
			ax3 = fig.add_subplot(gs[1:3, :])
			ax3.plot(dataset[:,0], dataset[:,1])
			# if self.fitresult_sp is not None:
			ax3.plot(Fit_7para_Func(params.params, samplefreqs), Fit_7para_Func(params.params, samplefreqs))
			ax3.set_aspect('equal')
			ax4 = fig.add_subplot(gs[4:6, :])
			ax4.plot(samplefreqs, t)
			ax4.plot(samplefreqs, Fit_SkewedLorentizian_Func(params1.params, samplefreqs))
			ax5 = fig.add_subplot(gs[7:9, :])
			ax5.plot(samplefreqs, Fit_SkewedLorentizian_Func(params1.params, samplefreqs, data=t))
			fig.tight_layout()
			plt.show()
			exit()

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

	def analyse_tod(self, prefix, filename, el=[],rpm=2,datarate=1000.0,starttime=0.0,numpix=4, plotlimit=0.0, quiet=False, dofft=False, plottods=True, plotmap=True, dopol=False, plotcombination=True, numfiles=50,swp_params=[],centerfreqs=[],lo=4.99e9,kidparams=[]):

		logfile = open(self.outdir+prefix+"/_log.txt","w")
		logfile.write(str(swp_params))
		# Sort out the output directories
		if prefix[-1] != '/':
			prefix = prefix + "/"
		plotext = 'plots/'
		print(self.outdir+prefix)
		ensure_dir(self.outdir+prefix)
		ensure_dir(self.outdir+prefix+plotext)

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
			centerfreqs = kids['kids_freqs']
			lo = kids['sg_freq']
			on_channels=kids['kids']
			print('Frequencies:')
			print(kids['kids_freqs'])
			print(kids['blinds_freqs'])
			numpix = len(kids['kids_freqs'])

		print(datarate)
		print(centerfreqs)
		print(starttime)


		print(centerfreqs)
		freqs = np.asarray(centerfreqs)
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
		# print(dataset)
		startsync = -1
		i=0
		while startsync == -1:
			startsync = dataset[i,-2]
			i+=1
		precount = i
		print(startsync)
		endsync = dataset[-1,-2]
		print(endsync)

		startsynctime = find_az_synctime(self.azdir, starttime-100, startsync)
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
			dataset, az = merge_az(dataset[:,0:-2], dataset[:,-2], dataset[:,-1], aztimes, azdata, azrot, azoffset)
			timestream = dataset[:,0]

			# This is the old code, still in use for now.
			# az = np.zeros(len(timestream))
			# aztimes_pos = 0
			# for i in range(0,len(timestream)):
			# 	while(aztimes[aztimes_pos] < timestream[i] and aztimes_pos < len(aztimes)-1):
			# 		aztimes_pos += 1
			# 	az[i] = azdata[aztimes_pos-1]

		# # Fetch the dome data
		# domet, domef = fetch_domedata(self.domedir, min(timestream), max(timestream))
		flag = np.ones(len(timestream))
		# if len(domet) >= 1:
		# 	# We can use the dome opening/closing to do an initial cut of the data
		# 	dometimes_pos = 0
		# 	for i in range(0,len(timestream)):
		# 		while(domet[dometimes_pos] < timestream[i] and dometimes_pos < len(domet)-1):
		# 			dometimes_pos += 1
		# 		if dometimes_pos == len(domet):
		# 			flag[i] = int(domef[dometimes_pos])
		# 		else:
		# 			flag[i] = int(domef[dometimes_pos-1])

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

		t = Time(timestream, format='unix', scale='utc')

		# Plot the azimuth, elevation, and flags
		plot_tod(az,self.outdir+prefix+plotext+'/plot_az.png')
		azdiff = np.diff(az)
		azdiff[azdiff < -300] += 360.0
		plot_tod(azdiff,self.outdir+prefix+plotext+'/plot_az_diff.png')
		plot_tod(el,self.outdir+prefix+plotext+'/plot_el.png')
		plot_tod(flag,self.outdir+prefix+plotext+'/plot_flag.png')

		if swp_params != []:
			for pix in range(0,numpix):

				newamp = subtract_blindtone_JS(dataset[:,pix*4+1]+1j*dataset[:,pix*4+2], dataset[:,pix*4+3]+1j*dataset[:,pix*4+4])

				newamp = gao_rewind(freqs[pix],newamp, swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])

				dataset[:,pix*4+3] = abs(newamp)
				# Need to correct the phase of the angle to avoid negatives
				angletmp = -1.0*angle(newamp)
				angletmp = angletmp - np.sign(angletmp)*np.pi
				# And apply the correction for linearity
				dataset[:,pix*4+4] = 2.0*np.tan(angletmp/2.0)

				# rewound = gao_rewind(freqs[pix],dataset[:,pix*4+1]+1j*dataset[:,pix*4+2], swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
				# # rewound_offpeak = gao_rewind(dataset[:,pix*4]+3, dataset[:,pix*4]+4, swp_params[pix]['arga'], swp_params[pix]['absa'], swp_params[pix]['tau'], swp_params[pix]['fr'], swp_params[pix]['Qr'], swp_params[pix]['Qc'], swp_params[pix]['phi0'])
				# dataset[:,pix*2+1] = abs(rewound)
				# dataset[:,pix*2+2] = angle(rewound)
				# # dataset[:,pix*4+2] = abs(rewound_offpeak)
				# # dataset[:,pix*4+3] = angle(rewound_offpeak)

		# sos = signal.butter(10, 0.333, 'hp', fs=1000, output='sos')

		# Calculating moon positions
		print('Moon calculations...')
		print(len(t))

		numinterp = 1000
		moonpos_az, moonpos_el = get_moon_azel(self.telescope, t, numinterp)
		print(len(moonpos_az))
		# moonpos = get_moon(t,self.telescope)
		# moonpos_azel = moonpos.transform_to(AltAz(location=self.telescope,obstime=t))
		print('Moon calculations done.')
		for pix in range(0,numpix*4):
			print(pix)
			logfile.write('\n\nPixel '+str(pix)+'\n')
			start = 100
			end = -100

			# This will need to be calculated for each pixel, but for now just calculate it the first time.
			# if pix%4 == 0:
			if pix == 0:
				# Calculate the positions for this pixel
				skypos = self.calc_positions(az, el, t.jd)
				healpix_pixel, centralpos = self.calc_healpix_pixels(skypos)
				plot_tod(skypos.ra,self.outdir+prefix+plotext+'/plot_ra.png')
				plot_tod(skypos.dec,self.outdir+prefix+plotext+'/plot_dec.png')
				# plot_tod(pa,self.outdir+prefix+plotext+'/plot_pa.png')


			plot_tod(dataset[:,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'.png')

			plot_ground(az[start:end],dataset[start:end,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'_az_binned.png')


			# High-pass filter the data
			# dataset[:,pix+1] = signal.sosfilt(sos, dataset[:,pix+1])
			num_to_average = 100
			if 'GB02' in prefix:
				num_to_average = 500
			dataset[:,pix+1] = subtractbaseline(dataset[:,pix+1],option=0,navg=num_to_average)

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

			plot_tod(dataset[:,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'_rescale.png')
			plot_val_tod(az[start:end],dataset[start:end,pix+1],self.outdir+prefix+plotext+'/plot_'+str(pix+1)+'_az.png')

			try:
				maxpix = np.argmax(np.abs(dataset[start:end,pix+1]))
			except:
				logfile.write('Something odd happened with getting argmax of the dataset!')
				return 0
			logfile.write('At array index ' + str(maxpix)+'\n')
			logfile.write('Maximum value is ' + str(dataset[maxpix,pix+1])+'\n')
			logfile.write('At azimuth ' + str(az[maxpix])+'\n')
			logfile.write('At elevation ' + str(el[maxpix])+'\n')
			moonpos = get_moon(t[maxpix],self.telescope)
			logfile.write(str(moonpos)+'\n')
			logfile.write(str(moonpos_az[maxpix])+'\n')
			logfile.write(str(moonpos_el[maxpix])+'\n')
			az_correction = az[maxpix]-moonpos_az[maxpix]#.az.deg
			el_correction = el[maxpix]-moonpos_el[maxpix]#.alt.deg
			logfile.write('Difference is ' + str(az_correction) + ' in azimuth and ' + str(el_correction) + ' in elevation.')
			# exit()

			skymap = np.zeros(self.npix, dtype=np.float)
			hitmap = np.zeros(self.npix, dtype=np.float)
			for i in range(0,len(healpix_pixel)):
				if flag[i] == 1:
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

			skymap[skymap != hp.pixelfunc.UNSEEN] -= np.median(skymap[skymap != hp.pixelfunc.UNSEEN])
			std = np.std(skymap[skymap != hp.pixelfunc.UNSEEN])
			hp.mollview(skymap,max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix+1)+'_skymap_std.png')
			plt.clf()
	
			hp.gnomview(skymap,rot=centralpos,reso=20,title='Channel ' + str(pix+1),max=3.0*std,min=-3.0*std)
			plt.savefig(self.outdir+prefix+'/skymap_'+str(pix+1)+'view_std.png')
			plt.clf()

			domoon = True
			if domoon == True:
				# Calculate the positions for this pixel, including the offset calculated from the moon
				# skypos2 = self.calc_positions(az+az_correction, el+el_correction, t.jd)
				azdist = az+az_correction - moonpos_az#.az.deg
				eldist = el+el_correction - moonpos_el#.alt.deg
				healpix_pixel2 = hp.ang2pix(self.nside, (np.pi/2)-eldist*np.pi/180.0, azdist*np.pi/180.0)
				skymap = np.zeros(self.npix, dtype=np.float)
				hitmap = np.zeros(self.npix, dtype=np.float)
				for i in range(0,len(healpix_pixel)):
					if flag[i] == 1:
						skymap[healpix_pixel2[i]] = skymap[healpix_pixel2[i]] + dataset[i,pix+1]
						hitmap[healpix_pixel2[i]] = hitmap[healpix_pixel2[i]] + 1
				for i in range(0,len(skymap)):
					if hitmap[i] >= 1:
						skymap[i] = skymap[i]/hitmap[i]
					else:
						skymap[i] = hp.pixelfunc.UNSEEN

				self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_moon_'+str(pix+1)+'.fits')
				self.write_healpix_map(hitmap,prefix,self.outdir+prefix+'/hitmap_moon_'+str(pix+1)+'.fits')


		for i in range(0,len(healpix_pixel)):
			skymap[healpix_pixel[i]] = az[i]
		self.write_healpix_map(skymap,prefix,self.outdir+prefix+'/skymap_az.fits')
		logfile.close()
		# exit()
		# # Write out the TOD to disk
		# ra_col = fits.Column(name='ra',format='E',array=np.array(Angle(skypos.ra).degree))
		# dec_col = fits.Column(name='dec',format='E',array=np.array(Angle(skypos.dec).degree))
		# jd_col = fits.Column(name='mjd',format='E',array=t.jd-2400000.5)
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

	# def get_azimuth_data(self,start,end):
	# 	# Work out which azimuth files we need to read in.
	# 	startdir = datetime.utcfromtimestamp(int(start)).strftime('%Y/%m/%d')
	# 	try:
	# 		filelist = os.listdir(self.azdir + startdir)
	# 	except:
	# 		print('No azimuth data found!')
	# 		return []

	# 	enddir = datetime.utcfromtimestamp(int(end)).strftime('%Y/%m/%d')
	# 	if enddir != startdir:
	# 		# Also need to append these
	# 		try:
	# 			filelist.append(os.listdir(self.azdir + enddir))
	# 		except:
	# 			print('No azimuth data found for the second day!')

	# 	print(filelist)
	# 	return filelist


	# def get_azimuth_data_from_dir(self,dir):
	# 	results = []
	# 	for file in os.listdir(self.aself.datadir):
	# 		if searchtext in file:
	# 			obsname = "-".join(file.split("-")[:-1])
	# 			if obsname not in results:
	# 				results.append(obsname)
	# 	return results
