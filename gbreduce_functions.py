#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to gbreduce, separate from the class
#
# Version history:
#
# 30-Sep-2019  M. Peel       Duplicate from tfgi_functions

import numpy as np
import healpy as hp
import astropy as ap
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.fftpack
from scipy import signal, optimize
import astropy.units as u
from scipy.optimize import curve_fit
import os
from lmfit import minimize, Parameters, fit_report
import gc

def calc_beam_area(beam):
	return (np.pi * (beam*np.pi/180.0)**2)/(4.0*np.log(2.0))

def calc_Tsys(std,B,t):
	return std*np.sqrt(B*t)

def linfit(x, A, B):
	return A*x+B

def fit_kneefreq(freq, param):
	sigma, fknee, alpha = param
	return sigma**2 * (1 + (fknee / freq)**alpha)

def compute_residuals(param, data, freq):
	model = fit_kneefreq(freq, param)
	residual = np.log(data / model)
	return residual

def fit_skydip(el, param):
	return param[0] + param[1]/(np.sin(el*np.pi/180.0))# + param[2]*el

def compute_residuals_skydip(param, el, data):
	model = fit_skydip(el, param)
	residual = data - model
	return residual

def ensure_dir(f):
	os.makedirs(f, exist_ok=True)

def fit_gaussian(x, param):
	return param[0]*np.exp(-(x*x)/(2*param[1]*param[1]))+param[2]

def compute_residuals_gaussian(param, x, y):
	return y - fit_gaussian(x,param)

# From mkid_pylibs.fitters.py by Honda
def gao_rewind(x, y, arga, absa, tau, fr, Qr, Qc, phi0):
	tmp = y/absa/np.exp(-1j*(2*np.pi*x*tau-arga))
	return (tmp-1)*Qc/Qr + 0.5

# From mkid_pylibl.modIQ.py by Honda
def subtract_blindtone_JS(mainIQ,blindIQ):
    nz_on = 1.0/np.average(np.abs(mainIQ))*mainIQ
    nz_off = 1.0/np.average(np.abs(blindIQ))*blindIQ
    namp_corrected = np.abs(nz_on) - (np.abs(nz_off) - np.average(np.abs(nz_off)))
    amp_corrected  = namp_corrected * np.average(np.abs(mainIQ))
    nphs_corrected = np.angle(nz_on) - (np.angle(nz_off) - np.average(np.angle(nz_off)))
    I_dash = amp_corrected*np.cos(nphs_corrected)
    Q_dash = amp_corrected*np.sin(nphs_corrected)
    return I_dash + 1j*Q_dash

def fit_mkids(f, param):
	# a, w, t, c, fr, qr, qc, p0 = param
	# return a * np.exp(-2.0*np.pi*1j*(f*t-w)) * (1 + c*(f-fr) - (qr/(qc*np.exp(1j*p0)))/(1.0+2*1j*qr*((f-fr)/fr)))
	absa, arga, tau, c, fr, Qr, Qc, phi0 = param
	# print('absa = ' + str(absa))
	# print('arga = ' + str(arga))
	# print('tau = ' + str(tau))
	# print('c = ' + str(c))
	# print('fr = ' + str(fr))
	# print('Qr = ' + str(Qr))
	# print('Qc = ' + str(Qc))
	# print('phi0 = ' + str(phi0))
	# exit()
	x = f
	# print(x)
	# exit()
	return (absa * np.exp(-1j*(2*np.pi*x*(tau*1e-7) - arga))*(1+(c*1e-8)*(x-(fr*1e9))-(Qr*1e4)/(Qc*1e4)*np.exp(1j*phi0)/(1+2*1j*(Qr*1e4)*((x-fr*1e9)/(fr*1e9)))))

def compute_residuals_mkids(param, x, y):
	diff = real_to_complex(y) - fit_mkids(x,param)
	# y1d = np.zeros(y.size*2, dtype = np.float64)
	# y1d[0:y1d.size:2] = diff.real
	# y1d[1:y1d.size:2] = diff.imag
	return complex_to_real(diff)

def lnlike(p, x, y, yerr):
	model = fit_mkids(x,p)
	inv_sigma2 = 1.0/(real_to_complex(yerr)**2) # + model**2*np.exp(2*lnf))
	value = -0.5*(np.sum((real_to_complex(y)-model)**2*inv_sigma2 - np.log(inv_sigma2)))
	return value

def lnprior(p, bounds):
	for i in range (0,len(p)):
		if p[i] < bounds[i][0] or p[i] > bounds[i][1]:
			return -np.inf
	# If we've got to here, we're within the bounds.
	return 0.0

def lnprob(p, bounds, x, y, yerr):
	lp = lnprior(p, bounds)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(p, x, y, yerr)

def real_to_complex(z):      # real vector of length 2n -> complex of length n
	return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z):      # complex vector of length n -> real of length 2n
	return np.concatenate((np.real(z), np.imag(z)))


def plot_tod(data, outputname,formatstr='b.'):
	plt.plot(data,formatstr)
	plt.xlabel('Samples')
	plt.ylabel('Power')
	plt.savefig(outputname)
	plt.clf()
	plt.close()
	gc.collect()
	return

# Plot the TODs against a given set of vals, e.g. az or el.
def plot_val_tod(val, data, outputname,formatstr='b.'):
	plt.plot(val,data,formatstr)
	plt.savefig(outputname)
	plt.clf()
	plt.close()
	gc.collect()
	return



## Fitting parameters for swp data, from Satoru, 20190925_main_v11.py
# fitter_nitta SkewedLorentizian (first fitting)
def Fit_SkewedLorentizian_initparams(f, t):
	A1_init = (t[0] + t[1])/2
	A2_init = 0.0
	A3_init = np.amin(t) - A1_init
	A4_init = 0.0
	fr_init = f[np.argmin(t)]
	Qr_init = 1000
	initparams = {'A1': A1_init, 'A2': A2_init, 'A3': A3_init, 'A4': A4_init, 'fr': fr_init, 'Qr': Qr_init}
	return initparams

def Fit_SkewedLorentizian_Func(pars, f, data=None):
	vals = pars.valuesdict()
	A1 = vals['A1']
	A2 = vals['A2']
	A3 = vals['A3']
	A4 = vals['A4']
	fr = vals['fr']
	Qr = vals['Qr']
	model = (A1 + A2*(f-fr) + (A3 + A4*(f-fr))/(1 + 4*Qr*Qr*((f-fr)/fr)**2))
	# return the residual
	if data is None:
		return model
	return model - data

def Fit_SkewedLorentizian(f, t):
	initparams = Fit_SkewedLorentizian_initparams(f, t)
	fit_params = Parameters()
	fit_params.add('A1', value = initparams['A1'], min = 100)  # Background
	fit_params.add('A2', value = initparams['A2'])#, min = -0.5, max = 0.5)  # Slope
	fit_params.add('A3', value = initparams['A3'], max = 0)  # Lowest point
	fit_params.add('A4', value = initparams['A4'])
	fit_params.add('fr', value = initparams['fr'], min = initparams['fr'] - 0.2e6, max = initparams['fr'] + 0.2e6)
	fit_params.add('Qr', value = initparams['Qr'], min = 0, max = 1e6)
	fit_params.add('Qi', expr = 'Qr*sqrt(A1/(A1 + A3))')
	fit_params.add('Qc', expr = '1/(1/Qr - 1/Qi)')
	out = minimize(Fit_SkewedLorentizian_Func, fit_params, args=(f,), kws={'data': t}, nan_policy='omit')
	# print the fitting result
	print('Fit_SkewedLorentizian')
	print(fit_report(out))
	return out

# fitter_nitta 7para (second fitting)
def Fit_7para_initparams(f, t, fitresult=None):
	t_abs = np.abs(t)
	t_arg = np.angle(t)
	t_unwrap = np.unwrap(t_arg)
	if fitresult is not None:
		a_init = np.sqrt(fitresult.params['A1'].value)
		alpha_init = t_unwrap[0]
		tau_init = (t_unwrap[0] - t_unwrap[-1])/((2*np.pi)*(f[-1] - f[0]))
		phi0_init = t_unwrap[0] - 2*np.pi*fitresult.params['fr'].value*tau_init
		fr_init = fitresult.params['fr'].value
		Qr_init = fitresult.params['Qr'].value
		Qc_init = fitresult.params['Qc'].value
	else:
		a_init = (t_abs[0] + t_abs[-1])/2
		alpha_init = t_unwrap[0]
		tau_init = (t_unwrap[0] - t_unwrap[-1])/((2*np.pi)*(f[-1] - f[0]))
		phi0_init = t_unwrap[0] - 2*np.pi*np.mean(f)*tau_init
		fr_init = np.mean(f)
		Qr_init = 1000
		Qc_init = 1000.1
	initparams = {'a': a_init, 'alpha': alpha_init, 'tau': tau_init, 'phi0': phi0_init, 'fr': fr_init, 'Qr': Qr_init, 'Qc': Qc_init}
	return initparams

def Fit_7para_Func(pars, f, data=None):
	vals = pars.valuesdict()
	a = vals['a']
	alpha = vals['alpha']
	tau = vals['tau']
	phi0 = vals['phi0']
	fr = vals['fr']
	Qr = vals['Qr']
	Qc = vals['Qc']
	model = a * np.exp(1j*alpha) * np.exp(-2*np.pi*1j*f*tau) * (1 - (Qr/Qc*np.exp(1j*phi0))/(1 + 2*1j*Qr*(f-fr)/fr))
	# return the residual
	if data is None:
		return model
	resid = model - data
	return resid.view(np.float)

def Fit_7para(f, t, **kwargs):
	initparams = Fit_7para_initparams(f, t, **kwargs)
	fit_params = Parameters()
	fit_params.add('a', value = initparams['a'], min = 0)#, max = 50)
	fit_params.add('alpha', value = initparams['alpha'])#, min = -np.pi, max = np.pi)
	fit_params.add('tau', value = initparams['tau'])#, min = -5e-9, max = 5e-9)
	fit_params.add('phi0', value = initparams['phi0'])#, min = -np.pi, max = np.pi)
	fit_params.add('fr', value = initparams['fr'])
	fit_params.add('Qr', value = initparams['Qr'], min = 1e2, max = 1e8)
	fit_params.add('Qc', value = initparams['Qc'], min = 1.01e2, max = 1e8)
	print('hi')
	print(fit_params)
	print(1/(1/initparams['Qr']-1/initparams['Qc']))
	fit_params.add('Qi', expr = '1/(1/Qr - 1/Qc)')
	print('hi')
	print(fit_params)
	out = minimize(Fit_7para_Func, fit_params, args=(f,), kws={'data': t}, nan_policy='omit')
	# print the fitting result
	print('Fit_7para')
	print(fit_report(out))
	return out

def subtractbaseline(data, option=0, navg=1500):
	# Crude baseline removal
	npoints = len(data)
	# print(npoints)
	for i in range(0,npoints//navg):
		start = i*navg
		end = ((i+1)*navg)
		# print('start:' + str(start) + ', end: ' + str(end))
		if option == 0:
			data[start:end] = data[start:end] - np.median(data[start:end])
		else:
			xline = range(0,len(data[start:end]))
			if len(xline) != 0:
				A,B=curve_fit(linfit,xline,data[start:end])[0]
				# print(A,B)
				data[start:end] = data[start:end] - (A*xline+B)

	if option == 0:
		data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - np.median(data[(npoints//navg)*navg-1:])
	else:
		xline = range(0,len(data[(npoints//navg)*navg-1:]))
		if len(xline) != 0:
			A,B=curve_fit(linfit,xline,data[(npoints//navg)*navg-1:])[0]
			# print(A,B)
			data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - (A*xline+B)

	return data

def plot_ground(az, data, outfile,nazpoints=360):
	# Define a histogram of the azimuth range
	azbins = np.arange(0,360.0,nazpoints/360)
	azvals = np.zeros(nazpoints)
	azhits = np.zeros(nazpoints)
	numdata = len(data)
	if len(az) < numdata:
		numdata = len(az)
	# Run through the data to populate it
	# Need to handle cases where az=360 better!
	for i in range(0,numdata):
		try:
			azvals[int(np.round(az[i])*360/nazpoints)] += data[i]
			azhits[int(np.round(az[i])*360/nazpoints)] += 1
		except:
			null = 0

	# Average it
	azvals /= azhits

	# Plot it
	plt.plot(azbins,azvals)
	plt.xlim(0,360.0)
	plt.xlabel('Azimuth')
	plt.ylabel('Amplitude')
	plt.savefig(outfile)
	plt.clf()
	plt.close()
	np.savetxt(outfile+'.txt',[azbins, azvals])

	del azbins
	del azvals
	del azhits
	gc.collect()
	return


def get_body_azel_radec(location, times, numinterp=1000,body='moon'):
	bodypos = ap.coordinates.get_body(body,times[0::numinterp],location)
	bodypos_azel = bodypos.transform_to(ap.coordinates.AltAz(location=location,obstime=times[0::numinterp]))
	az = np.interp(times.mjd, times[0::numinterp].mjd, bodypos_azel.az.deg)
	el = np.interp(times.mjd, times[0::numinterp].mjd, bodypos_azel.alt.deg)
	ra = np.interp(times.mjd, times[0::numinterp].mjd, bodypos.ra.deg)
	dec = np.interp(times.mjd, times[0::numinterp].mjd, bodypos.dec.deg)
	del bodypos_azel
	del bodypos
	gc.collect()
	return az, el, ra, dec
