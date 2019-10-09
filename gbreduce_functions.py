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
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.fftpack
from scipy import signal, optimize
import astropy.units as u
from scipy.optimize import curve_fit
import os


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
	plt.close()
	plt.clf()
	return

# Plot the TODs against a given set of vals, e.g. az or el.
def plot_val_tod(val, data, outputname,formatstr='b.'):
	plt.plot(val,data,formatstr)
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return