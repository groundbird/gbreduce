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


def plot_tod(data, outputname,formatstr='b.'):
	plt.plot(data,formatstr)
	plt.xlabel('Samples')
	plt.ylabel('Power')
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return

# Plot the TODs against a given set of vals, e.g. az or el.
def plot_val_tod(val, data, outputname):
	plt.plot(val,data,'b.')
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return