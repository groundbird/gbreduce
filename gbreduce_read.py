#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to reading GroundBIRD data, separate from the class
# 
# Version history:
#
# 30-Sep-2019  M. Peel       Started

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.fftpack
from scipy import signal, optimize
import os
from rhea_comm.packet_reader import read_file, get_length
import time
import datetime

# From rhea_comm.reader_tod.header_read()
def read_rhea_header(fname):
	time = 0
	data = 0
	test = read_file(fname, length = 1)
	print(test)
	for val in test:
		print(val)
	# for t, n_r, s_off, d in read_file(fname, length = 1):
	# 	time = t
	# 	data = [v for i, v in enumerate(d) if i % 2 == 0]
	# 	pass
	return time , data

def parse_filename(filename,quiet=True):
	if '/' in filename:
		filename = filename.split('/')[-1]
	fileformat = ''
	freqs = []
	freqwidth = 0
	freqstep = 0
	date = 0
	samplespeed = 0
	if 'swp' in filename:
		fileformat = 'swp'
		# mulswp_+003.000MHzWidth_+000.001MHzStep_-082.740MHz_-080.740MHz_+016.162MHz_+018.162MHz_+029.276MHz_+031.276MHz_+078.775MHz_+080.775MHz_2019-0921-053919.rawdata
		testfilename = filename.replace('mulswp_','').replace('.rawdata','')
		if not quiet:
			print(testfilename)
		test = testfilename.split('_')
		if not quiet:
			print(test)
		for line in test:
			if 'Width' in line:
				if 'MHz' in line:
					freqwidth = float(line.replace('MHz','').replace('Width',''))*1e6
				else:
					print("Unknown frequency unit for width! Line: " + str(line))
			elif 'Step' in line:
				if 'MHz' in line:
					freqstep = float(line.replace('MHz','').replace('Step',''))
				else:
					print("Unknown frequency unit for step! Line: " + str(line))*1e6
			elif 'MHz' in line:
				freqs.append(float(line.replace('MHz',''))*1e6)
		date = time.mktime(datetime.datetime.strptime(test[-1], "%Y-%m%d-%H%M%S").timetuple())
		if not quiet:
			print(fileformat)
			print(freqwidth)
			print(freqstep)
			print(freqs)
			print(date)
		return [fileformat, freqs, freqwidth, freqstep, date]
	elif 'tod' in filename:
		fileformat = 'tod'
		#tod_-082.664MHz_-080.664MHz_+016.174MHz_+018.174MHz_+029.302MHz_+031.302MHz_+078.873MHz_+080.873MHz_0001kSPS_2019-0921-052531.rawdata
		testfilename = filename.replace('tod_','').replace('.rawdata','')
		if not quiet:
			print(testfilename)
		test = testfilename.split('_')
		if not quiet:
			print(test)
		for line in test:
			if 'kSPS' in line:
				samplespeed = float(line.replace('kSPS',''))*1e3
			elif 'MHz' in line:
				freqs.append(float(line.replace('MHz',''))*1e6)
		date = time.mktime(datetime.datetime.strptime(test[-1], "%Y-%m%d-%H%M%S").timetuple())
		if not quiet:
			print(fileformat)
			print(samplespeed)
			print(freqs)
			print(date)
		return [fileformat, freqs, samplespeed, date]
	return []

def read_rhea_swp_data(fname, length=None, offset=0):
	inputdata = read_file(fname, length = length, offset = offset, sync=True)
	freqset = []
	dataset = []
	i = 0
	chan_accumulation = []
	for datum in inputdata:
		if datum[0] == 0:
			# We have frequency data
			# print('Frequencies:')
			# print(datum)
			freqset.append(datum[1])
			# print(len(datum[1]))
			if i != 0:
				dataset.append(chan_accumulation / 10)
			chan_accumulation = np.zeros(len(datum[1]))
		else:
			# We have accumulation data
			# print(datum)
			chan_accumulation = chan_accumulation + datum[1]
		i += 1
		# if i > 20:
		# 	exit()
	# We need one last dataset appending.
	dataset.append(chan_accumulation / 10)
	return (np.asarray(freqset,dtype=np.float64), np.asarray(dataset,dtype=np.float64))

def read_rhea_data(fname, length=None, offset=0):
	inputdata = read_file(fname, length = length, offset = offset, sync=True)
	dataset = []
	i = 0
	for datum in inputdata:
		# print(datum)
		# i += 1
		# if i > 20:
		# 	exit()

		# Skip the first data point
		# if i == 1:
			# continue
		# print(np.shape(datum))
		# print(datum)
		# print(len(datum))
		datum2 = datum[1].copy()
		datum2.insert(0,datum[0])
		dataset.append(datum2)

		# datum2 = np.asarray(datum[1])
		# datum2 = np.insert(datum2,0,datum[0])
		# print(datum2)
		# exit(0)
		# dataset = np.append(datum2)
		# print(dataset)
	return np.asarray(dataset,dtype=np.float64)
	# return [time, data, tmp, tmp1]
