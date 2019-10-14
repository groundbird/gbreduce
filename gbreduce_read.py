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
