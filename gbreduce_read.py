#!/usr/bin/env python3
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
import packet_reader
import time
import datetime
import json
import eldata
import glob
from enum import Enum
from pathlib import Path
import scipy.stats as stats
import RotLog_file
import RotLog_file_oldversion
import pytz     ## pip install pytz
import gc
from itertools import compress
from struct import pack, unpack
import struct
# import multiprocessing
from numba import njit

# Make sure we are using UTC!
# mytz = pytz.timezone('UTC')             ## Set your timezone
# datetime = mytz.normalize(mytz.localize(datetime, is_dst=False))  ## Set is_dst accordingly

def get_pixinfo(filename, focallength, angle):
	npix,mod,mod_pix,x,y,theta,phi,omt1,omt2, pixelfreq = np.loadtxt(filename,unpack=True)
	pixelfreq2 = pixelfreq.tolist()
	x = (x/focallength)*180.0/np.pi
	y = (y/focallength)*180.0/np.pi
	new_x = x*np.cos(angle*np.pi/180.0) - y*np.sin(angle*np.pi/180.0)
	new_y = x*np.sin(angle*np.pi/180.0) + y*np.cos(angle*np.pi/180.0)

	return new_x, new_y, pixelfreq2

def get_az_data(filename):
	rotlog = RotLog_file.RotLog_file(filename)
	unixtime = float(rotlog.start_time)
	datalength = rotlog.length
	times, nr, roff, vals = rotlog.read_file()
	# for i in range(0,1000):
	# 	print(str(times[i]) + '-' + str(nr[i]) + ' - ' + str(roff[i]) + ' - ' + str(vals[i]))
	# exit()
	# numbytes, version, unixtime, headertxt = tmped.get_header()
	starttime = float(times[0])
	times_new = (np.asarray(times).astype(float) - starttime)*0.001 + unixtime
	# times_new = []
	# for i in range(len(times)):
	# 	times_new.append((float(times[i]) - starttime)*0.001 + unixtime)
	# print(filename)
	# print(times[0])
	# print(times[-1])
	return list(times_new), vals, nr, roff


def compress_and_plot_azdata(indir,out_prefix,reload=False,savedata=True):
	if reload:
		times, vals = np.loadtxt(indir+out_prefix+'.txt.gz',unpack=True)
	else:
		files = sorted(Path(indir).rglob('*.dat.xz'))
		# print(files)
		times = []
		vals = []
		times_sync = []
		vals_sync = []
		for filename in files:
			print(filename)
			timeset, valset = get_az_data(filename)
			times = times + timeset
			vals = vals + valset
		# Calculate zenith distance, and save to disk
		if savedata:
			np.savetxt(indir+out_prefix+'.txt.gz',np.transpose([times, vals]))

	times = np.asarray(times)
	vals = np.asarray(vals)

	date = datetime.datetime.utcfromtimestamp(times[0])
	daystart = datetime.datetime(date.year,date.month,date.day,0,0,0,tzinfo=pytz.timezone('UTC')).timestamp()
	# Plot the encoder data
	plt.figure(figsize=(8.0, 5.0), dpi=100)
	plt.grid(True)
	plt.xlabel('Unix time')
	plt.ylabel('Azimuth value')
	plt.plot((times-daystart)/(60.0*60.0), vals,'b-')
	plt.savefig(indir+out_prefix+'_az.png', dpi=1000)
	plt.clf()

	return



def fetch_azdata(indir, starttime, endtime, compressed=False):
	# We should have unix time as an input, so convert that to Year/Month/day
	# We'll assume that we'll only ever need 3 days of data at most.
	if compressed:
		ext = '.txt.gz'
	else:
		ext = '.dat.xz'

	prevdir = datetime.datetime.utcfromtimestamp(int(starttime)-24*60*60).strftime('%Y/%m/%d')
	# print(indir)
	# print(prevdir)
	try:
		inputlist = os.listdir(indir + prevdir)
		filelist = [prevdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No azimuth data found for the previous day!')
		# return []


	startdir = datetime.datetime.utcfromtimestamp(int(starttime)-1).strftime('%Y/%m/%d')
	# print(indir)
	# print(startdir)
	try:
		inputlist = os.listdir(indir + startdir)
		filelist = [startdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No azimuth data found for the day!')
		# return []

	enddir = datetime.datetime.utcfromtimestamp(int(starttime)+24*60*60).strftime('%Y/%m/%d')
	# print(startdir)
	# print(enddir)
	if enddir != startdir:
		# print('Hi')
		# Also need to append these
		try:
			inputlist = os.listdir(indir + enddir)
			filelist = filelist + [enddir+'/'+f for f in inputlist if ext in f]
		except:
			print('No azimuth data found for the second day!')
	# print(filelist)
	timeset = []
	az = []
	rot = []
	offset = []
	if compressed:
		# We can just read in all the data and return it, since there isn't that much
		for file in filelist:
			if 'sync' not in file:
				times, azimuths = np.loadtxt(indir+file,unpack=True)
				timeset = list(timeset) + list(times)
				az = list(az) + list(azimuths)
	else:
		# We need to be more picky to save on RAM, only return data relevant to the observation.
		# First calculate the start times for each file
		starttimes = []
		use_files = []
		for file in filelist:
			time = file[-18:-12]
			date = file[0:10]
			# print(date[0:4]+' ' + date[5:7] + ' ' + date[8:10] +" " + time[0:2] + ' ' + time[2:4] + ' ' + time[4:6])
			timestamp = datetime.datetime(int(date[0:4]),int(date[5:7]),int(date[8:10]),int(time[0:2]),int(time[2:4]),int(time[4:6]),tzinfo=pytz.timezone('UTC')).timestamp()
			starttimes.append(timestamp)
			# print(starttime)
			# print(timestamp)
			# print(endtime)
			if timestamp > starttime and timestamp < endtime:
				if len(use_files) != 0:
					use_files[-1] = 1
				use_files.append(1)
			else:
				use_files.append(0)
			if timestamp > endtime and np.sum(use_files) == 0:
				# We have the case that the entire dataset is in one file, just use that one.
				if len(use_files) > 1:
					use_files[-2] = 1

		use_files = np.asarray(use_files).astype(int)
		starttimes = np.asarray(starttimes)
		filelist = np.asarray(filelist)
		# print(starttimes[use_files>0])

		for file in filelist[use_files > 0]:
			print(file)
			times, azimuths, nrot, roff = get_az_data(indir+file)
			timeset = list(timeset) + list(times)
			az = list(az) + list(azimuths)
			rot = list(rot) + list(nrot)
			offset = list(offset)+list(roff)

	return timeset, az, rot, offset

# Find the sync time from the azimuth data. Requires an approximate start time to start searching from.
def find_az_synctime(indir, starttime, synctime):
	# We should have unix time as an input, so convert that to Year/Month/day
	# We'll assume that we'll only ever need 3 days of data at most.
	ext = '.dat.xz'
	filelist=[]
	prevdir = datetime.datetime.utcfromtimestamp(int(starttime)-24*60*60).strftime('%Y/%m/%d')
	print(prevdir)
	try:
		inputlist = os.listdir(indir + prevdir)
		filelist = filelist+[prevdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No azimuth data found for previous day!')
		# return []

	startdir = datetime.datetime.utcfromtimestamp(int(starttime)-1).strftime('%Y/%m/%d')
	try:
		inputlist = os.listdir(indir + startdir)
		filelist = filelist+[startdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No azimuth data found for the first day!')
		return []

	enddir = datetime.datetime.utcfromtimestamp(int(starttime)+24*60*60).strftime('%Y/%m/%d')
	if enddir != startdir:
		try:
			inputlist = os.listdir(indir + enddir)
			filelist = filelist + [enddir+'/'+f for f in inputlist if ext in f]
		except:
			print('No azimuth data found for the second day!')

	# We need to be more picky to save time.
	# First calculate the start times for each file
	starttimes = []
	use_files = []
	for file in filelist:
		time = file[-18:-12]
		date = file[0:10]
		# print(date[0:4]+' ' + date[5:7] + ' ' + date[8:10] +" " + time[0:2] + ' ' + time[2:4] + ' ' + time[4:6])
		timestamp = datetime.datetime(int(date[0:4]),int(date[5:7]),int(date[8:10]),int(time[0:2]),int(time[2:4]),int(time[4:6]),tzinfo=pytz.timezone('UTC')).timestamp()
		starttimes.append(timestamp)
		# print(starttime)
		# print(timestamp)
		# print(endtime)
		if timestamp > starttime:
			if len(use_files) != 0:
				use_files[-1] = 1
			use_files.append(1)
		else:
			use_files.append(0)
	if np.sum(use_files) == 0:
		# We have the case that the entire dataset is in one file, just use that one.
		use_files[-2] = 1

	use_files = np.asarray(use_files).astype(int)
	starttimes = np.asarray(starttimes)
	filelist = np.asarray(filelist)
	# print(starttimes[use_files>0])

	synctimestamp = 0
	for file in filelist[use_files > 0]:
		# print(file)
		rotlog = RotLog_file.RotLog_file(indir+file)
		unixtime = float(rotlog.start_time)
		datalength = rotlog.length
		times, nr, roff, vals = rotlog.read_file()
		if synctime in nr:
			for i in range(0,len(times)):
				if nr[i] == synctime:
					synctimestamp = (float(times[i]) - float(times[0]))*0.001 + unixtime
					break
			break
		else:
			continue

	return synctimestamp


def get_el_data(filename):
	tmped = eldata.ElData(filename)
	numbytes, version, unixtime, headertxt = tmped.get_header()
	times = []
	vals = []
	for i in range(0,tmped._length):
		try:
			data = tmped.get_data(i)
			if data[2].name == 'DATA':
				if starttime == 0:
					starttime = data[0]
				times.append((data[0] - starttime)*0.001 + unixtime)
				vals.append(data[1])
		except:
			# If we find a problem, move onto the next one
			continue

	return times, vals

def get_el_data_condensed(filename):
	tmped = eldata.ElData(filename)
	numbytes, version, unixtime, headertxt = tmped.get_header()
	times = []
	vals = []
	lastval = np.nan # Set to an invalid value to start with
	starttime = 0
	times_sync = []
	vals_sync = []
	length = tmped.get_length()
	# print(length)
	for i in range(0,tmped.get_length()):
		try:
			data = tmped.get_data(i)
			# print(data)
			if data[2].name == 'DATA':
				if starttime == 0:
					starttime = data[0]
				if data[1] != lastval:
					times.append((data[0] - starttime)*0.001 + unixtime)
					vals.append(data[1])
					lastval = data[1]
			elif data[2].name == 'SYNC':
				times_sync.append((data[0] - starttime)*0.001 + unixtime)
				vals_sync.append(data[1])
		except:
			# If we find a problem, move onto the next one
			continue
	# print(len(times))
	# print(len(vals))
	return times, vals, times_sync, vals_sync

def compress_and_plot_eldata(indir,out_prefix,reload=False,savedata=True):
	if reload:
		times, vals, zenithdistance = np.loadtxt(indir+out_prefix+'.txt.gz',unpack=True)
	else:
		files = sorted(Path(indir).rglob('*.dat.xz'))
		print(files)
		times = []
		vals = []
		times_sync = []
		vals_sync = []
		for filename in files:
			print(filename)
			timeset, valset, timesyncset, valsyncset = get_el_data_condensed(filename)
			times = times + timeset
			vals = vals + valset
			times_sync = times_sync + timesyncset
			vals_sync = vals_sync + valsyncset
		# Calculate zenith distance, and save to disk
		zenithdistance = (np.asarray(vals).copy()-9113.0)/900.0
		if savedata:
			np.savetxt(indir+out_prefix+'.txt.gz',np.transpose([times, vals, zenithdistance]))
			np.savetxt(indir+out_prefix+'_sync.txt.gz',np.transpose([times_sync, vals_sync]))

	times = np.asarray(times)
	vals = np.asarray(vals)

	# Plot the encoder data
	plt.figure(figsize=(8.0, 5.0), dpi=100)
	plt.plot(times, vals,'b-')
	plt.savefig(indir+out_prefix+'_encoder.png', dpi=1000)
	plt.clf()

	# Plot the zenith distance
	plt.plot(times, zenithdistance,'b-')
	plt.grid(True)
	plt.xlabel('Unix time')
	plt.ylabel('Zenith distance')
	plt.savefig(indir+out_prefix+'_zenith.png', dpi=1000)
	plt.clf()

	return

def fetch_eldata(indir, starttime, endtime, compressed=False):
	# We should have unix time as an input, so convert that to Year/Month/day
	# We'll assume that we'll only ever need 3 days of data at most.
	if compressed:
		ext = '.txt.gz'
	else:
		ext = '.dat.xz'
	filelist = []
	prevdir = datetime.datetime.utcfromtimestamp(int(starttime)-24*60*60).strftime('%Y/%m/%d')
	# print(indir)
	# print(prevdir)
	try:
		inputlist = os.listdir(indir + prevdir)
		filelist = filelist+[prevdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No elevation data found for the previous day!')
		# return []


	startdir = datetime.datetime.utcfromtimestamp(int(starttime)-1).strftime('%Y/%m/%d')
	# print(indir)
	# print(startdir)
	try:
		inputlist = os.listdir(indir + startdir)
		filelist = filelist+[startdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No elevation data found for the day!')
		# return []

	enddir = datetime.datetime.utcfromtimestamp(int(starttime)+24*60*60).strftime('%Y/%m/%d')
	# print(startdir)
	# print(enddir)
	if enddir != startdir:
		# print('Hi')
		# Also need to append these
		try:
			inputlist = os.listdir(indir + enddir)
			filelist = filelist + [enddir+'/'+f for f in inputlist if ext in f]
		except:
			print('No elevation data found for the second day!')

	timeset = []
	valset = []
	zenith = []
	if compressed:
		# We can just read in all the data and return it, since there isn't that much
		for file in filelist:
			if 'sync' not in file:
				times, vals, zenithdistance = np.loadtxt(indir+file,unpack=True)
				timeset = list(timeset) + list(times)
				valset = list(valset) + list(vals)
				zenith = list(zenith) + list(zenithdistance)
	else:
		# We need to be more picky to save on RAM, only return data relevant to the observation.
		# First calculate the start times for each file
		starttimes = []
		use_files = []
		for file in filelist:
			# print(file)
			time = file[-18:-12]
			date = file[0:10]
			# print(date[0:4]+' ' + date[5:7] + ' ' + date[8:10] +" " + time[0:2] + ' ' + time[2:4] + ' ' + time[4:6])
			timestamp = datetime.datetime(int(date[0:4]),int(date[5:7]),int(date[8:10]),int(time[0:2]),int(time[2:4]),int(time[4:6]),tzinfo=pytz.timezone('UTC')).timestamp()
			starttimes.append(timestamp)
			# print(timestamp)
			# print(starttime)
			# print(endtime)
			if timestamp > starttime and timestamp < endtime:
				if len(use_files) != 0:
					use_files[-1] = 1
				use_files.append(1)
			elif timestamp > starttime and np.sum(use_files) == 0:
				# Catch the case where the observation is entirely within 1 file
				if len(use_files) != 0:
					use_files[-1] = 1
				use_files.append(0)
			else:
				use_files.append(0)
		use_files = np.asarray(use_files).astype(int)
		starttimes = np.asarray(starttimes)
		filelist = np.asarray(filelist)

		print('\n')
		for file in filelist[use_files > 0]:
			print(file)
			# try:
			times, vals, timesync, valsync = get_el_data_condensed(indir+file)
			# except:
				# continue
			timeset = list(timeset) + list(times)
			zenith = list(zenith) + list((np.asarray(vals).copy()-9113.0)/900.0)

	return timeset, zenith


def fetch_domedata(indir, starttime, endtime):
	# We should have unix time as an input, so convert that to Year/Month/day
	# We'll assume that we'll only ever need 2 days of data at most.
	ext = '.dat'
	startdir = datetime.datetime.utcfromtimestamp(int(starttime)-1).strftime('%Y/%m/')
	print(indir)
	print(startdir)
	try:
		inputlist = os.listdir(indir + startdir)
		filelist = [startdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No dome data found!')
		return []

	enddir = datetime.datetime.utcfromtimestamp(int(endtime)+1).strftime('%Y/%m/')
	# print(startdir)
	# print(enddir)
	if enddir != startdir:
		# print('Hi')
		# Also need to append these
		try:
			inputlist = os.listdir(indir + enddir)
			for file in inputlist:
				if file not in filelist:
					filelist.append(enddir+'/'+file)
		except:
			print('No dome data found for the second day!')
	print(filelist)
	timeset = []
	valset = []
	# We can just read in all the data and return it, since there isn't that much
	for file in filelist:
		print(indir+file)
		data = np.loadtxt(indir+file,unpack=False,dtype={'names':('date','unix','val1','val2','val3','val4','water'),'formats':('S1','d','S6','S6','S6','S6','S3')})
		# print(data)
		# print(data.size)
		if data.size > 1:
			for line in data:
				timeset.append(line['unix'])
				# print(line['val1'].decode('ascii'))
				if line['val1'].decode('ascii') == 'Interm' or line['val1'].decode('ascii') == 'Closed':
					# print(0)
					valset.append(0)
				else:
					# print(1)
					valset.append(1)

	return timeset, valset

def fetch_tempdata(indir, starttime, endtime):
	# We should have unix time as an input, so convert that to Year/Month/day
	# We'll assume that we'll only ever need 2 days of data at most.
	ext = 'detector.cal'
	startdir = datetime.datetime.utcfromtimestamp(int(starttime)-1).strftime('%Y/%m/')
	print(indir+startdir)
	try:
		inputlist = os.listdir(indir + startdir)
		# print(inputlist)
		filelist = [startdir+'/'+f for f in inputlist if ext in f]

	except:
		print('No temperature data found!')
		return []

	enddir = datetime.datetime.utcfromtimestamp(int(endtime)+1).strftime('%Y/%m/')
	# print(startdir)
	# print(enddir)
	if enddir != startdir:
		# print('Hi')
		# Also need to append these
		try:
			inputlist = os.listdir(indir + enddir)
			for file in inputlist:
				if file not in filelist:
					filelist.append(enddir+'/'+file)
		except:
			print('No temperature data found for the second day!')
	print(filelist)
	timeset = []
	valset = []
	# We can just read in all the data and return it, since there isn't that much
	for file in filelist:
		data = np.loadtxt(indir+file,unpack=False,dtype={'names':('date','unix','t1','t2','t3'),'formats':('S1','d','f','f','f')})
		for line in data:
			timeset.append(line['unix'])
			valset.append([line['t1'],line['t2'],line['t3']])

	return timeset, valset


# From rhea_comm.reader_tod.header_read()
def read_rhea_header(fname):
	time = 0
	data = 0
	test = packet_reader.read_file(fname, length = 1)
	print(test)
	for val in test:
		print(val)
	# for t, n_r, s_off, d in packet_reader.read_file(fname, length = 1):
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

def read_kidslist(filename):
	with open(filename) as f:
		kiddict = json.load(f)
	# Rescale frequencies so they are all in Hz
	kiddict['kids_freqs'] = np.asarray(kiddict['kids_freqs']) * 1e6
	if 'blinds_freqs' in kiddict.keys():
		kiddict['blinds_freqs'] = np.asarray(kiddict['blinds_freqs']) * 1e6
	elif 'blinds_relfreqs' in kiddict.keys():
		kiddict['blinds_freqs'] = ((np.asarray(kiddict['kids_freqs'])/1e6) + np.asarray(kiddict['blinds_relfreqs'])) * 1e6
	return kiddict

def read_rhea_swp_data(fname, length=None, offset=0):
	inputdata = packet_reader.read_file(fname, length = length, offset = offset, sync=True)
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
	dataset.append(chan_accumulation / 10)
	del inputdata
	del chan_accumulation
	gc.collect()
	# We need one last dataset appending.
	return (np.asarray(freqset,dtype=np.float64), np.asarray(dataset,dtype=np.float64))

# This was running but really slow
# def read_rhea_data(filename, length=None, offset=0):
# 	buffsize = 4096
# 	header = 0xff
# 	header_sgsync = 0xaa
# 	header_sync = 0xf5
# 	footer = 0xee
# 	packet_size = packet_reader.get_packet_size(filename)
# 	buff = b''
# 	cnt = 0
# 	readcnt = -1
# 	step = 1
# 	f = open(filename, 'rb')
# 	## HEADER
# 	buff += f.read(buffsize*2)
# 	inputdata = np.asarray([packet_reader.read_iq_packet(buff[0:packet_size])])
#
# 	## BODY
# 	buff = b''
# 	n_rot,sync_off,offset = packet_reader.seek_sync(f,packet_reader.read_iq_packet,packet_size,offset)
# 	f.seek(packet_size * offset)
# 	while True:
# 		if type(length) == int and cnt >= length: break
# 		if len(buff) < packet_size: buff += f.read(buffsize)
# 		if len(buff) < packet_size: break
# 		if buff[0] == header_sync:
# 			n_rot, sync_off = packet_reader.read_sync_packet(buff[0:packet_size])
# 		else:
# 			readcnt += 1
# 			if readcnt % step ==0:
# 				cnt += 1
# 				inputdata = np.concatenate((inputdata, [packet_reader.read_iq_packet(buff[0:packet_size], n_rot, sync_off)]),axis=1)
# 		buff = buff[packet_size:]
#
# 	return inputdata

# Currently the best solution
def read_rhea_data(fname, length=None, offset=0):
	inputdata = packet_reader.read_file(fname, length = length, offset = offset, sync=True)#,packet_size=735)
	return np.array([x for x in inputdata],dtype=np.float64)
# This also works - but runs longer than the above.
	# dataset = []
	# for datum in inputdata:
	# 	print(datum)
	# 	exit()
	# 	datum[1].insert(0,datum[0])
	# 	datum[1].append(datum[2])
	# 	datum[1].append(datum[3])
	# 	dataset.append(datum[1])
	# 	del datum
	# inputdata.close()
	# del inputdata
	# gc.collect()
	# return np.asarray(dataset,dtype=np.float64)

def read_rhea_data_new(filename, length=None, offset=0):
	buffsize = 4096
	header = 0xff
	header_sgsync = 0xaa
	header_sync = 0xf5
	footer = 0xee
	packet_size = packet_reader.get_packet_size(filename)
	buff = b''
	cnt = 0
	readcnt = -1
	step = 1
	struct_b = struct.Struct("b")
	struct_h = struct.Struct("<H")
	struct_i = struct.Struct(">I")
	f = open(filename, 'rb')
	## HEADER
	buff += f.read(buffsize*2)
	inputdata = np.asarray([packet_reader.read_iq_packet(buff[0:packet_size])])

	## BODY
	buff = b''
	n_rot,sync_off,offset = packet_reader.seek_sync(f,packet_reader.read_iq_packet,packet_size,offset)
	f.seek(packet_size * offset)
	bytes_read = f.read()
	rows = list(bytes_read[i:i+packet_size] for i in range(0, len(bytes_read), packet_size))
	print(len(bytes_read))
	numpackets = len(bytes_read)/packet_size
	numvals = (packet_size - 7) / 7

	# Get the first column
	firstcol = [x[0] for x in rows]
	# Get the list of sync packets - use list comprehension
	sync = [n for n,x in enumerate(firstcol) if x==header_sync]
	numsync = len(sync)
	print(numsync)
	numdata = numpackets-numsync
	print(numdata)
	# Define a mask of whether they are data (True) or sync signals (False)
	mask = [False if x==header_sync else True for x in firstcol]
	print(len(mask))
	# print(mask)

	# Get the sync values
	justsync = list(compress(rows, [not x for x in mask]))
	n_rot = [(struct_b.unpack(x[1:2])[0] << 32) + struct_i.unpack(x[2:6])[0] for x in justsync]
	sync_off = [(struct_b.unpack(x[6:7])[0] << 48) + (struct_h.unpack(x[7:9])[0] << 32) + struct_i.unpack(x[9:13])[0] for x in justsync]
	del justsync

	# Now look at the data
	justdata = list(compress(rows, mask))
	data_out = np.zeros((int(numdata), int(numvals)+3))
	# Time
	data_out[:,0] = [(struct_b.unpack(x[1:2])[0] << 32) + struct_i.unpack(x[2:6])[0] for x in justdata]
	# Copy over the sync signals
	for i in range(0,len(sync)-1):
		data_out[sync[i]:sync[i+1],-2] = n_rot[i]
		data_out[sync[i]:sync[i+1],-1] = sync_off[i]
	print(numvals)
	# KIDs data
	for i in range(int(numvals)):
		print(i)
		firstcol = [x[6+7*i:7+7*i] for x in rows]
		secondcol = [x[7+7*i:9+7*i] for x in rows]
		thirdcol = [x[9+7*i:13+7*i] for x in rows]
		# firstcol = struct_b.unpack(firstcol << 48)[0]
		# secondcol = struct_h.unpack(firstcol << 32)[0]
		# thirdcol = struct_i.unpack(firstcol)[0]
		#
		# data_out[:,i+1] = firstcol + secondcol + thirdcol

		data_out[:,i+1] = [(struct_b.unpack(x[6+7*i:7+7*i])[0] << 48) + (struct_h.unpack(x[7+7*i:9+7*i])[0] << 32) + struct_i.unpack(x[9+7*i:13+7*i])[0] for x in justdata]

	return data_out
