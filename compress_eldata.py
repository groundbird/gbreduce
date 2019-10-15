#!/usr/bin/env python3
import eldata
import glob
from enum import Enum
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def get_el_data(indir, filename):
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

def get_el_data_condensed(indir, filename):
	tmped = eldata.ElData(filename)
	numbytes, version, unixtime, headertxt = tmped.get_header()
	times = []
	vals = []
	lastval = 0
	starttime = 0
	for i in range(0,tmped._length):
		try:
			data = tmped.get_data(i)
			if data[2].name == 'DATA':
				if starttime == 0:
					starttime = data[0]
				if data[1] != lastval:
					times.append((data[0] - starttime)*0.001 + unixtime)
					vals.append(data[1])
					lastval = data[1]
		except:
			# If we find a problem, move onto the next one
			continue

	return times, vals

def compress_and_plot_eldata(indir,out_prefix,reload=False,savedata=True):
	if reload:
		times, vals, zenithdistance = np.loadtxt(indir+out_prefix+'.txt.gz',unpack=True)
	else:
		files = sorted(Path(indir).rglob('*.dat'))
		print(files)
		times = []
		vals = []
		for filename in files:
			print(filename)
			timeset, valset = get_el_data_condensed(indir, filename)
			times = times + timeset
			vals = vals + valset

		# Calculate zenith distance, and save to disk
		zenithdistance = (np.asarray(vals).copy()-9113.0)/900.0
		if savedata:
			np.savetxt(indir+out_prefix+'.txt.gz',np.transpose([times, vals, zenithdistance]))

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
	