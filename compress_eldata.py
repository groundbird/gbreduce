#!/usr/bin/env python3
import eldata
import glob
from enum import Enum
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

class DataType(Enum):
    DATA = 1
    SYNC = 2
    UART = 3


def getdata(indir, filename):
	tmped = eldata.ElData(filename)
	numentries = tmped._length
	# numentries = 100

	# print(tmped.parse_all())
	times = []
	vals = []
	for i in range(0,numentries):
		try:
			times.append(tmped.get_data(i)[0])
			vals.append(tmped.get_data(i)[1])
		except:
			# If we find a problem, move onto the next one
			continue

	return times, vals


def getdata_condensed(indir, filename):
	tmped = eldata.ElData(filename)
	numentries = tmped._length
	# numentries = 100

	# print(tmped.parse_all())
	times = []
	vals = []
	lastval = 0
	for i in range(0,numentries):
		try:
			data = tmped.get_data(i)
			if data[2].name == 'DATA':
				if data[1] != lastval:
					times.append(data[0])
					vals.append(data[1])
					lastval = data[1]
		except:
			# If we find a problem, move onto the next one
			continue

	return times, vals

def compress_and_plot_eldata(indir,out_prefix,reload=False):
	if reload:
		times, vals, zenithdistance = np.loadtxt(indir+out_prefix+'.txt',unpack=True)
	else:
		files = sorted(Path(indir).rglob('*.dat'))
		print(files)
		times = []
		vals = []
		for filename in files:
			print(filename)
			timeset, valset = getdata_condensed(indir, filename)
			times = times + timeset
			vals = vals + valset
		# print(times)
		# print(vals)
		# print(stats.mode(vals,axis=None))
		# mode = stats.mode(vals,axis=None)
		# mode = int(mode[0])
		# print(len(vals))
		# print(mode)
		# print(vals[vals != mode])
		# print(len(vals[vals != mode]))
		zenithdistance = (np.asarray(vals).copy()-9113.0)/900.0
		np.savetxt(indir+out_prefix+'.txt',np.transpose([times, vals, zenithdistance]))
	times = np.asarray(times)
	vals = np.asarray(vals)

	plt.figure(figsize=(8.0, 5.0), dpi=100)
	# print(vals)
	print(min(times[times > 100.0]))
	# plt.plot(times[vals != mode], vals[vals != mode])
	plt.plot(times[times > 100.0], vals[times >= 100.0],'b-')
	plt.savefig(indir+out_prefix+'_encoder.png', dpi=1000)
	plt.clf()
	# mask = np.ones(len(times))
	# mask[times < 3.5e8] = 0
	# mask[times > 3.52e8] = 0
	# plt.plot(times[mask == 1], vals[mask == 1],'b.')
	# plt.savefig(indir+'elevation_zoom.png', dpi=1000)
	# plt.clf()


	# startel = 0.0
	# endel = -20.0
	# elrange=startel-endel
	# startencode = vals[0]
	# print(startencode)
	# endencode = vals[-1]
	# print(endencode)
	# encode_range = startencode-endencode
	# print(encode_range)
	# encode_perdeg = encode_range/(startel-endel)
	# print(encode_perdeg)
	# print('Expected:' + str((81000*4)/360.0))

	plt.plot(times, zenithdistance,'b-')
	plt.grid(True)
	plt.xlabel('Unix time')
	plt.ylabel('Zenith distance')
	plt.savefig(indir+'_zenith.png', dpi=1000)
	plt.clf()


	# plt.plot(times[mask==1], (vals[mask==1]-9113.0)/900.0,'b-')
	# plt.grid(True)
	# plt.xlabel('Unix time')
	# plt.ylabel('Zenith distance')
	# plt.savefig(indir+'elevation_corr_zoom.png', dpi=1000)
	# plt.clf()

	return
	