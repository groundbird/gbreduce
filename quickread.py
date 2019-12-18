import matplotlib.pyplot as plt
import numpy as np
from gbreduce_read import *
from gbreduce_functions import *
import gbreduce

indir = '/Users/mpeel/Desktop/temp/'
indir = '/Volumes/iOmega/GroundBIRD/data/kiddata/20191218/data_044428/'
# filenames = ['swp_16.rawdata','swp_32.rawdata']
filenames = ['tod_16.rawdata','tod_32.rawdata']
# kidfiles = ['kids_16.list','kids_32.list']
run = gbreduce.gbreduce(outdir=indir,\
	datadir=indir,\
	nside = 512, use_mkidpylibs=True)


for filename in filenames:
	if 'swp' in filename:
		# data = read_rhea_swp_data(indir+filename)
		# data = data[1][:][:]
		# print(np.shape(data))
		run.analyse_swp(filename.replace('.','_'),indir+filename,kidparams=indir+filename.replace('rawdata','list').replace('swp','kids'))
	else:
		data = read_rhea_data(indir+filename)
		print(np.shape(data))
		# data -= np.mean(data,axis=0)
		std = np.std(data[:,3])
		for pix in range(0,len(data[0])):
			data[:,pix] = subtractbaseline(data[:,pix],option=0,navg=100) + pix * std * 3
		data = data[100:-100,:]
		fig = plt.gcf()
		fig.set_size_inches(10, 8)
		plt.plot(data)
		plt.savefig(indir+filename+'.png')
		plt.clf()