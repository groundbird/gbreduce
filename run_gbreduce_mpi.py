#!/usr/bin/env python
# -*- coding: utf-8  -*-
# 
# This is an example run script for gbreduce

import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gbreduce
from gbreduce_read import *
import datetime
import pytz
from mpi4py import MPI

from rhea_comm.lib_read_rhea import *
from rhea_comm.packet_reader import read_file, get_length

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# This is needed to start the MPI run
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
status = MPI.Status()
# Constants
MASTER_PROCESS = 0
WORK_TAG = 1
DIE_TAG = 2
print("-"*78)
print(" Running on %d cores" % comm.size)
print(" Rank %d" % rank)
print("-"*78)

# Set this to where you have a copy of the data
basedir = '/Volumes/Toshiba5TB2/GroundBIRD/data/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Volumes/Toshiba5TB2/GroundBIRD/analysis/'
# Set this to the directory for KID data
datadir = basedir+'kiddata/'
# Set this to the location of the Azimuth encoder data
azdir = basedir+'logbdata/az_enc/'
# Set this to the location of the elevation encoder data (in compressed format if elcompressed=True)
eldir = basedir+'logbdata/el_enc/'
elcompressed=False
domedir = basedir+'logdata/dome/'
tempdir = basedir+'logdata/thermo/'
pixinfo = 'gb_pixinfo_20191217.txt'
nside = 512

# Settings for the run
ext = '_swp_poscor'

# Set this to the start and end directories, the code will automatically find and run the subdirectories between them.
startdir = ''#'20200103'
enddir = ''#'20200131'
skipfirst = 0 # Note: this already excludes folders containing files with 'KSPS' in the filenames


# You shouldn't need to change anything below!
# Start the class
run = gbreduce.gbreduce(outdir=outdir,\
	datadir=basedir,\
	azdir=azdir,\
	eldir=eldir,\
	elcompressed=elcompressed,\
	domedir=domedir,\
	tempdir=tempdir,\
	pixinfo=pixinfo,\
	nside=nside, use_mkidpylibs=True)


### Master Process ###
if rank == MASTER_PROCESS:
	num_processors = comm.size
	print("Master process found " + str(num_processors) + " worker processors.")

	# Work out the subfolders to process
	folderlist = os.listdir(datadir)
	work_array = []
	trip = 0
	count = 0
	if startdir == '':
		trip = 1
	for folder in sorted(folderlist):
		if trip or startdir in folder:
			work_array.append([count, folder])
			trip = 1
			count += 1
		if enddir != '':
			if enddir in folder:
				trip = 0

	print(work_array)
	print(len(work_array))
	work_size = len(work_array)
	# Dispatch jobs to worker processes
	work_index = 0
	num_completed = 0

	# Start all worker processes
	for i in range(1, min(num_processors, work_size+1)):
		comm.send(work_index, i, tag=WORK_TAG)
		comm.send(work_array[work_index], i)
		print("Sent work index " + str(work_index) + " to processor " + str(i))
		work_index += 1

	# Receive results from each worker, and send it new data
	for i in range(num_processors, work_size+1):
		results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		index = status.tag
		proc = status.source
		num_completed += 1
		comm.send(work_index, proc, tag=WORK_TAG)
		comm.send(work_array[work_index], proc)
		print("Sent work index " + str(work_index) + " to processor " + str(proc))
		work_index += 1

	# Get results from remaining worker processes
	while num_completed < work_size-1:
		results = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
		num_completed += 1

	# Shut down worker processes
	for proc in range(1, num_processors):
		print("Stopping worker process " + str(proc))
		comm.send(-1, proc, tag=DIE_TAG)

else:
	### Worker Processes ###
	continue_working = True
	while continue_working:

		work_index =  comm.recv(source=MASTER_PROCESS, tag=MPI.ANY_TAG, status=status)

		if status.tag == DIE_TAG:
			continue_working = False
		else:
			work_array = comm.recv(source=MASTER_PROCESS, tag=MPI.ANY_TAG, status=status)
			work_index = status.tag

			# We want to pause for a random time between 0 and 3 minutes
			# to avoid all jobs asking for data at the same time
			# time.sleep(random.random_integers(low=0, high=180))

			print("Running job " + str(work_array[0]) + " on process " + str(rank) + ", subfolder " + work_array[1])
			run.runset(subdir=datadir+work_array[1],ext=ext,skipfirst=skipfirst)

			comm.send(work_array, dest=MASTER_PROCESS, tag=work_index)


# times, vals, times_sync, vals_sync = get_el_data_condensed(eldir+'2020/01/13/el_2020-0113-005556+0000.dat.xz')
# print(vals)
# print(np.median(vals))
# exit()
# Run a batch job
# subdir = 'kiddata/20191111/'
# skipfirst = 0 # Note: this already excludes folders containing files with 'KSPS' in the filenames
# run.runset(subdir=subdir,ext=ext,skipfirst=skipfirst)

# # Or run a single job
# folder = 'kiddata/20200109/data_005248/'
# starttime = datetime(2019,12,17,6,25,34).timestamp()
# name='20201009_data_005248'
# swpfile = 'swp.rawdata'
# todfile = 'tod.rawdata'
# kidparams = 'kids.list'
# postfixes = ['_16', '_32']

# # swp_params = run.analyse_swp(name+postfix,basedir+folder+swpfile.replace('.',postfix+'.'),kidparams=basedir+folder+kidparams.replace('.',postfix+'.'))
# # swp_params = []
# tod_analysis = run.analyse_tod(name+postfix,basedir+folder+todfile.replace('.',postfix+'.'),kidparams=basedir+folder+kidparams.replace('.',postfix+'.'),swp_params=swp_params,starttime=starttime)

# Close
MPI.Finalize()

# EOF