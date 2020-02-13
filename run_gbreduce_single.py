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

from rhea_comm.lib_read_rhea import *
from rhea_comm.packet_reader import read_file, get_length

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

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


# Run a batch job
subdir = 'kiddata/20200204/'
skipfirst = 0 # Note: this already excludes folders containing files with 'KSPS' in the filenames
run.runset(subdir=subdir,ext=ext,skipfirst=skipfirst)

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


# EOF