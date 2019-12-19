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
# basedir = '/Users/mpeel/Documents/GroundBIRD/data/'
# basedir = '/Volumes/proyectos/cosmology/groundbird/data/'
basedir = '/Volumes/iOmega/GroundBIRD/data/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Volumes/iOmega/GroundBIRD/analysis/'
# Set this to the location of the Azimuth encoder data
azdir = basedir+'logbdata/az_enc/'
# Set this to the location of the elevation encoder data (in compressed format if elcompressed=True)
eldir = basedir+'logbdata/el_enc/'
elcompressed=False
domedir = basedir+'logdata/dome/'
tempdir = basedir+'logdata/thermo/'

# Start the class
run = gbreduce.gbreduce(outdir=outdir,\
	datadir=basedir,\
	azdir=azdir,\
	eldir=eldir,\
	elcompressed=elcompressed,\
	domedir=domedir,\
	tempdir=tempdir,\
	nside = 512, use_mkidpylibs=True)

subdir = 'kiddata/20191216/'
ext = '_swp_newaz'
skipfirst = 0 # Note: this already excludes folders containing files with 'KSPS' in the filenames
run.runset(subdir=subdir,ext=ext,skipfirst=skipfirst)

