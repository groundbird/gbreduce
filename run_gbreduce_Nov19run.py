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

#rsync -avz -e "ssh -A mpeel@moa.gb.iac.es ssh"  mpeel@dodo.gb.iac.es:/data/gb/ /net/nas/proyectos/cosmology/groundbird/data


# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# Set this to where you have a copy of the data
# basedir = '/Users/mpeel/Documents/GroundBIRD/data/'
# basedir = '/net/nas/proyectos/groundbird/data/'
basedir = '/Volumes/iOmega/GroundBIRD/data/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
# outdir = '/Users/mpeel/Documents/git/gbreduce/output/'
#outdir = '/Volumes/iOmega/GroundBIRD/analysis/'
outdir=basedir+'../analysis/'
# Set this to the location of the Azimuth encoder data
# azdir = '/Volumes/proyectos/cosmology/groundbird/data/logdata/rotary/'
azdir = basedir+'logbdata/az_enc/'
# Set this to the location of the elevation encoder data (in compressed format if elcompressed=True)
# eldir = '/Users/mpeel/Documents/groundbird/data/elevation/'
eldir = basedir+'logbdata/el_enc/'
elcompressed=True

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

toreduce = []
folder='kiddata/'
# Dome closed during SWP
# toreduce.append(['20191106/data_050316',datetime(2019,11,6,5,8,12).timestamp()])
# toreduce.append(['20191106/data_051016',datetime(2019,11,6,5,10,38).timestamp()])
# toreduce.append(['20191106/data_061040',datetime(2019,11,6,6,11,1).timestamp()])
# toreduce.append(['20191106/data_071103',datetime(2019,11,6,7,11,25).timestamp()])
# toreduce.append(['20191106/data_081127',datetime(2019,11,6,8,11,49).timestamp()])
# toreduce.append(['20191106/data_091151',datetime(2019,11,6,9,11,51).timestamp()])
# # Next one was a skydip
# # toreduce.append(['20191106/data_100028',datetime(2019,11,6,10,0,28).timestamp()])
# toreduce.append(['20191106/data_102531',datetime(2019,11,6,10,25,53).timestamp()])
# toreduce.append(['20191106/data_112555',datetime(2019,11,6,11,26,16).timestamp()])
# # Regeneration
# # toreduce.append(['20191106/data_122618',datetime(2019,11,6,12,26,39).timestamp()])
# # toreduce.append(['20191106/data_132641',datetime(2019,11,6,13,27,7).timestamp()])
# # toreduce.append(['20191106/data_142709',datetime(2019,11,6,14,27,31).timestamp()])
# # Dome closed during SWP
# # toreduce.append(['20191106/data_152733',datetime(2019,11,6,15,27,56).timestamp()])
# toreduce.append(['20191106/data_162758',datetime(2019,11,6,16,28,18).timestamp()])
# toreduce.append(['20191106/data_172820',datetime(2019,11,6,17,28,42).timestamp()])
# toreduce.append(['20191106/data_182844',datetime(2019,11,6,18,29,5).timestamp()])
# toreduce.append(['20191106/data_192907',datetime(2019,11,6,19,29,29).timestamp()])
# toreduce.append(['20191106/data_202931',datetime(2019,11,6,20,29,53).timestamp()])
# # Dome closed during SWP
# # toreduce.append(['20191106/data_212955',datetime(2019,11,6,21,30,16).timestamp()])
# toreduce.append(['20191106/data_223018',datetime(2019,11,6,22,30,29).timestamp()])
# toreduce.append(['20191106/data_233041',datetime(2019,11,6,23,31,3).timestamp()])

# # Dome closed during SWP
# #toreduce.append(['20191107/data_003105',datetime(2019,11,7,0,31,27).timestamp()])
# toreduce.append(['20191107/data_013129',datetime(2019,11,7,1,31,51).timestamp()])
# toreduce.append(['20191107/data_023153',datetime(2019,11,7,2,32,15).timestamp()])
# toreduce.append(['20191107/data_033217',datetime(2019,11,7,3,32,38).timestamp()])
# toreduce.append(['20191107/data_043240',datetime(2019,11,7,4,32,58).timestamp()])
# toreduce.append(['20191107/data_053300',datetime(2019,11,7,5,33,18).timestamp()])
# toreduce.append(['20191107/data_063320',datetime(2019,11,7,6,33,42).timestamp()])
# toreduce.append(['20191107/data_073344',datetime(2019,11,7,7,34,5).timestamp()])
# # Regeneration
# # toreduce.append(['20191107/data_083407',datetime(2019,11,7,8,34,31).timestamp()])
# # toreduce.append(['20191107/data_093434',datetime(2019,11,7,9,34,59).timestamp()])
# # toreduce.append(['20191107/data_103501',datetime(2019,11,7,10,35,23).timestamp()])
# toreduce.append(['20191107/data_113525',datetime(2019,11,7,11,35,46).timestamp()])
# toreduce.append(['20191107/data_123549',datetime(2019,11,7,12,36,10).timestamp()])
# toreduce.append(['20191107/data_133612',datetime(2019,11,7,13,36,34).timestamp()])
# toreduce.append(['20191107/data_143636',datetime(2019,11,7,14,36,36).timestamp()])
# toreduce.append(['20191107/data_153701',datetime(2019,11,7,15,37,22).timestamp()])
# toreduce.append(['20191107/data_163724',datetime(2019,11,7,16,37,36).timestamp()])
# toreduce.append(['20191107/data_171818',datetime(2019,11,7,17,18,18).timestamp()])

# toreduce.append(['20191111/data_121505',datetime(2019,11,11,12,15,20).timestamp()])
# toreduce.append(['20191111/data_124031',datetime(2019,11,11,12,40,43).timestamp()])
# toreduce.append(['20191111/data_134055',datetime(2019,11,11,13,41,16).timestamp()])

# toreduce.append(['20191111/data_144118',datetime(2019,11,11,14,41,40).timestamp()])
# toreduce.append(['20191111/data_154142',datetime(2019,11,11,15,42,3).timestamp()])

# toreduce.append(['20191111/data_164205',datetime(2019,11,11,16,42,27).timestamp()])
# toreduce.append(['20191111/data_174229',datetime(2019,11,11,17,42,50).timestamp()])
# toreduce.append(['20191111/data_184252',datetime(2019,11,11,18,43,13).timestamp()])
# toreduce.append(['20191111/data_194315',datetime(2019,11,11,19,43,37).timestamp()])
# toreduce.append(['20191111/data_204339',datetime(2019,11,11,20,44,1).timestamp()])
# toreduce.append(['20191111/data_214403',datetime(2019,11,11,21,44,24).timestamp()])
# toreduce.append(['20191111/data_224426',datetime(2019,11,11,22,44,48).timestamp()])
# toreduce.append(['20191111/data_234450',datetime(2019,11,11,23,45,12).timestamp()])

toreduce.append(['20191112/data_004514',datetime(2019,11,12,0,45,35).timestamp()])
toreduce.append(['20191112/data_014537',datetime(2019,11,12,1,45,59).timestamp()])
toreduce.append(['20191112/data_024601',datetime(2019,11,12,2,46,23).timestamp()])
toreduce.append(['20191112/data_034625',datetime(2019,11,12,3,46,47).timestamp()])
toreduce.append(['20191112/data_044649',datetime(2019,11,12,4,47,10).timestamp()])
toreduce.append(['20191112/data_054712',datetime(2019,11,12,5,47,34).timestamp()])
toreduce.append(['20191112/data_064736',datetime(2019,11,12,6,47,59).timestamp()])
toreduce.append(['20191112/data_074801',datetime(2019,11,12,7,48,20).timestamp()])

skymaps = []
hitmaps = []
for item in toreduce:
	skymaps.append(outdir+item[0].replace('/','_')+'/skymap_1.fits')
	hitmaps.append(outdir+item[0].replace('/','_')+'/hitmap_1.fits')
print(skymaps)
print(hitmaps)
prefix='combine'
outputname='201911_longrun'
# run.combine_sky_maps(skymaps,hitmaps,prefix,outputname,centralpos=(0,0),plotlimit=0.0)
# exit()





# # folder = 'kiddata/20191018/data_022435/'
# # starttime = datetime(2019,10,18,2,24,35).timestamp()
# folder = 'kiddata/20191018/data_031254/'
# name='20191018_data_031254'

# folder = 'kiddata/20191018/data_043458/'
# starttime = datetime(2019,10,18,4,35,43,tzinfo=pytz.timezone('UTC')).timestamp()
# name='20191018_data_043458'
for item in toreduce:
	folder = 'kiddata/'+item[0]+'/'
	starttime = item[1]
	name=item[0].replace('/','_')+'_swp'
	swpfile = 'swp.rawdata'
	todfile = 'tod.rawdata'
	kidparams = 'kids.list'
	print(starttime)
	# swp_params = run.analyse_swp(name,basedir+folder+swpfile,kidparams=basedir+folder+kidparams)
	swp_params=[]
	tod_analysis = run.analyse_tod(name,basedir+folder+todfile,kidparams=basedir+folder+kidparams,swp_params=swp_params,starttime=starttime)

