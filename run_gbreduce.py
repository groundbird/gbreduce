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

from rhea_comm.lib_read_rhea import *
from rhea_comm.packet_reader import read_file, get_length

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# Set this to where you have a copy of the data
basedir = '/Users/mpeel/Documents/GroundBIRD/data/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Users/mpeel/Documents/git/gbreduce/output/'



# Start the class
run = gbreduce.gbreduce(outdir=outdir,\
	datadir=basedir,\
	nside = 512)

# testfile = basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test2/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-054509.rawdata'

# print(read_rhea_swp(testfile))

# print(read_rhea_tod(testfile))

# swp_params = run.analyse_swp('20190921_sron_kid000_04_rot_open_test2',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test2/mulswp_+003.000MHzWidth_+000.001MHzStep_-082.740MHz_-080.740MHz_+016.162MHz_+018.162MHz_+029.276MHz_+031.276MHz_+078.775MHz_+080.775MHz_2019-0921-053919.rawdata',freqrange=3.0,freqstep=0.001,centerfreqs=[-082.740,-080.740,+016.162,+018.162,+029.276,+031.276,+078.775,+080.775])
# swp_params = []
# tod_analysis = run.analyse_tod('20190921_sron_kid000_04_rot_open_test2',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test2/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-054509.rawdata',el=70,rpm=2.045,starttime=2458747.5,swp_params=swp_params,centerfreqs=[-082.669,-080.669,+016.173,+018.173,+029.299,+031.299,+078.865,+080.865],lo=4.99e9)
swp_params = run.analyse_swp('20190924_el80last',basedir+'20190924/moon/El_80_last/data_mulch__2019-0924-102753/swpmul__el80_+006.000MHzWidth_+000.010MHzStep_-082.501MHz_0016.304MHz_0029.452MHz_0079.066MHz_2019-0924-102753.rawdata',freqrange=6.0,freqstep=0.01,centerfreqs=[-082.501,0016.304,0029.452,0079.066])
# swp_params = []
tod_analysis = run.analyse_tod('20190924_el80last',basedir+'20190924/moon/El_80_last/data_mulch__2019-0924-102753/tod__el80_+2.0MHzBlindTone_1kSPS_-082.497MHz_0016.303MHz_0029.451MHz_0079.074MHz_2019-0924-102825.rawdata',el=80,rpm=2.045,starttime=2458750.5,swp_params=swp_params,centerfreqs=[082.497,0016.303,0029.451,0079.074],lo=4.99e9,numpix=2)

exit()

# print(read_rhea_header(testfile))
# print(read_rhea_data(testfile))

run.analyse_tod('20190921_sron_kid000_04_rot_open_test1a',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test1/tod_-082.664MHz_-080.664MHz_+016.174MHz_+018.174MHz_+029.302MHz_+031.302MHz_+078.873MHz_+080.873MHz_0001kSPS_2019-0921-052531.rawdata',el=70,rpm=2.045,starttime=2458747.5)
exit()
run.analyse_tod('20190921_sron_kid000_04_rot_open_test1b',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test1/tod_-082.664MHz_-080.664MHz_+016.174MHz_+018.174MHz_+029.302MHz_+031.302MHz_+078.873MHz_+080.873MHz_0001kSPS_2019-0921-053018.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test2',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test2/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-054509.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test3',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test3/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-064155.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test4',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test3/tod_-082.671MHz_-080.671MHz_+016.172MHz_+018.172MHz_+029.298MHz_+031.298MHz_+078.860MHz_+080.860MHz_0001kSPS_2019-0921-083503.rawdata',el=70,rpm=2.045,starttime=2458747.5)


exit()