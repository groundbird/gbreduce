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

# print(read_rhea_header(testfile))
# print(read_rhea_data(testfile))

run.analyse_tod('20190921_sron_kid000_04_rot_open_test1a',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test1/tod_-082.664MHz_-080.664MHz_+016.174MHz_+018.174MHz_+029.302MHz_+031.302MHz_+078.873MHz_+080.873MHz_0001kSPS_2019-0921-052531.rawdata',el=70,rpm=2.045,starttime=2458747.5)
run.analyse_tod('20190921_sron_kid000_04_rot_open_test1b',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test1/tod_-082.664MHz_-080.664MHz_+016.174MHz_+018.174MHz_+029.302MHz_+031.302MHz_+078.873MHz_+080.873MHz_0001kSPS_2019-0921-053018.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test2',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test2/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-054509.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test3',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test3/tod_-082.669MHz_-080.669MHz_+016.173MHz_+018.173MHz_+029.299MHz_+031.299MHz_+078.865MHz_+080.865MHz_0001kSPS_2019-0921-064155.rawdata',el=70,rpm=2.045,starttime=2458747.5)

run.analyse_tod('20190921_sron_kid000_04_rot_open_test4',basedir+'20190921/sron/sg_4990MHz/KID00_04/rot_open/test3/tod_-082.671MHz_-080.671MHz_+016.172MHz_+018.172MHz_+029.298MHz_+031.298MHz_+078.860MHz_+080.860MHz_0001kSPS_2019-0921-083503.rawdata',el=70,rpm=2.045,starttime=2458747.5)


exit()