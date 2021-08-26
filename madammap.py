# Based on https://github.com/hpc4cmb/libmadam/blob/master/python/tests/test_libmadam.py
from mpi4py import MPI
import os
import sys
import shutil
import numpy as np
import healpy as hp
import libmadam_wrapper as madam


# Configuration
indirectory = '/home/mpeel/test/'
infile = indirectory+'20200110_data_001944_swp_poscor55_GB02_0_tod.fits'
nside = 64
npix = 12 * nside ** 2
fsample = 32.5
nnz = 1  # number or non zero pointing weights, typically 3 for IQU

pars = {}
pars["base_first"] = 1.0
pars["fsample"] = fsample
pars["nside_map"] = nside
pars["nside_cross"] = nside // 2
pars["nside_submap"] = nside // 4
pars["write_map"] = True
pars["write_binmap"] = True
pars["write_matrix"] = True
pars["write_wcov"] = True
pars["write_hits"] = True
pars["write_leakmatrix"] = True
pars["kfilter"] = True
pars["diagfilter"] = 0
pars["file_root"] = ""
pars["path_output"] = indirectory
pars["iter_max"] = 100
pars["nsubchunk"] = 2
pars["allreduce"] = True
pars["info"] = 0
pars["survey"] = ['test']
pars["bin_subsets"] = True

dets = ["GB020"]

# Running
comm = MPI.COMM_WORLD
itask = comm.Get_rank()
ntask = comm.Get_size()

hdul = fits.open(infile)
nsamp = len(hdul[0].data[0])

ndet = len(dets)
weights = np.ones(ndet, dtype=np.float64)
timestamps = np.zeros(nsamp, dtype=madam.TIMESTAMP_TYPE)
timestamps[:] = np.arange(nsamp) + itask * nsamp

# Convert from RA/dec to Healpix pixels
pixels = np.zeros(ndet * nsamp, dtype=madam.PIXEL_TYPE)
pixel[:] = hp.ang2pix(nside, (np.pi/2)-hdul[1].data.field(1)*np.pi/180.0, hdul[1].data.field(0)*np.pi/180.0)

pixweights = np.zeros(ndet * nsamp * nnz, dtype=madam.WEIGHT_TYPE)
pixweights[:] = 1

signal = np.zeros(ndet * nsamp, dtype=madam.SIGNAL_TYPE)
signal[:] = hdul[1].data.field(6)

nperiod = 1  # number of pointing periods

periods = np.zeros(nperiod, dtype=np.int64)

npsd = np.ones(ndet, dtype=np.int64)
npsdtot = np.sum(npsd)
psdstarts = np.zeros(npsdtot)
npsdbin = 10
psdfreqs = np.arange(npsdbin) * fsample / npsdbin
npsdval = npsdbin * npsdtot
psdvals = np.ones(npsdval)


comm.Reduce(hmap, hmap_tot, op=MPI.SUM, root=0)
comm.Reduce(bmap, bmap_tot, op=MPI.SUM, root=0)

madam.destripe(
	comm,
	pars,
	dets,
	weights,
	timestamps,
	pixels,
	pixweights,
	signal,
	periods,
	npsd,
	psdstarts,
	psdfreqs,
	psdvals,
)

# if itask == 0 and hp is not None:
# 	good = hmap_tot != 0
# 	bmap_tot[good] /= hmap_tot[good]

# 	hmap = hmap_tot.astype(np.int32)
# 	bmap = bmap_tot.astype(np.float32)

# 	# hp.write_map("hits.fits", hmap, nest=True, overwrite=True)
# 	# hp.write_map("binned.fits", bmap, nest=True, overwrite=True)

# 	madam_hmap = hp.read_map("pymaps/madam_pytest_hmap.fits", nest=True)
# 	madam_bmap = hp.read_map("pymaps/madam_pytest_bmap.fits", nest=True)

# 	npt.assert_allclose(madam_hmap, hmap)
# 	npt.assert_allclose(madam_bmap, bmap)

# 	madam_hmap1 = hp.read_map("pymaps/madam_pytest_hmap_sub1of2.fits", nest=True)
# 	madam_hmap2 = hp.read_map("pymaps/madam_pytest_hmap_sub2of2.fits", nest=True)

# 	npt.assert_allclose(madam_hmap, madam_hmap1 + madam_hmap2)

print("Done")
