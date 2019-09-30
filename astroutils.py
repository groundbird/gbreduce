import os
import numpy as np
import healpy as hp
from math import pi
import matplotlib.pyplot as plt

# Ensure that a directory exists
# Function from http://stackoverflow.com/questions/273192/in-python-check-if-a-directory-exists-and-create-it-if-necessary
# Upgraded to python v3, 14-May-2019
def ensure_dir(f):
	os.makedirs(f, exist_ok=True)
#def ensure_dir(f):
#	if f == '':
#		print('No directory path passed to ensure_dir, not checking or creating one!')
#		return
#	else:
#		d = os.path.dirname(f)
#		if not os.path.exists(d):
#			os.makedirs(d)


# Plot the maps
def plotmap(mapname,outputname):
	mapdata_temp = hp.read_map(mapname)
	fig = plt.figure(1)
	hp.mollview(mapdata_temp,fig=1)
	plt.savefig(outputname)
	plt.close()

# Create a healpix mask given an nside, min/max longitudes and latitudes, and a coordinate system
def healpixmask(nside, long_min, long_max, lat_min, lat_max, coordsystem='G'):

	npix = hp.nside2npix(nside)
	masked_map = np.zeros(npix)

	for i in range(0,npix):
		pos = hp.pixelfunc.pix2ang(nside, i)
		phi = pos[1]*(180.0/pi)
		theta = 90.0-(pos[0]*180.0/pi)

		# Check for negative longitude ranges
		if (phi <= long_max and phi >= long_min and theta <= lat_max and theta >= lat_min):
			masked_map[i] = 1
		if (long_max < 0): 
			# Assuming that this means that long_min is also less than zero
			if (phi-360.0 <= long_max and phi-360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 1
		if (long_min < 0): 
			# Longitudes higher than 0 will have been dealt with above, deal with the negative ones here.
			if (phi-360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 1
		if (long_max > 360): 
			# Assuming that this means that long_min is also less than zero
			if (phi+360.0 <= long_max and phi+360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 1
		if (long_min > 360): 
			# Longitudes higher than 0 will have been dealt with above, deal with the negative ones here.
			if (phi+360.0 >= long_min and theta <= lat_max and theta >= lat_min):
				masked_map[i] = 1


	return masked_map

# Create a mask of the Galactic plane. You can do this using healpixmask, but this is a bit more condensed.
# Original version was galacticmask.pro
def galacticmask(nside, degrees):
	npix = hp.nside2npix(nside)
	mask = np.ones(npix)
	for i in range(0,npix):
		pos = hp.pixelfunc.pix2ang(nside, i)
		if (abs(90.0-(pos[0]*180.0/pi)) <= degrees):
			mask[i] = 0
	return mask

# Compute factor to convert from brightness temp to thermodynamic temp
# INPUTS:
#   nu_ghz -    frequency in GHz
#
# OUTPUTS:
#   <value> -   conversion factor
#
# COMMENTS:
#   This conversion factor corresponds to the PLNCKCOR FITS header
#   keyword found in COBE/DMR data files produced by NASA.
#   For comparison, their results for 31.5, 53, and 90 GHz are
#
#   PLNCKCOR=             1.025724 /  Thermodynamic temperature = 
#   PLNCKCOR=             1.074197 /  Thermodynamic temperature =  
#   PLNCKCOR=             1.225941 /  Thermodynamic temperature =
#   COMMENT                        /   PLNCKCOR * antenna temperature  
#
# REVISION HISTORY:
#   Written by D. Finkbeiner, 10 March, 1999 - Berkeley
#   Converted to Python by M. Peel, 14 May 2019
def planckcorr(nu_ghz):
  k_b = 1.3806E-23              # J/K
  h = 6.6262E-34                # J*s
  T_cmb = 2.725                  # K	
  
  nu = nu_ghz*1.0E9             # Hz
  x = h * nu / (k_b * T_cmb)
  
  result = (np.exp(x)-1.)**2. / (x**2. * np.exp(x)) # thermodynamic uK  

  return result

# From http://home.fnal.gov/~stoughto/build/ARCONS/html/_modules/util/robust_sigma.html
def robust_sigma(in_y, zero=0):
   """
   Calculate a resistant estimate of the dispersion of
   a distribution. For an uncontaminated distribution,
   this is identical to the standard deviation.

   Use the median absolute deviation as the initial
   estimate, then weight points using Tukey Biweight.
   See, for example, Understanding Robust and
   Exploratory Data Analysis, by Hoaglin, Mosteller
   and Tukey, John Wiley and Sons, 1983.

   .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.

   :History:
       * H Freudenreich, STX, 8/90
       * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
       * Converted to Python by P. L. Lim, 11/2009

   Examples
   --------
   >>> result = robust_sigma(in_y, zero=1)

   Parameters
   ----------
   in_y: array_like
       Vector of quantity for which the dispersion is
       to be calculated

   zero: int
       If set, the dispersion is calculated w.r.t. 0.0
       rather than the central value of the vector. If
       Y is a vector of residuals, this should be set.

   Returns
   -------
   out_val: float
       Dispersion value. If failed, returns -1.

   """
   # Flatten array
   y = in_y.reshape(in_y.size, )

   eps = 1.0E-20
   c1 = 0.6745
   c2 = 0.80
   c3 = 6.0
   c4 = 5.0
   c_err = -1.0
   min_points = 3

   if zero:
       y0 = 0.0
   else:
       y0 = np.median(y)

   dy    = y - y0
   del_y = abs( dy )

   # First, the median absolute deviation MAD about the median:

   mad = np.median( del_y ) / c1

   # If the MAD=0, try the MEAN absolute deviation:
   if mad < eps:
       mad = np.mean( del_y ) / c2
   if mad < eps:
       return 0.0

   # Now the biweighted value:
   u  = dy / (c3 * mad)
   uu = u*u
   q  = np.where(uu <= 1.0)
   count = len(q[0])
   if count < min_points:
       print('ROBUST_SIGMA: This distribution is TOO WEIRD! Returning', c_err)
       return c_err

   numerator = np.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
   n    = y.size
   den1 = np.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
   siggma = n * numerator / ( den1 * (den1 - 1.0) )

   if siggma > 0:
       out_val = np.sqrt( siggma )
   else:
       out_val = 0.0

   return out_val


# Extract the pixels in an ellipse of given dimensions.
# Falls back to using query_disc when abratio=1.0.
# KNOWN BUG: does not work perfectly at high latitudes, and does not work when overlapping the poles.
# 
# Mike Peel, 19 March 2012 - started
# Mike peel, 14 May 2019 - converted to python
def query_ellipse(nside, lon, lat, radius, abratio, angle, nest=False):
	phi = lon*np.pi/180.0
	theta = np.pi/2.0-lat*np.pi/180.0
	radius_corr = radius * np.pi/180.0
	vec0 = hp.ang2vec(theta, phi)

	# Get all of the pixels that might be in the ellipse
	listpix_temp = hp.query_disc(nside, vec0, radius*np.pi/180.0,nest=nest)
	# Convert pixel to position
	thetanew, phinew = hp.pix2ang(nside, listpix_temp, nest=nest)

	if abratio != 1.0:
		# Translate and rotate coordinates
		x = (thetanew-theta) * np.cos(angle) + (phinew-phi) * np.sin(angle)*np.cos(lat*np.pi/180.)
		y = -(thetanew-theta) * np.sin(angle) + (phinew-phi) * np.cos(angle)*np.cos(lat*np.pi/180.)
		# check whether they lie within the ellipse
		listpix = listpix_temp[(x**2)/((radius_corr*abratio)**2) + (y**2)/(radius_corr**2) < 1.0]
	else:
		listpix = listpix_temp

	return listpix


# Do aperture photometry on a Healpix map with a circle of a given
# radius, and subtracting the background in an annulus
# Units are assumed to be mK_RJ and are converted to Jy
# 
# INPUTS
#
# inmap               Healpix fits file or map array 
#                     If an array, it assumed to be 'RING' ordering
#                     if ordering keyword not set 
#                  
# freq                Frequency (GHz) to allow conversion to Jy
#
# res_arcmin          Angular resolution (arcmin FWHM) of the data
#                     This is needed for the error estimation
#                     and only needs to be approximate 
#                     ( for noise_model=1 only )
#
# lon                 Longitude of aperture (deg)
#
# lat                 Latitude of aperture (deg)
#
# aper_inner_radius   Inner radius of aperture (arcmin)
#
# aper_outer_radius1  1st Outer radius of aperture beteween aperture
#                     and  b/g annulus (arcmin)
#                     Make this the same as the inner radius for
#                     single annulus with no gap
#
# aper_outer_radius2  2nd Outer radius of aperture for b/g annulus (arcmin)
#
# units               String defining units in the map
#                     Options include ['K','K_RJ', 'K_CMB', 'MJy/sr','Jy/pixel'] 
#                     m for milli and u for micro can also be used for
#                     K e.g. 'mK_RJ' or 'mK'
#                     Default is 'K_RJ'
#
# OPTIONAL:-
#
# column              This can be set to any integer which represents
#                     the column of the map (default is column=0)
#          
# dopol               If this keyword is set, then it will calculate the
#                     polarized intensity from columns 1 and 2 (Q and
#                     U) as PI=sqrt(Q^2+U^2) with *no noise bias* correction          
#                     N.B. This overrides the column option
#
# nested              If set, then the ordering is NESTED 
#                     (default is to assume RING if not reading in a file)
#
# noise_model         Noise model for estimating the uncertainty
#                     0 = approx uncertainty for typical/bg annulus aperture
#                     sizes (only approximate!).
#                     1 (DEFAULT) = assumes white uncorrelated noise (exact)
#                     and may under-estimate in some cases with real backgrounds!
#
# abratio             To do elliptical photometry, set this to the ratio of minor to major axis.
#                     Will be applied to all aperture radii.
#
# angle               To rotate the ellipse, set the angle here. 0 = elongated in direction of glat.
# 
# silent              If set, minimal output will appear
#
# Returns:
# fd                  Integrated flux density in aperture after b/g
#                     subtraction (Jy)
#
# fd_err              Estimate of error on integrated flux density (Jy)
#
# fd_bg               Background flux density estimate (Jy)
#
# HISTORY
#
# 26-Jun-2010  C. Dickinson   1st go
# 25-Jul-2010  C. Dickinson   Added conversion option from T_CMB->T_RJ
# 25-Jul-2010  C. Dickinson   Added 2nd outer radius
# 26-Aug-2010  C. Dickinson   Added generic 'unit' option
# 02-Sep-2010  C. Dickinson   Tidied up the code a little
# 19-Oct-2010  C. Dickinson   Use ang2vec rather than long way around
#                             via the pixel number
# 20-Oct-2010  M. Peel        Fix MJy/Sr conversion; add aliases for
#							  formats excluding underscores
# 10-Feb-2011  C. Dickinson   Added column option to allow polarization
# 16-Feb-2011  C. Dickinson   Added /dopol keyword for doing polarized intensity
# 12-Mar-2011  C. Dickinson   Added flexibility of reading in file or array
# 01-Jun-2011  C. Dickinson   Added noise_model option
# 10-Nov-2011  M. Peel        Adding AVG unit option.
# 02-Feb-2012  C. Dickinson   Added silent option
# 04-Mar-2012  C. Dickinson   Changed default noise model to 1
# 19-Mar-2012  M. Peel        Added support for elliptical annuli
# 06-Dec-2012  M. Peel		  Fix bug with noise_model always being set to 1
# 14-May-2019  M. Peel        Convert to Python
#
def haperflux(inmap, freq, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=0, abratio=1.0, angle=0.0, silent=False):

	# read in data
	if isinstance(inmap, str) and not silent:
		print ('***' + inmap)
		map = hp.read_map(inmap)
	else:
		map = inmap

	# Handle the case where we only have one map
	if len(map) > 100:
		map = [map]

	# Check dimensions
	nside = hp.get_nside(map[0])
	if not isinstance(nside,int):
		print('Not a standard Healpix map...')
		return
	npix = hp.nside2npix(nside)
	ncolumn = len(map)

	# set column number and test to make sure there is enough columns in the file!
	if column+1 > ncolumn and not dopol:
		print('Column number requested larger than the number of columns in the file!')
		return

	# check for dopol keyword for calculating polarized intensity
	if dopol:
		column=1
		if ncolumn < 3:
			print('To calculate polarized intensity (PI), requires polarization data with 3 columns or more...')
			return

	# get pixels in aperture
	innerpix = query_ellipse(nside, lon, lat, aper_inner_radius/60.0, abratio, angle,nest=nested)
	outerpix1 = query_ellipse(nside, lon, lat, aper_outer_radius1/60.0, abratio, angle,nest=nested)
	outerpix2 = query_ellipse(nside, lon, lat, aper_outer_radius2/60.0, abratio, angle,nest=nested)

	# do not include pixels that are bad
	innerpix = innerpix[map[column][innerpix] > -1.e30]
	innerpix = innerpix[map[column][innerpix] != 0.0000]
	outerpix1 = outerpix1[map[column][outerpix1] > -1.e30]
	outerpix1 = outerpix1[map[column][outerpix1] != 0.0000]
	outerpix2 = outerpix2[map[column][outerpix2] > -1.e30]
	outerpix2 = outerpix2[map[column][outerpix2] != 0.0000]

	ninnerpix = len(innerpix)
	nouterpix1 = len(outerpix1)
	nouterpix2 = len(outerpix2)
	# print(ninnerpix,nouterpix1,nouterpix2)

	if ninnerpix == 0 or nouterpix1 == 0 or nouterpix2 == 0:
		print, '*** No good pixels inside aperture! ***'
		fd = float('nan')
		fd_err = float('nan')
		fd_bg = float('nan')
		return fd,fd_err,fd_bg

	# find pixels in the annulus (between outerradius1 and outeradius2) 
	outerpix = list(set(outerpix2)-set(outerpix1))
	nouterpix = len(outerpix)

	bg_zero = 0
	if len(outerpix) <= 1:
		print('***No good pixels in the background annulus!***')
		outerpix = -1
		nouterpix = 0
		fd_bg = 0.
		bg_zero = 1


	# get conversion from to Jy/pix
	pix_area = 4.0*np.pi / npix
	factor = 1.0
	if units == 'K' or units == 'K_RJ' or units == 'KRJ':
		factor = 2.*1381.*(freq*1.0e9)**2/(2.997e8)**2 * pix_area
	elif units == 'mK' or units == 'mK_RJ' or units == 'mKRJ':
		factor = 2.*1381.*(freq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e3
	elif units == 'uK' or units == 'uK_RJ' or units == 'uKRJ':
		factor = 2.*1381.*(freq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e6
	elif units == 'K_CMB' or units == 'KCMB':
		factor = 2.*1381.*(freq*1.0e9)**2/(2.997e8)**2 * pix_area / planckcorr(freq)
	elif units == 'mK_CMB' or units == 'mKCMB':
		factor = 2.*1381.*(freq*1.0e9)**2/(2.997e8)**2 * pix_area / 1.0e3 / planckcorr(freq)
	elif units == 'uK_CMB' or units == 'uKCMB':
		factor = 2.*1381.*(freq*1.0e9)^2/(2.997e8)^2 * pix_area / 1.0e6 / planckcorr(freq)
	elif units == 'MJy/sr' or units == 'MJY/SR' or units == 'MjySr':
		factor = pix_area * 1.0e6
	elif units == 'Jy/pixel' or units == 'JY/PIXEL' or units == 'JY/PIX' or units == 'JyPix':
		factor = 1.0
	elif units == 'average' or units == 'avg' or units == 'Average' or units == 'AVG':
		factor = 1.0 / float(ninnerpix)
	# print(factor)
	# ; override columns if /dopol keyword is set
	# IF (dopol EQ 1) THEN ncalc = 2 ELSE ncalc = 1

	# FOR i=1L, ncalc DO BEGIN

	# IF (dopol EQ 1) THEN column=i

	# Get pixel values in inner radius, converting to Jy/pix
	fd_jypix_inner = map[column][innerpix] * factor
	# sum up integrated flux in inner
	fd_jy_inner = np.sum(fd_jypix_inner)

	# subtract background
	if bg_zero != 1:
		# same for outer radius but take a robust estimate and scale by area
		fd_jy_outer = np.median(map[column][outerpix]) * factor * float(ninnerpix)
	else:
		fd_bg = 0.
		fd_jy_outer = 0.

	fd_bg = fd_jy_outer
	fd = fd_jy_inner - fd_bg
	fd_err = 0
	# ; estimate error based on robust sigma of the background annulus
	# IF (size(noise_model,/TYPE)) EQ 0 THEN noise_model = 1 ELSE noise_model = noise_model

	# ; new version (2-Dec-2010) that has been tested with simulations for
	# ; Planck early paper and seems to be reasonable for many applications
	if noise_model == 0:
		npoints = (pix_area*ninnerpix) /    (1.13*(float(res_arcmin)/60. * np.pi/180.)**2)
		npoints_outer = (pix_area*nouterpix) /    (1.13*(float(res_arcmin)/60. * np.pi/180.)**2)
		fd_err = np.std(map[column][outerpix])  * factor * ninnerpix / np.sqrt(npoints)

	# (works exactly for white uncorrelated noise only!)
	if noise_model == 1:
		k = np.pi/2.
		fd_err = factor * np.sqrt(float(ninnerpix) + (k * float(ninnerpix)**2/nouterpix)) * robust_sigma(map[column][outerpix])
	
	# ; if dopol is set, then store the Q estimate the first time only
	# IF (dopol EQ 1 AND i EQ 1) THEN BEGIN
	# fd1 = fd
	# fd_err1 = fd_err
	# fd_bg1 = fd_bg
	# ENDIF

	# ENDFOR


	# ; if dopol is set, combine Q and U to give PI
	# IF (dopol EQ 1) THEN BEGIN
	# fd = sqrt(fd1^2 +fd^2)
	# fd_err = sqrt( 1./(fd1^2 + fd^2) * (fd1^2*fd_err1^2 + fd^2*fd_err^2))
	# fd_bg = sqrt(fd_bg1^2 + fd_bg^2)
	# ENDIF


	# ; display results
	# IF (keyword_set(silent) EQ 0) THEN BEGIN
	# print, ''
	# print, 'Coordinates: ', lon, lat
	# IF (dopol EQ 0) THEN print, 'Column: ', column ELSE print, 'Column: ', 'Pol (sqrt(col1^2 + col2^2))'
	# print, 'Noise model: ', noise_model
	# IF (noise_model EQ 0) THEN print, 'N.B. errors are only approximate - user should test to see if they are reasonable for their purposes'
	# IF (noise_model EQ 1) THEN print, 'N.B.errors will only be strictly exact for white uncorrelated noise with no background'

	# ;print, 'Angular resolution (arcmin): ', res_arcmin
	# ;print, 'Inner: ', fd_jy_inner
	# ;print, 'B/g: ', fd_bg
	# ;print, 'Npixels (inner):', ninnerpix
	# ;print, 'Npixels (outer):', nouterpix
	# ;print, 'Npoints: ', npoints
	# print, 'Integrated flux: ', fd
	# print, 'Integrated flux error: ', fd_err
	# print, ''

	# ;print, '***', freq, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units

	# ENDIF
	print(freq,fd,fd_err,fd_bg)
	return fd, fd_err, fd_bg
