import healpy as hp
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from gbreduce_functions import *
from gbreduce_read import *
import os
import matplotlib as mpl
from astropy.time import Time

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

outputdir='/Volumes/Toshiba5TB2/GroundBIRD/analysis/'
subdir='_maps/'

basedir = '/Volumes/Toshiba5TB2/GroundBIRD/data/'
domedir = basedir+'logdata/dome/'

folderlist = os.listdir(outputdir)
todolist = []

use_dome_info = False

# startdir = '20200204_data_171138_swp_poscor55_GB02'
startdir = '20191106_data_050316_swp_poscor55'
enddir = '20191112_data_151536_swp_poscor55'
trip = 0
if startdir == '':
	trip = 1

for folder in folderlist:
	if (trip or startdir in folder) and not '.DS' in folder:
		trip = 1
		subfolderlist = os.listdir(outputdir+folder)
		search = 'log.txt'
		search2 = '_tod.fits'
		test = [f for f in subfolderlist if search in f]
		test2 = [f for f in subfolderlist if search2 in f]
		if test != [] and test2 != []:
			todolist.append([outputdir+folder+'/', test2])
	if enddir in folder:
		break
print(todolist)

nside = 64
npix=12*nside*nside
frequencies = []
skymaps = []
hitmaps = []
flag = []
filesdone = 0
for file in todolist:
	filesetnum = 0
	for fitsfile in file[1]:
		filesetnum += 1
		print(file[0]+fitsfile)
		hdul = fits.open(file[0]+fitsfile)
		print(hdul)
		print(hdul[0].header['CENFREQ'])
		print(frequencies)
		matchfreq = int(float(hdul[0].header['CENFREQ'])/1e6)
		print(matchfreq)
		try:
			index = frequencies.index(matchfreq)
		except:
			try:
				index = frequencies.index(matchfreq-1)
			except:
				try:
					index = frequencies.index(matchfreq+1)
				except:
					print('Not found, adding new array')
					frequencies.append(matchfreq)
					skymaps.append([np.zeros(npix), np.zeros(npix), np.zeros(npix), np.zeros(npix)])
					hitmaps.append(np.zeros(npix))
					print(np.shape(skymaps))
					index = len(frequencies)-1
		print(hdul[1].data.field(0))

		# Convert to Healpix pixels
		healpix_pixel = hp.ang2pix(nside, (np.pi/2)-hdul[1].data.field(1)*np.pi/180.0, hdul[1].data.field(0)*np.pi/180.0)
		print(len(healpix_pixel))

		if len(np.unique(healpix_pixel)) > 10:

			# Check the dome
			if filesetnum == 1 and use_dome_info:
				starttime = Time(np.min(hdul[1].data.field(2)), format='mjd', scale='utc')
				endtime = Time(np.max(hdul[1].data.field(2)), format='mjd', scale='utc')

				domet, domef = fetch_domedata(domedir, starttime.unix, endtime.unix)
				domet = Time(domet, format='unix', scale='utc').mjd
				flag = np.zeros(len(healpix_pixel))
				# print(domet,domef)
				if len(domet) >= 1:
					# We can use the dome opening/closing to do an initial cut of the data
					dometimes_pos = 0
					for i in range(0,len(hdul[1].data.field(2))):
						# print(domet[dometimes_pos])
						# print(hdul[1].data.field(2)[i])
						while(domet[dometimes_pos] < hdul[1].data.field(2)[i] and dometimes_pos < len(domet)-1):
							dometimes_pos += 1
						if dometimes_pos == len(domet):
							if int(domef[dometimes_pos]) == 0:
								flag[i] = 1
						else:
							if int(domef[dometimes_pos-1]) == 0:
								flag[i] = 1
			if use_dome_info:
				healpix_pixel[flag == 1] = -1
				print('Dome flagged:')
				print(np.sum(healpix_pixel==-1))

			num_to_average = 200
			col1 = subtractbaseline(hdul[1].data.field(3),option=0,navg=num_to_average)
			col2 = subtractbaseline(hdul[1].data.field(4),option=0,navg=num_to_average)
			col3 = subtractbaseline(hdul[1].data.field(5),option=0,navg=num_to_average)
			col4 = subtractbaseline(hdul[1].data.field(6),option=0,navg=num_to_average)

			print(index)
			for pixnum in np.unique(healpix_pixel):
				if pixnum != -1:
					skymaps[index][0][pixnum] = skymaps[index][0][pixnum] + np.sum(col1[healpix_pixel==pixnum])
					skymaps[index][1][pixnum] = skymaps[index][1][pixnum] + np.sum(col2[healpix_pixel==pixnum])
					skymaps[index][2][pixnum] = skymaps[index][2][pixnum] + np.sum(col3[healpix_pixel==pixnum])
					skymaps[index][3][pixnum] = skymaps[index][3][pixnum] + np.sum(col4[healpix_pixel==pixnum])
					hitmaps[index][pixnum] = hitmaps[index][pixnum] + np.sum(healpix_pixel==pixnum)

		hdul.close()
	filesdone += 1
	# if filesdone > 10:
	# 	break

for index in range(0,len(frequencies)):
	if np.sum(hitmaps[index]) > 10.0:
		for i in range(0,npix):
			if hitmaps[index][i] >= 1:
				skymaps[index][0][i] = skymaps[index][0][i]#/hitmaps[index][i]
				skymaps[index][1][i] = skymaps[index][1][i]#/hitmaps[index][i]
				skymaps[index][2][i] = skymaps[index][2][i]#/hitmaps[index][i]
				skymaps[index][3][i] = skymaps[index][3][i]#/hitmaps[index][i]
			else:
				skymaps[index][0][i] = hp.pixelfunc.UNSEEN
				skymaps[index][1][i] = hp.pixelfunc.UNSEEN
				skymaps[index][2][i] = hp.pixelfunc.UNSEEN
				skymaps[index][3][i] = hp.pixelfunc.UNSEEN


		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_0.fits',skymaps[index][0],overwrite=True)
		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_1.fits',skymaps[index][1],overwrite=True)
		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_2.fits',skymaps[index][2],overwrite=True)
		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_3.fits',skymaps[index][3],overwrite=True)
		hp.write_map(outputdir+subdir+'hits_'+str(nside)+'_'+str(frequencies[index])+'.fits',hitmaps[index],overwrite=True)
		hp.mollview(skymaps[index][0])
		plt.savefig(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_0.png')
		hp.mollview(hitmaps[index])
		plt.savefig(outputdir+subdir+'hits_'+str(nside)+'_'+str(frequencies[index])+'.png')


