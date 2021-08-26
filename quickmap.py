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
print(folderlist)
todolist = []

use_flag = True
moon_mindist = -1#5.0 # Set to -1 to include moon data
sun_mindist = -1#30.0 # Set to -1 to include sun data

# startdir = '20200204_data_171138_swp_poscor55_GB02'
startdir = ''#'20191106_data_050316_swp_poscor55'
enddir = ''#20191112_data_151536_swp_poscor55'
trip = False
if startdir == '':
	trip = True

print(trip)
for folder in folderlist:
	print(folder)
	if (trip or startdir in folder) and not '.' in folder:
		trip = True
		subfolderlist = os.listdir(outputdir+folder)
		print(subfolderlist)
		search = '_log.txt'
		search2 = '_tod.fits'
		test = [f for f in subfolderlist if search in f]
		test2 = [f for f in subfolderlist if search2 in f]
		if test != [] and test2 != []:
			todolist.append([outputdir+folder+'/', test2])
	if enddir != '' and enddir in folder:
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
		print(hdul[0])
		# print(hdul[0].header)
		# exit()
		# print(hdul[0].header['CENFREQ'])
		# print(frequencies)
		# matchfreq = int(float(hdul[0].header['CENFREQ'])/1e6)
		if 'GB03' in file[0]+fitsfile:
			arr = 'GB03'
		else:
			arr = 'GB04'
		pixname = arr+str(hdul[0].header['PIX'])
		# try:
		# 	index = frequencies.index(matchfreq)
		# except:
		# 	try:
		# 		index = frequencies.index(matchfreq-1)
		# 	except:
		# 		try:
		# 			index = frequencies.index(matchfreq+1)
		# 		except:
		# 			print('Not found, adding new array')
		# 			frequencies.append(matchfreq)
		# 			skymaps.append([np.zeros(npix), np.zeros(npix)])#, np.zeros(npix), np.zeros(npix)])
		# 			hitmaps.append(np.zeros(npix))
		# 			print(np.shape(skymaps))
		# 			index = len(frequencies)-1
		try:
			index = frequencies.index(pixname)
		except:
			print('Not found, adding new array')
			frequencies.append(pixname)
			skymaps.append([np.zeros(npix), np.zeros(npix)])#, np.zeros(npix), np.zeros(npix)])
			hitmaps.append(np.zeros(npix))
			print(np.shape(skymaps))
			index = len(frequencies)-1
		print(hdul[1].data.field(0))
		print(min(hdul[1].data.field(-2)))
		print(max(hdul[1].data.field(-2)))
		print(min(hdul[1].data.field(-1)))
		print(max(hdul[1].data.field(-1)))

		# Convert to Healpix pixels
		healpix_pixel = hp.ang2pix(nside, (np.pi/2)-hdul[1].data.field(1)*np.pi/180.0, hdul[1].data.field(0)*np.pi/180.0)
		print(len(healpix_pixel))

		if len(np.unique(healpix_pixel)) > 10:

			if use_flag:
				healpix_pixel[hdul[1].data.field('flag') != 1] = -1
				print('Dome flagged:')
				print(np.sum(healpix_pixel==-1))
			if moon_mindist > 0:
				healpix_pixel[hdul[1].data.field('moondist') < moon_mindist] = -1
				print('+ Moon flagged:')
				print(np.sum(healpix_pixel==-1))
			if sun_mindist > 0:
				healpix_pixel[hdul[1].data.field('sundist') < sun_mindist] = -1
				print('+ Sun flagged:')
				print(np.sum(healpix_pixel==-1))

			num_to_average = 200
			col1 = subtractbaseline(hdul[1].data.field(3),option=0,navg=num_to_average)
			col2 = subtractbaseline(hdul[1].data.field(4),option=0,navg=num_to_average)
			# col3 = subtractbaseline(hdul[1].data.field(5),option=0,navg=num_to_average)
			# col4 = subtractbaseline(hdul[1].data.field(6),option=0,navg=num_to_average)

			print(index)
			for pixnum in np.unique(healpix_pixel):
				if pixnum != -1:
					skymaps[index][0][pixnum] = skymaps[index][0][pixnum] + np.sum(col1[healpix_pixel==pixnum])
					skymaps[index][1][pixnum] = skymaps[index][1][pixnum] + np.sum(col2[healpix_pixel==pixnum])
					# skymaps[index][2][pixnum] = skymaps[index][2][pixnum] + np.sum(col3[healpix_pixel==pixnum])
					# skymaps[index][3][pixnum] = skymaps[index][3][pixnum] + np.sum(col4[healpix_pixel==pixnum])
					hitmaps[index][pixnum] = hitmaps[index][pixnum] + np.sum(healpix_pixel==pixnum)

		hdul.close()
	filesdone += 1
	# if filesdone > 10:
	# 	break

for index in range(0,len(frequencies)):
	if np.sum(hitmaps[index]) > 10.0:
		for i in range(0,npix):
			if hitmaps[index][i] >= 1:
				skymaps[index][0][i] = skymaps[index][0][i]/hitmaps[index][i]
				skymaps[index][1][i] = skymaps[index][1][i]/hitmaps[index][i]
				# skymaps[index][2][i] = skymaps[index][2][i]#/hitmaps[index][i]
				# skymaps[index][3][i] = skymaps[index][3][i]#/hitmaps[index][i]
			else:
				skymaps[index][0][i] = hp.pixelfunc.UNSEEN
				skymaps[index][1][i] = hp.pixelfunc.UNSEEN
				# skymaps[index][2][i] = hp.pixelfunc.UNSEEN
				# skymaps[index][3][i] = hp.pixelfunc.UNSEEN


		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_0.fits',skymaps[index][0],overwrite=True)
		hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_1.fits',skymaps[index][1],overwrite=True)
		# hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_2.fits',skymaps[index][2],overwrite=True)
		# hp.write_map(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_3.fits',skymaps[index][3],overwrite=True)
		hp.write_map(outputdir+subdir+'hits_'+str(nside)+'_'+str(frequencies[index])+'.fits',hitmaps[index],overwrite=True)
		hp.mollview(skymaps[index][0])
		plt.savefig(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_0.png')
		hp.mollview(skymaps[index][1])
		plt.savefig(outputdir+subdir+'map_'+str(nside)+'_'+str(frequencies[index])+'_1.png')
		hp.mollview(hitmaps[index])
		plt.savefig(outputdir+subdir+'hits_'+str(nside)+'_'+str(frequencies[index])+'.png')
