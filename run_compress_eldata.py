from gbreduce_read import *
#from RotLog_file import *
# from RotLog_file_v2 import *
import numpy as np
import os

# indir = '/Users/mpeel/Desktop/aztest/'
# compress = True

# rotlog = RotLog_file(indir+'rot_2019-1011-000002.dat.xz', compress=compress)
# print(rotlog)
# print(rotlog.start_time)
# print(rotlog.data)
# exit()

# indir = '/Users/mpeel/Desktop/aztest/18/'
# rotlog = RotLog_file(indir+'az_2019-1018-033805+0000.dat.xz')
# test = compress_and_plot_azdata(indir,'aztest')
# az_t, az_val = get_az_data(indir, 'az_2019-1018-033805+0000.dat.xz')
# print(az_t)
# print(az_val)
# print(rotlog.start_time)
# ret_t, ret_nr, ret_off, ret_d = rotlog.read_file()
# print(min(ret_d))
# print(max(ret_d))
# print(np.mean(ret_d))
# exit()
# print(rotlog.start_time)
# print(rotlog.length)
# for _ts, _d, _nrot, _sync_off in rotlog:
#     print(_ts, _d, _nrot, _sync_off)

# print(rotlog)
# print(rotlog._d)
# print(rotlog._nrot)
# print(rotlog._sync_off)

# print(rotlog)
# print(rotlog.start_time)
# print(data)
# print(rotlog.data)
# exit()

# reload = False
# indir = '/Users/mpeel/Desktop/eltest/'
# compress_and_plot_eldata(indir,'el_2019-10-15_rot',reload=True)

basedir = '/net/nas/proyectos/cosmology/groundbird/'
# basedir = '/Volumes/iOmega/GroundBIRD/'

indir = basedir+'data/logbdata/el_enc/2019/'
folders = [f.path for f in os.scandir(indir) if f.is_dir() ]
redo_all = True
print(folders)
for folder in folders:
	# print(folder)
	folders2 = [f.path for f in os.scandir(folder) if f.is_dir() ]
	for folder2 in folders2:
		# print(folder2)
		# print((folder2+'/')[-3:-1])
		# print(folder2[-5:-3])
		out_prefix = 'el_'+str(2019)+'-'+str(folder2[-5:-3])+'-'+str((folder2+'/')[-3:-1])
		if redo_all == True or (redo_all == False and not os.path.exists(folder2+'/'+out_prefix+'.txt.gz')):
			print(folder2)
			compress_and_plot_eldata(folder2+'/',out_prefix)



indir = basedir+'data/logbdata/az_enc/2019/'
folders = [f.path for f in os.scandir(indir) if f.is_dir() ]
redo_all = True
print(folders)
for folder in folders:
	# print(folder)
	folders2 = [f.path for f in os.scandir(folder) if f.is_dir() ]
	for folder2 in folders2:
		# print(folder2)
		# print((folder2+'/')[-3:-1])
		# print(folder2[-5:-3])
		out_prefix = 'az_'+str(2019)+'-'+str(folder2[-5:-3])+'-'+str((folder2+'/')[-3:-1])
		if redo_all == True or (redo_all == False and not os.path.exists(folder2+'/'+out_prefix+'.txt.gz')):
			print(folder2)
			compress_and_plot_azdata(folder2+'/',out_prefix,reload=True)
			# exit()
