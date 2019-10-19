from gbreduce_read import *
#from RotLog_file import *
# from RotLog_file_v2 import *
import numpy as np

indir = '/Users/mpeel/Desktop/aztest/'
compress = True

# rotlog = RotLog_file(indir+'rot_2019-1011-000002.dat.xz', compress=compress)
# print(rotlog)
# print(rotlog.start_time)
# print(rotlog.data)
# exit()

indir = '/Users/mpeel/Desktop/aztest/18/'
# rotlog = RotLog_file(indir+'az_2019-1018-033805+0000.dat.xz')
test = compress_and_plot_azdata(indir,'aztest')
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
exit()

# reload = False
# indir = '/Users/mpeel/Desktop/eltest/'
# compress_and_plot_eldata(indir,'el_2019-10-15_rot',reload=True)

indir = '/net/nas/proyectos/cosmology/groundbird/data/logbdata/el_enc/2019/'
# compress_and_plot_eldata(indir+'10/10/','el_2019-10-10')
# compress_and_plot_eldata(indir+'10/11/','el_2019-10-11')
# compress_and_plot_eldata(indir+'10/12/','el_2019-10-12')
# compress_and_plot_eldata(indir+'10/13/','el_2019-10-13')
# compress_and_plot_eldata(indir+'10/14/','el_2019-10-14')
# compress_and_plot_eldata(indir+'10/15/','el_2019-10-15')
# compress_and_plot_eldata(indir+'10/16/','el_2019-10-16')
# compress_and_plot_eldata(indir+'10/17/','el_2019-10-17')
compress_and_plot_eldata(indir+'10/18/','el_2019-10-18')
