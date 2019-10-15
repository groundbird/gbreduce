from compress_eldata import *
from RotLog_file import *

indir = '/Users/mpeel/Desktop/aztest/'
compress = True

rotlog = RotLog_file(indir+'rot_2019-1011-000002.dat.xz', compress=compress)
print(rotlog)
print(rotlog.start_time)
print(rotlog.data)
exit()

# reload = False
# indir = '/Users/mpeel/Desktop/eltest/'
# compress_and_plot_eldata(indir,'el_2019-10-15_rot',reload=True)

indir = '/net/nas/proyectos/cosmology/groundbird/logbdata/el_enc/2019/'
compress_and_plot_eldata(indir+'10/10/','el_2019-10-10')
compress_and_plot_eldata(indir+'10/11/','el_2019-10-11')
compress_and_plot_eldata(indir+'10/12/','el_2019-10-12')
compress_and_plot_eldata(indir+'10/13/','el_2019-10-13')
compress_and_plot_eldata(indir+'10/14/','el_2019-10-14')
compress_and_plot_eldata(indir+'10/15/','el_2019-10-15')
