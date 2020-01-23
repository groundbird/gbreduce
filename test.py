import datetime
import sys
from datetime import timezone
from gbreduce_read import *
import numpy as np
import matplotlib.pyplot as plt

files = ['az_2019-1216-000540+0000.dat.xz','az_2019-1216-002629+0000.dat.xz','az_2019-1216-004717+0000.dat.xz','az_2019-1216-010806+0000.dat.xz','az_2019-1216-012855+0000.dat.xz','az_2019-1216-014944+0000.dat.xz','az_2019-1216-021032+0000.dat.xz','az_2019-1216-023121+0000.dat.xz','az_2019-1216-025210+0000.dat.xz','az_2019-1216-031259+0000.dat.xz','az_2019-1216-033347+0000.dat.xz','az_2019-1216-035436+0000.dat.xz','az_2019-1216-041525+0000.dat.xz','az_2019-1216-043614+0000.dat.xz','az_2019-1216-045702+0000.dat.xz','az_2019-1216-051751+0000.dat.xz','az_2019-1216-053840+0000.dat.xz','az_2019-1216-055928+0000.dat.xz']
for file in files:
	test_t, test_az, test_nr, test_roff = get_az_data('/Volumes/iOmega/GroundBIRD/data/logbdata/az_enc/2019/12/16/'+file)
	print(file)
	print(test_az[0])
	print(test_az[-1])

	azdiff = np.diff(test_az)
	azdiff[azdiff < -300] += 360.0
	azdiff[azdiff > 300] -= 360.0
	plt.plot(azdiff)
	plt.savefig('test_'+file+'.png')
	plt.clf()
	plt.plot(test_az)
	plt.savefig('test2_'+file+'.png')
	print(np.max(azdiff))
	plt.clf()

exit()

current_time = datetime.datetime.now()
utime = current_time.timestamp()
print(utime)
utime_int = int(utime)
print(utime_int) # 4 bytes, integer part of the current time in unix time
print(int((utime - utime_int)*1e6)) # microseconds

