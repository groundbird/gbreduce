#!/usr/bin/env python3

import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

def merge_az(kid_data, kid_sync_n, kid_sync_off, az):
    sync_n   = np.array(kid_sync_n)
    sync_off = np.array(kid_sync_off) * 5e-9 # sec
    kiddata  = kid_data

    start = sync_n.searchsorted(0)
    sync_n   = sync_n[start:]
    sync_off = sync_off[start:]
    kiddata  = kiddata[start:]
    utime    = np.zeros(len(kiddata))
    azdata   = np.zeros(len(kiddata)) - 180.

    utime_az = az[:,1]
    print(utime_az)
    azdata_az = az[:,2]/np.pi*180.
    print(azdata_az)
    sync_n_az   = az[:,3].astype(int)
    print(sync_n_az)
    sync_off_az = az[:,4] * 40e-9 # sec
    print(sync_off_az)
    # exit()

    for sn in range(sync_n[0], sync_n[-1]+1):
        assert sn in sync_n_az

        start   = sync_n.searchsorted(sn)
        end     = sync_n.searchsorted(sn+1)

        print('\n')
        print(sn)
        print(start)
        print(end)

        offset_kid = sync_off[start]
        print(offset_kid)

        start_az = sync_n_az.searchsorted(sn)
        print(start_az)
        end_az   = start_az + (end - start)
        print(end_az)
        offset_az  = sync_off_az[start_az]
        print(offset_az)

        offset = offset_az - offset_kid
        print(offset)

        offset_sign = 1 if offset > 0 else -1
        offset_abs  = abs(offset)
        print(azdata_az[start_az])
        print(azdata_az[end_az])

        tmp_az_0 = azdata_az[start_az             : end_az            ]
        tmp_az_1 = azdata_az[start_az+offset_sign : end_az+offset_sign]
        tmp_az_1 = tmp_az_1 + np.round((tmp_az_0 - tmp_az_1) / 360.) * 360.

        tmp_az = tmp_az_0 * (1 - offset_abs/0.001) + tmp_az_1 * (offset_abs/0.001)
        azdata[start:end] = tmp_az % 360.

        tmp_ut_0 = utime_az[start_az             : end_az            ]
        tmp_ut_1 = utime_az[start_az+offset_sign : end_az+offset_sign]
        tmp_ut = tmp_ut_0 * (1 - offset_abs/0.001) + tmp_ut_1 * (offset_abs/0.001)
        utime[start:end] = tmp_ut

        if sn == 5123626:
            print(azdata[start-1])
            print(azdata[start:end])

            val = input('wait')

        pass
    time = np.array(list(map(dt.fromtimestamp, utime)))

    return time, utime, kiddata, azdata

def test1():
    ##########################
    ## FROM SHUNSUKE'S CODE ##
    ##########################
    print('read kid data')
    import mkid_pylibs as klib
    import glob

    ## SG FREQ
    lo = "4.99GHz"
    #dir0 = "/home/monitor/data_local/kiddata/20191107/data_171818/"
    # dir0 = "/data/gb/kiddata/20191112/data_145400/"
    dir0 = "/Volumes/iOmega/GroundBIRD/data/kiddata/20191218/data_033616/"
    ## SWPDATASET
    sfp0 = glob.glob(dir0+"/swp*.rawdata")[0]
    ## TODDATASET
    tfp0 = glob.glob(dir0+'/tod*.rawdata')[0]
    print(sfp0)
    print(tfp0)

    swpset = klib.readfile_swp('rhea',sfp0,-1,lo)
    todset = klib.readfile_tod('rhea',tfp0,-1,lo)

    k = [
        klib.kidana.KidAnalyzer(swp=swpset[0], tod=todset[0], ctod=todset[1]),
        klib.kidana.KidAnalyzer(swp=swpset[2], tod=todset[2], ctod=todset[3]),
        klib.kidana.KidAnalyzer(swp=swpset[4], tod=todset[4], ctod=todset[5]),
        # klib.kidana.KidAnalyzer(swp=swpset[6], tod=todset[6], ctod=todset[7]),
    ]

    for ik in k:
        ik.fitIQ(nfwhm=3)

    tod = k[0].tod
    kid_data = tod.rwmdata.corphase
    kid_sync_n = tod.info['n_rot']
    kid_sync_off = tod.info['sync_off']


    #######################
    ## READ AZIMUTH DATA ##
    #######################
    print('read azimuth data')
    from groundbird.az_read.RotLog import RotLog

    rot = RotLog('2019-1218-032554',
                 '2019-1218-050958',
                 '%Y-%m%d-%H%M%S',
                 angle_output_by_radian = True)
    az = np.array(list(rot))


    ################
    ## MERGE THEM ##
    ################
    print('merge them')
    time, utime, kiddata, azdata = merge_az(kid_data, kid_sync_n, kid_sync_off, az)

    print('time', len(time))
    print(time[0:5])
    print('utime', len(utime))
    print(utime[0:5])
    print('kiddata', len(kiddata))
    print(kiddata[0:5])
    print('azdata', len(azdata))
    print(azdata[0:5])
    plt.plot(azdata)
    plt.savefig('merge_az_after.png')
    plt.clf()
    plt.plot(az[:,2]/np.pi*180.)
    plt.savefig('merge_az_before.png')
    plt.clf()

if __name__ == '__main__':
    test1()
