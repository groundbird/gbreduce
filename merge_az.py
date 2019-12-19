#!/usr/bin/env python3
# Code by Shugo, tweaks by Mike

import numpy as np
from datetime import datetime as dt

def merge_az(kid_data, kid_sync_n, kid_sync_off, utime_az, azdata_az, sync_n_az, sync_off_az):
    sync_n   = np.array(kid_sync_n,dtype=np.int)
    # print(sync_n)
    # print(len(sync_n))
    sync_off = np.array(kid_sync_off) * 5e-9 # sec
    kiddata  = kid_data

    start = sync_n.searchsorted(0)
    # print(start)
    sync_n   = sync_n[start:]
    sync_off = sync_off[start:]
    kiddata  = kiddata[start:]
    utime    = np.zeros(len(kiddata))
    azdata   = np.zeros(len(kiddata))# - 180.

    utime_az = np.asarray(utime_az)
    # print(utime_az)
    azdata_az = np.asarray(azdata_az)
    # print(azdata_az)
    sync_n_az   = np.asarray(sync_n_az.copy()).astype(int)

    # sync_n_az   = np.asarray(sync_n_az,dtype=np.int)
    # print(sync_n_az)
    sync_off_az = np.asarray(sync_off_az) * 40e-9 # sec
    # print(sync_off_az)

    for sn in range(sync_n[0], sync_n[-1]+1):
        assert sn in sync_n_az
        # print(sn)
        # print(sn in sync_n_az)
        start   = sync_n.searchsorted(sn)
        end     = sync_n.searchsorted(sn+1)
        # print('\n')
        # print(sn)
        # print(start)
        # print(end)
        offset_kid = sync_off[start]
        # print(offset_kid)
        vals = np.where(sync_n_az==sn)
        start_az = np.min(vals)
        # end_az = np.max(vals)+1
        # start_az = sync_n_az.searchsorted(int(sn))
        # print(start_az)
        end_az   = start_az + (end - start)
        # print(end_az)
        offset_az  = sync_off_az[start_az]
        # print(offset_az)
        offset = offset_az - offset_kid
        # print(offset)
        offset_sign = 1 if offset > 0 else -1
        offset_abs  = abs(offset)
        # print(azdata_az[start_az])
        # print(azdata_az[end_az])
        tmp_az_0 = azdata_az[start_az             : end_az            ]
        tmp_az_1 = azdata_az[start_az+offset_sign : end_az+offset_sign]
        tmp_az_1 = tmp_az_1 + np.round((tmp_az_0 - tmp_az_1) / 360.) * 360.

        tmp_az = tmp_az_0 * (1 - offset_abs/0.001) + tmp_az_1 * (offset_abs/0.001)
        azdata[start:end] = tmp_az % 360.

        # if sn == 5123626:
        #     print(azdata[start-1])
        #     print(azdata[start:end])

        #     val = input('wait')

        # if np.abs(azdata[start-1]-azdata[start]) > 10.0:
        #     print(azdata[start-1])
        #     print(azdata[start:end])
        #     val = input('huh?')
        # print(len(tmp_az))
        tmp_ut_0 = utime_az[start_az             : end_az            ]
        tmp_ut_1 = utime_az[start_az+offset_sign : end_az+offset_sign]
        tmp_ut = tmp_ut_0 * (1 - offset_abs/0.001) + tmp_ut_1 * (offset_abs/0.001)
        # print(len(tmp_ut))
        kiddata[start:end,0] = tmp_ut
        pass
    # time = np.array(list(map(dt.fromtimestamp, utime)))

    return kiddata, azdata

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
    dir0 = "/data/gb/kiddata/20191112/data_145400/"
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
        klib.kidana.KidAnalyzer(swp=swpset[6], tod=todset[6], ctod=todset[7]),
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
    from az_read.RotLog import RotLog

    rot = RotLog('2019-1112-145000',
                 '2019-1112-151600',
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


if __name__ == '__main__':
    test1()
