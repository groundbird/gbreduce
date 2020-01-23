#!/usr/bin/env python3
""" read a rotary encoder log file"""

from enum import Enum
from struct import unpack
from datetime import datetime
import numpy as np
import lzma

class PacketError(Exception):
    """Exception"""
    def __init__(self, key):
        super().__init__(key)
        self.key = key
    def print(self):
        print('packet error. (' + str(self.key) +')')

class PacketType(Enum):
    DATA = 0
    SYNC = 1

class DataType(Enum):
    timestamp = 0
    n_rotate = 1
    sync_off = 2
    az_angle = 3

def read_packet(buff):
    try:
        if buff[0] == 0xff:
            ptype = PacketType.DATA
        elif buff[0] == 0xf5:
            ptype = PacketType.SYNC
        else: raise PacketError('header')
        time = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
        data = unpack('>H', buff[6:8])[0]
    except PacketError as _e:
        _e.print()
    return time, data, ptype

# Originally -55.0, changed to -42.0 on 17-Jan-2020
def conv_angle(_v):
    return (-(float(360.0 * _v / 8192.))-42.0) % 360.0

# def conv_angle(_v):
#     return float(360.0 * _v / 8192.)

class RotLog_file():
    """read file"""
    def __init__(self, _path, n_rotate=-1, sync_off=0):
        self._headersize = 256
        self._packetsize = 8
        self._path = _path
        try:
            self._f = lzma.open(self._path, mode='rb')
            self._f.read(1024)
            self._f.seek(0)
        except:
            self._f = open(self._path, 'rb')
            pass
        self._length = None
        self._nrot_init = n_rotate
        self._sync_off_init = sync_off
        self._nrot = n_rotate
        self._sync_off = sync_off
        self._read_header()
        assert self.version >= 2019101401

    def __iter__(self):
        return self

    def __next__(self):
        # read file
        buff = self._f.read(self._packetsize)
        if len(buff) != self._packetsize: raise StopIteration()
        timestamp, data, ptype = read_packet(buff)
        if ptype == PacketType.DATA:
            return timestamp, data, self._nrot, self._sync_off
        elif ptype == PacketType.SYNC:
            self._nrot = timestamp
            self._sync_off = data
            return self.__next__()
        pass

    @property
    def length(self):
        if self._length is None:
            pos = self._f.tell()
            self._f.seek(0, 2)
            endpos = self._f.tell()
            self._f.seek(pos)
            self._length = (endpos - self._headersize) // self._packetsize
            pass
        return self._length

    def seek(self, data_number = 0, set_nrot = True):
        assert data_number < self.length
        if set_nrot: self._set_nrot(data_number)
        self._seek(data_number)
        pass

    def _set_nrot(self, data_number):
        # check packet type
        self._seek(data_number)
        buff = self._f.read(self._packetsize)
        timestamp, data, ptype = read_packet(buff)
        if ptype == PacketType.SYNC: return
        # check if sync signal exists
        self._nrot = None
        it_end   = data_number
        ts_end   = timestamp
        self._seek(0)
        timestamp, data, nrot, sync_off = self.__next__()
        it_begin = self._tell() - 1
        ts_begin = timestamp
        if it_end - it_begin == ts_end - ts_begin:
            if nrot == None:
                self._nrot     = self._nrot_init
                self._sync_off = self._sync_off_init
                pass
            return
        # binary search
        while it_begin < it_end:
            self._nrot = None
            it = (it_begin + it_end) // 2
            self._seek(it)
            timestamp, data, nrot, sync_off = self.__next__()
            it = self._tell() - 1
            ts = timestamp
            if it_end - it == ts_end - ts:
                if nrot is not None: break
                it_end = it
                ts_end = ts
            else:
                it_begin = it
                ts_begin = ts
                pass
            pass
        return

    def _read_header(self):
        self._f.seek(0)
        self._header = self._f.read(self._headersize)
        assert self._header[0:4] == b'%3d\n' % self._headersize
        self.version = unpack('<I', self._header[4:8])[0]
        self.start_time_sec = unpack('<I', self._header[ 8:12])[0]
        self.start_time_us  = unpack('<I', self._header[12:16])[0]
        self.start_time = self.start_time_sec + self.start_time_us * 1.e-6
        self.description = self._header[17:]
        return

    def _seek(self, data_number):
        seek_pos = self._headersize + self._packetsize * data_number
        self._f.seek(seek_pos)
        return

    def _tell(self):
        seek_pos = self._f.tell()
        return (seek_pos - self._headersize) // self._packetsize

    def read_file(self):
        import pandas as pd
        ret_t = []
        ret_nr = []
        ret_off = []
        ret_d = []
        self.seek(0)
        for _ts, _d, _nrot, _sync_off in self:
            ret_t.append(_ts)
            ret_nr.append(self._nrot)
            ret_off.append(self._sync_off)
            ret_d.append(conv_angle(_d))
        # data = pd.DataFrame({DataType.timestamp:ret_t, DataType.n_rotate:ret_nr,
                             # DataType.sync_off:ret_off, DataType.az_angle:ret_d})
        return ret_t, ret_nr, ret_off, ret_d
    pass

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('path', action='store', type=str)
    parser.add_argument('length', nargs='?', type=int, default=10)
    args = parser.parse_args()
    rotlog = RotLog_file(args.path)
    print('unixtime in header:', rotlog.start_time)
    print('datetime in header:', datetime.fromtimestamp(rotlog.start_time))
    print('data size:', rotlog.length)
    print('== start %d data ==' % args.length)
    cnt = 0
    for _ts, _d, _nrot, _sync_off in rotlog:
        if cnt >= args.length: break
        print(_ts, _d, _nrot, _sync_off)
        cnt += 1
        pass
    print('== last %d data ==' % args.length)
    rotlog.seek(rotlog.length - args.length, set_nrot = True)
    for _ts, _d, _nrot, _sync_off in rotlog:
        print(_ts, _d, _nrot, _sync_off)
        pass
    pass

def main2():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('path', action='store', type=str)
    args = parser.parse_args()
    rotlog = RotLog_file(args.path)
    print('unixtime in header:', rotlog.start_time)
    print('datetime in header:', datetime.fromtimestamp(rotlog.start_time))
    print('data size:', rotlog.length)
    for _ts, _d, _nrot, _sync_off in rotlog:
        print(_ts, _d, _nrot, _sync_off)
        pass

if __name__ == '__main__':
    main()
    #main2()
