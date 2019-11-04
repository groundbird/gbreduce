""" read a rotary encoder log file"""
#!/usr/bin/env python3
# From https://gitlab.emperor.riken.jp/taku/fpga_gb_rotary_ctrl_12/blob/master/comm/RotLog_file.py

from enum import Enum
from struct import unpack
from datetime import datetime
import pandas as pd
import numpy as np
import lzma

class PacketError(Exception):
    """Exception"""
    def __init__(self, key):
        super().__init__(key)
        self.key = key
    def print(self):
        print(f'packet error. ({self.key})')

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

def conv_angle(_v):
    return float(2. * np.pi * _v / 8192.)

class RotLog_file_old():
    """read file"""
    def __init__(self, _path, n_rotate=-1, sync_off=0, compress=True):
        self._packetsize = 8
        self._path = _path
        if compress is True:
            with lzma.open(self._path, mode='rb') as f:
                self.rawdata = f.read()
                pass
        else:
            with open(self._path, 'rb') as f:
                self.rawdata = f.read()
                pass
        self._length = int(len(self.rawdata)/self._packetsize)
        self._i = 0
        self._nrot = n_rotate
        self._sync_off = sync_off
        self.start_time = self._read_header()
        self.data = self._read_file()

    def __iter__(self):
        return self

    def __next__(self):
        if self._i == self._length:
            raise StopIteration()
        buff = self.rawdata[self._i*self._packetsize:(self._i+1)*self._packetsize]
        timestamp, data, ptype = read_packet(buff)
        self._i += 1
        if ptype == PacketType.DATA:
            return timestamp, data
        elif ptype == PacketType.SYNC:
            self._nrot = timestamp
            self._sync_off = data
            return self.__next__()

    def _read_header(self):
        buff = self.rawdata[self._i*self._packetsize:(self._i+1)*self._packetsize]
        _t, _d, _p = read_packet(buff)
        self._i += 1
        return np.float64(_t + _d / 10000.)

    def _read_file(self):
        ret_t = []
        ret_nr = []
        ret_off = []
        ret_d = []
        for _ts, _d in self:
            ret_t.append(_ts)
            ret_nr.append(self._nrot)
            ret_off.append(self._sync_off)
            ret_d.append(conv_angle(_d))
        data = pd.DataFrame({DataType.timestamp:ret_t, DataType.n_rotate:ret_nr,
                             DataType.sync_off:ret_off, DataType.az_angle:ret_d})
        return data

def main(path_str, compress=True):
    rotlog = RotLog_file(path_str, compress=compress)
    return rotlog.start_time, rotlog.data

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('path', action='store', type=str)
    args = parser.parse_args()
    HEADER, D_OUT = main(args.path)
    print(f'unixtime in header: {HEADER}')
    print(f'datetime in header: {datetime.fromtimestamp(HEADER)}')
    print(D_OUT)
