#!/usr/bin/env python
import os
from enum import Enum
import numpy as np
from argparse import ArgumentParser

class DataType(Enum):
    DATA = 1
    SYNC = 2
    UART = 3

def header_info(d_bytes):
    numbytes = d_bytes[0:4].decode('utf-8')
    version = int.from_bytes(d_bytes[4:8], 'little',signed=False)
    time = float(int.from_bytes(d_bytes[8:12], 'little',signed=False)) + float(int.from_bytes(d_bytes[12:16], 'little',signed=False))*1e-6
    headertxt = d_bytes[16:].decode('utf-8')

    return [numbytes, version, time, headertxt]

def parsebytes(d_bytes):
    if d_bytes[0:2] != b'\x07\x12':
        raise Exception('HEADER ERROR: {}'.format(d_bytes[0:2]))
    if d_bytes[10:12] == b'z\xda': # DATA
        dt = DataType.DATA
    elif d_bytes[10:12] == b'\x0cW': # SYNC
        dt = DataType.SYNC
    elif d_bytes[10:12] == b'H ': # UART
        dt = DataType.UART
    else:
        raise Exception('FOOTER ERROR ', d_bytes[10:12])
    return int.from_bytes(d_bytes[2:6], 'little'), int.from_bytes(d_bytes[6:10], 'little', signed=True), dt

class ElData:
    def __init__(self, path):
        self._path = path
        try:
            self._fd = lzma.open(path, mode='rb')
            self._fd.read(1024)
            self._fd.seek(0)
        except:
            self._fd = open(path, 'rb')
            pass
        self._i = 0
        self._length = int((int(os.path.getsize(self._path)) - 256)/12)
        self.sync = []

    def __del__(self):
        if not self._fd.closed:
            self._fd.close()

    def __iter__(self):
        return self

    def __next__(self):
        if self._i == self._length:
            raise StopIteration()
        self._i += 1            
        stamp, data, d_type = parsebytes(self._fd.read(12))
        if d_type == DataType.DATA:
            return stamp, data
        elif d_type == DataType.SYNC:
            self.sync.append((stamp, data))
            return self.__next__()
    
    def parse_all(self):
        ret_list = []
        for stmp, data in self:
            ret_list.append([stmp, data])
        return np.array(ret_list)

    def get_data(self, cur):
        self._fd.seek(256+12*cur)
        return parsebytes(self._fd.read(12))

    def get_header(self):
        self._fd.seek(0)
        return header_info(self._fd.read(256))

def main():
    parser = ArgumentParser()
    parser.add_argument('path', action='store', type=str)
    args = parser.parse_args()

    eldata = ElData(args.path)
    for stamp, data in eldata:
        print(f'{stamp} {data}')


if __name__ == '__main__':
    main()
