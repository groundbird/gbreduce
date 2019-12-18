#!/usr/bin/env python
import os
import numpy as np
import datetime
import lzma

from datetime import timezone
from enum import Enum
from argparse import ArgumentParser
from pathlib import Path

PACKET_LENGTH = 12

class DataType(Enum):
    DATA = 1
    SYNC = 2
    UART = 3

def header_info(d_bytes):
    numbytes = d_bytes[0:4].decode('utf-8')
    version = int.from_bytes(d_bytes[4:8], 'little', signed=False)
    time = float(int.from_bytes(d_bytes[8:12], 'little', signed=False)) + \
        float(int.from_bytes(d_bytes[12:16], 'little',signed=False))*1e-6
    headertxt = d_bytes[16:].decode('utf-8')

    return numbytes, version, time, headertxt

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
        self._path = Path(path)
        if self._path.suffix == '.xz':
            self._fd = lzma.open(self._path, 'rb')
            self._isxz = True
        else:
            self._fd = open(self._path, 'rb')
            self._isxz = False
        self._i = 0

        # First 4 bytes: header length in bytes
        self._hlen = int(self._fd.read(4).decode('utf-8'))
        self._fd.seek(0)
        
        # Read and parse header
        header = self._fd.read(self._hlen)
        _, version, utime, txt = header_info(header)
        self.file_version = version
        self.c_utime = utime # creation unix time
        self.c_dt = datetime.datetime.fromtimestamp(utime, tz=timezone.utc) # creation datetime
        self.header_info = txt

        if self._isxz:
            self._length = None
        else:
            self._length = int((int(os.path.getsize(self._path)) - self._hlen)/PACKET_LENGTH)
        self.sync = []

    @property
    def length(self):
        if self._length is None:
            cur_pos = self._fd.tell()
            pos = self._fd.seek(0, whence=2)
            self._length = int((pos-self._hlen)/PACKET_LENGTH)
            self._fd.seek(cur_pos)
        return self._length

    def __del__(self):
        if not self._fd.closed:
            self._fd.close()

    def __iter__(self):
        self._fd.seek(self._hlen)
        return self

    def __next__(self):
        if self._i == self.length:
            raise StopIteration()
        self._i += 1            
        stamp, data, d_type = parsebytes(self._fd.read(12))
        if d_type == DataType.DATA:
            return stamp, data
        elif d_type == DataType.SYNC:
            self.sync.append((stamp, data))
            return self.__next__()
        else:
            return self.__next__()
    
    def parse_all(self):
        ret_list = []
        for stmp, data in self:
            ret_list.append([stmp, data])
        return np.array(ret_list)

    def get_data(self, cur):
        self._fd.seek(self._hlen + PACKET_LENGTH*cur)
        return parsebytes(self._fd.read(PACKET_LENGTH))

    def get_header(self):
        self._fd.seek(0)
        return header_info(self._fd.read(self._hlen))

    def get_length(self):
        cur_pos = self._fd.tell()
        pos = self._fd.seek(0, whence=2)
        length = int((pos-self._hlen)/PACKET_LENGTH)
        self._fd.seek(cur_pos)
        return length

def main():
    parser = ArgumentParser()
    parser.add_argument('path', action='store', type=str)
    args = parser.parse_args()

    eldata = ElData(args.path)
    for stamp, data in eldata:
        print(f'{stamp} {data}')


if __name__ == '__main__':
    main()
