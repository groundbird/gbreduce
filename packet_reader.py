# From Shunsuke rhea_comm
# Modified to return things in a simple array rather than nested arrays
#!/usr/bin/env python3

from struct import pack, unpack
from numpy import median
from sys import stderr

buffsize = 4096
# buffsize = 32768
header = 0xff
header_sgsync = 0xaa
header_sync = 0xf5
footer = 0xee

class packet_reader_error(Exception): pass


def get_packet_size(filename):
	ret = []
	offset = 0

	f = open(filename, 'rb')
	buff = f.read(buffsize)
	f.close()


	if buff[0] == header:
		ret += [0]
	else:
		raise packet_reader_error('error : get_packet_size.header')

	while True:
		p = buff.find(footer)
		if p == len(buff) - 1: break
		if p == -1: break
		buff = buff[p+1:]
		offset += p+1
		if buff[0] in [header, header_sync]: ret += [offset]
		pass

	ret_diff = []
	p = ret[0]
	for n in ret[1:]:
		ret_diff += [n-p]
		p = n
		pass

	return int(median(ret_diff))


def get_length(filename, packet_size = None):
	if packet_size is None: packet_size = get_packet_size(filename)
	f = open(filename, 'rb')

	buff = b''
	cnt = 0
	try:
		while True:
			buff += f.read(buffsize)
			if len(buff) < packet_size: break

			while True:
				if buff[0] not in [header, header_sync]:
					raise packet_reader_error('error : get_length.header', cnt)
				if buff[packet_size - 1] != footer:
					raise packet_reader_error('error : get_length.footer', cnt)
				cnt += 1
				buff = buff[packet_size:]
				if len(buff) < packet_size: break
				pass

			pass
	except packet_reader_error as e:
		print(e)
		print(f'packet error: {cnt}', file=stderr)
		pass

	f.close()
	return cnt


def read_packet_in_swp(buff):
	dlen = (len(buff) - 7) / 7
	t = -1
	if len(buff) != dlen * 7 + 7:
		raise packet_reader_error('error : read_packet_in_swp.size')
	if buff[0] not in [header, header_sgsync, header_sync]:
		raise packet_reader_error('error : read_packet_in_swp.header')
	if buff[-1] != footer:
		raise packet_reader_error('error : read_packet_in_swp.footer')
	if buff[0] == header:
		t = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
	return t

def read_iq_packet(buff, n_rot = -1, sync_off = 0):
	dlen = (len(buff) - 7) / 7
	if len(buff) != dlen * 7 + 7:
		raise packet_reader_error('error : read_iq_packet.size')
	if buff[0] not in [header, header_sgsync]:
		raise packet_reader_error('error : read_iq_packet.header')
	if buff[-1] != footer:
		raise packet_reader_error('error : read_iq_packet.footer')
	t = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
	data = []
	for dn in range(int(dlen)):
		d1 = unpack('b',  buff[6 + 7 * dn : 6 + 7 * dn + 1])[0]
		d2 = unpack('>H', buff[6 + 7 * dn + 1 : 6 + 7 * dn + 3])[0]
		d3 = unpack('>I', buff[6 + 7 * dn + 3 : 6 + 7 * dn + 7])[0]
		data += [(d1 << (8 * 6)) + (d2 << (8 * 4)) + d3]
		pass
	# Next line modified to return in a different format.
	# return t, data, n_rot, sync_off
	retvals = [t] + data + [n_rot] + [sync_off]
	del data
	return retvals


def read_sync_packet(buff):
	dlen = (len(buff) - 7) / 7
	if len(buff) != dlen * 7 + 7:
		raise packet_reader_error('error : read_sync_packet.size')
	if buff[0] != header_sync:
		raise packet_reader_error('error : read_sync_packet.header')
	if buff[-1] != footer:
		raise packet_reader_error('error : read_sync_packet.footer')
	n_rot = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0] # ts = n_rot
	d1 = unpack('b',  buff[6  :6+1])[0]
	d2 = unpack('>H', buff[6+1:6+3])[0]
	d3 = unpack('>I', buff[6+3:6+7])[0]
	sync_off = (d1 << (8 * 6)) + (d2 << (8 * 4)) + d3 # I_data = sync_offset
	return n_rot, sync_off

def read_snap_packet(buff, n_rot = -1, sync_off = 0):
	if len(buff) != 15:
		raise packet_reader_error('error : read_snap_packet.size')
	if buff[0] != header:
		raise packet_reader_error('error : read_snap_packet.header')
	if buff[-1] != footer:
		raise packet_reader_error('error : read_snap_packet.footer')
	t = (unpack('b', buff[1:2])[0] << (8 * 4)) + unpack('>I', buff[2:6])[0]
	d1 = unpack('>i', buff[ 6:10])[0]
	d2 = unpack('>i', buff[10:14])[0]
	return t, [d1, d2]

def seek_sync(fd, read_packet, packet_size, offset=0):
	buff = b''

	# Find real offset
	it_begin = 1
	it_end = offset + 1
	ts_begin = None
	while True:
		fd.seek(packet_size * it_begin)
		buff = fd.read(packet_size)
		if buff[0] == header:
			# ts_begin, _, _, _ = read_packet(buff)
			ts_begin = read_packet(buff)[0]
			break
		else:
			it_begin += 1
			it_end   += 1

	ts_end = ts_begin
	while ts_end - ts_begin != offset:
		fd.seek(packet_size * it_end)
		buff = fd.read(packet_size)
		if buff[0] == header:
			# ts_end, _, _, _ = read_packet(buff)
			ts_end = read_packet(buff)[0]
			it_end += offset - (ts_end - ts_begin)
		else:
			it_end += 1
	cnt = it_end

	# Find n_rot, sync_off
	n_rot = -1
	sync_off = 0

	if it_end - it_begin == offset: # no sync lines before `offset`
		return n_rot, sync_off, cnt

	# binary search
	while it_begin < it_end:
		n_rot = None
		it = (it_begin + it_end) // 2
		## check sync lines
		while it < it_end:
			fd.seek(packet_size * it)
			buff = fd.read(packet_size)
			if buff[0] == header:
				# ts, _, _, _ = read_packet(buff)
				ts = read_packet(buff)[0]
				break
			elif buff[0] == header_sync:
				n_rot, sync_off = read_sync_packet(buff)
			it += 1

		if it_end - it == ts_end - ts:
			if n_rot is not None: break
			it_end = it
			ts_end = ts
		else:
			it_begin = it
			ts_begin = ts

	return n_rot, sync_off, cnt

def _read_file(filename, packet_size = None, length = None, offset = 0, read_packet = read_iq_packet, step=1):
	if packet_size is None: packet_size = get_packet_size(filename)
	buff = b''
	readcnt = -1
	cnt = 0
	f = open(filename, 'rb')
	## HEADER
	buff += f.read(buffsize*2)
	yield read_packet(buff[0:packet_size])

	## BODY
	buff = b''
	f.seek(packet_size * (offset + 1))

	while True:
		if type(length) == int and cnt >= length: break
		if len(buff) < packet_size: buff += f.read(buffsize)
		if len(buff) < packet_size: break
		readcnt += 1
		if readcnt % step ==0:
			cnt += 1
			yield read_packet(buff[0:packet_size])[0:2]
		buff = buff[packet_size:]
		pass

	pass

def _read_file_sync(filename, packet_size = None, length = None, offset = 0, read_packet = read_iq_packet, step=1):
	if packet_size is None: packet_size = get_packet_size(filename)

	buff = b''
	cnt = 0
	readcnt = -1
	f = open(filename, 'rb')
	## HEADER
	buff += f.read(buffsize*2)
	yield read_packet(buff[0:packet_size])

	## BODY
	buff = b''
	n_rot,sync_off,offset = seek_sync(f,read_packet,packet_size,offset)
	f.seek(packet_size * offset)
	try:
		while True:
			if type(length) == int and cnt >= length: break
			if len(buff) < packet_size: buff += f.read(buffsize)
			if len(buff) < packet_size: break
			if buff[0] == header_sync:
				n_rot, sync_off = read_sync_packet(buff[0:packet_size])
			else:
				readcnt += 1
				if readcnt % step ==0:
					cnt += 1
					yield read_packet(buff[0:packet_size], n_rot, sync_off)
					pass
			buff = buff[packet_size:]
			pass
		pass

	except packet_reader_error as e:
		print(e)

# def _read_file_sync(filename, packet_size = None, length = None, offset = 0, read_packet = read_iq_packet, step=1):
# 	if packet_size is None: packet_size = get_packet_size(filename)
#
# 	f = open(filename, 'rb')
# 	n_rot,sync_off,offset = seek_sync(f,read_packet,packet_size,offset)
# 	f.close()
# 	# buff = b''
# 	cnt = 0
# 	readcnt = -1
# 	with open(filename, "rb") as f:
# 		bytes_read = f.read()
# 	## HEADER
# 	# buff += f.read(buffsize*2)
# 	yield read_packet(bytes_read[0:packet_size])
# 	up_to = packet_size
# 	total_length = len(bytes_read)
# 	## BODY
# 	# buff = b''
# 	up_to += packet_size*offset
# 	# f.seek(packet_size * offset)
# 	try:
# 		while True:
# 			if type(length) == int and cnt >= length: break
# 			# if len(buff) < packet_size: buff += f.read(buffsize)
# 			if total_length - up_to < packet_size: break
# 			if bytes_read[up_to] == header_sync:
# 				n_rot, sync_off = read_sync_packet(bytes_read[up_to:up_to+packet_size])
# 			else:
# 				readcnt += 1
# 				if readcnt % step ==0:
# 					cnt += 1
# 					yield read_packet(bytes_read[up_to:up_to+packet_size], n_rot, sync_off)
# 					pass
# 			up_to += packet_size
# 			# buff = buff[packet_size:]
# 			pass
# 		pass
#
# 	except packet_reader_error as e:
# 		print(e)

def read_iq_file(filename, packet_size = None, length = None, offset = 0):
	return _read_file(filename, packet_size = packet_size,
					  length = length, offset = offset,
					  read_packet = read_iq_packet)


def read_snap_file(filename, length = None, offset = 0):
	return _read_file(filename, packet_size = 15,
					  length = length, offset = offset,
					  read_packet = read_snap_packet)


def read_file(filename, packet_size = None, length = None, offset = 0, sync=False, step=1):
	if packet_size is None: packet_size = get_packet_size(filename)
	read_packet = read_iq_packet
	if packet_size == 15: read_packet = read_snap_packet
	if not sync:
		return _read_file(filename, packet_size = packet_size,
						  length = length, offset = offset, step = step,
						  read_packet = read_packet)
	else:
		return _read_file_sync(filename, packet_size = packet_size,
							   length = length, offset = offset, step = step,
							   read_packet = read_packet)

# def read_iq_file_back(filename, packet_size = None, length = None, offset = 0):
#     if packet_size is None: packet_size = get_packet_size(filename)
#     if length is None:
#         length = get_length(filename, packet_size = packet_size)
#         length -= offset
#         pass

#     buff = b''
#     f = open(filename, 'rb')
#     f.seek(packet_size * offset)

#     for i in range(length):
#         if len(buff) < packet_size: buff += f.read(buffsize)
#         if len(buff) < packet_size: raise packet_reader_error('end_of_file')
#         yield read_iq_packet(buff[0:packet_size])
#         buff = buff[packet_size:]
#         pass

#     pass

# def read_snap_file_back(filename, length = None, offset = 0):
#     packet_size = 15
#     if length is None:
#         length = get_length(filename, packet_size = packet_size)
#         length -= offset
#         pass

#     buff = b''
#     f = open(filename, 'rb')
#     f.seek(packet_size * offset)

#     for i in range(length):
#         if len(buff) < packet_size: buff += f.read(buffsize)
#         if len(buff) < packet_size: raise packet_reader_error('end_of_file')
#         yield read_snap_packet(buff[0:packet_size])
#         buff = buff[packet_size:]
#         pass

#     pass


def test1(filename):
	print(get_packet_size(filename))
	return

def test2(filename):
	print(get_length(filename))
	return

def test3(fname, length = None, offset = 0):
	for t, d in read_file(fname,
						  length = length,
						  offset = offset):
		print(t, end=' ')
		for dd in d: print(dd, end=' ')
		print()
		pass
	return

if __name__ == '__main__':
	from sys import argv
	#test1(argv[1])
	#test2(argv[1])
	#test3(argv[1])
	#test3(argv[1], length = 20000)
	#test3(argv[1], length = 2000)
	#test3(argv[1], length = 20, offset = 1000)
	exit(0)
