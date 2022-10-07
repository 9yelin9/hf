# pyhf3/wan2lat.py : extract lattice information from wannier90_hr.dat

import re

def TransWan2Lat(input_path):
	f_wan = open('%s/wannier90_hr.dat' % (input_path), 'r')
	f_lat = open('%s/lattice.txt' % (input_path), 'w')

	pat_site = '[-]?\d+\s+'
	pat_obt = '[-]?\d+\s+'
	pat_t = '[-]?\d+[.]\d+\s+'

	pat = 3*pat_site + 2*pat_obt + 2*pat_t
	len = 0

	for line in f_wan:
		if re.search(pat, line): len += 1
	f_lat.write("%d\n" % (len))

	f_wan.seek(0, 0)
	for line in f_wan:
		if re.search(pat, line): f_lat.write(line) 

	f_wan.close()
	f_lat.close()
