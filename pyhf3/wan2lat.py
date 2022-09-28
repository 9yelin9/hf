# pyhf3/wan2lat.py : extract lattice information from wannier90_hr.dat

import re

class Wan2Lat:
	def Run(self, material):
		f_wan = open('%s/input/wannier90_hr.dat' % (material), 'r')
		f_lat = open('%s/input/lattice.txt' % (material), 'w')

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
