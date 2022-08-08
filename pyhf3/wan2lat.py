# pyhf3/wan2lat.py : extract lattice information from wannier90_hr.dat

import re

class Wan2Lat:
	def runWan2Lat(self, material):
		f_wan = open('%s/input/wannier90_hr.dat' % (material), 'r')
		f_lat = open('%s/input/lattice.txt' % (material), 'w')

		pat_site = '[-]?\d+\s+'
		pat_obt = '[-]?\d+\s+'
		pat_t = '[-]?\d+[.]\d+\s+'

		pat = 3*pat_site + 2*pat_obt + 2*pat_t

		for line in wan:
			if re.search(pat, line): f_lat.write(line) 

		f_wan.close()
		f_lat.close()
