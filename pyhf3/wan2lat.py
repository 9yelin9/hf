# Wan2Lat.py : extract lattice information from wannier90_hr.dat

import re

class Wan2Lat:
	def Extract(self):
		wan = open('input/wannier90_hr.dat', 'r')
		lat = open('input/lattice.txt', 'w')

		pat_site = '[-]?\d+\s+'
		pat_obt = '[-]?\d+\s+'
		pat_t = '[-]?\d+[.]\d+\s+'
		pat = pat_site + pat_site + pat_site + pat_obt + pat_obt + pat_t + pat_t

		for line in wan:
			if re.search(pat, line): lat.write(line) 

		wan.close()
		lat.close()
