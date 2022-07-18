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

	def SeperateLat(self):
		lat = open('input/lattice.txt', 'r')
		onsite = open('input/lattice_onsite.txt', 'w')
		super = open('input/lattice_super.txt', 'w')
		remain = open('input/lattice_remain.txt', 'w')

		pat_onsite = '\s+0\s+0\s+0\s+'
		pat_super = '\s+1\s+0\s+0\s+'

		for line in lat:
			if re.match(pat_onsite, line): onsite.write(line)
			elif re.match(pat_super, line): super.write(line)
			else: remain.write(line)

		lat.close()
		onsite.close()
		super.close()
		remain.close()
		
