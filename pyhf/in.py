# pyhf/in.py : class to generate input

import os
import re
import numpy as np
import pandas as pd

class In:
	def __init__(self, name):
		self.path_input = "input/%s" % name
		os.makedirs(self.path_input, exist_ok=True)

		self.lat_dict = {
			'sc' : [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
			'fcc': [[0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0.5, 0.5]],
		}
	
	def Wan2Lat(self):
		fwn = '%s/wannier90_hr.dat' % self.path_input
		fln = '%s/lat.txt'          % self.path_input

		pat_site = '[-]?\d+\s+'
		pat_obt  = '[-]?\d+\s+'
		pat_t    = '[-]?\d+[.]\d+\s+'
		pat      = 3*pat_site + 2*pat_obt + 2*pat_t

		with open(fwn, 'r') as f:
			f_len = sum(1 for line in f if re.search(pat, line))

		with open(fwn, 'r') as fw, open(fln, 'w') as fl:
			for line in fw:
				if re.search(pat, line): fl.write(line)

	def Wan2Path(self, Nk=1024):
		fwn = '%s/wannier90.win' % self.path_input
		fpn = '%s/path.txt'      % self.path_input

		a = 

