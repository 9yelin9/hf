# pyhf3/magspec.py : predict magnetic structure from spectrum

import os
import re
import sys
import h5py
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class MagSpec:
	def __init__(self, name='baoso3_dU0.2'):
		self.BINS_MAX = 1024
		self.Nkb = 1024
		self.Nbs = 12
		self.Nhsp = 4
		self.emin = -8
		self.emax = 0
		self.points = [0, 198, 396, 677, 1023]
		self.labels = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']
		self.type_dict = {'a':1, 'c':2, 'g':3}
		self.type_dict_r = {'1':'a', '2':'c', '3':'g'}
		self.params = ['type', 'JU', 'N', 'U', 'm', 'gap']
		self.path_output = 'output/%s/' % name
		self.path_save = 'output/%s/magspec/' % name
		os.makedirs(self.path_save, exist_ok=True)

	def GenEnergy(self):
		t0 = time.time()

		en = '%s/energy.h5' % self.path_save
		e = np.linspace(self.emin, self.emax, self.BINS_MAX)
		with h5py.File(en, 'w') as f: f.create_dataset('energy', data=e, dtype='d')

		t1 = time.time()
		print('GenEnergy(%s) : %fs' % (en, t1-t0))
	
	def ReadFn_(self, fn):
		fn_dict = {
			'type':  self.type_dict[re.sub('band_', '', re.search('band_[a-z]\d',      fn).group())[0]],
			'JU':    float(re.sub('JU',     '', re.search('JU\d+[.]\d+',         fn).group())),
			'N':     float(re.sub('_N',     '', re.search('_N\d+[.]\d+',         fn).group())),
			'U':     float(re.sub('_U',     '', re.search('_U\d+[.]\d+',         fn).group())),
			'n':     float(re.sub('_n',     '', re.search('_n\d+[.]\d+',         fn).group())),
			'm':     float(re.sub('_m',     '', re.search('_m[-]?\d+[.]\d+',     fn).group())),
			'e':     float(re.sub('_e',     '', re.search('_e[-]?\d+[.]\d+',     fn).group())),
			'fermi': float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fn).group())),
			'gap':   float(re.sub('_gap',   '', re.search('_gap[-]?\d+[.]\d+',   fn).group())),
		}

		return fn_dict

	def GroundOnly_(self, fn_list):
		params = ['type', 'JU', 'N', 'U', 'e']

		data = np.zeros(len(params))
		for fn in fn_list:
			fn_dict = self.ReadFn_(fn)
			data = np.vstack((data, [fn_dict[p] for p in params]))
		data = np.delete(data, 0, axis=0)

		df = pd.DataFrame(data, columns=params)
		df = df.sort_values(by=['JU', 'N', 'U', 'e'])
		df = df.drop_duplicates(subset=['JU', 'N', 'U', 'type'], keep='first')
		grd_list = df.index.to_list()

		return grd_list

	def GenFileList(self, dtype):
		t0 = time.time()

		pn = '%s/params_%s.txt' % (self.path_save, dtype)
		bn = '%s/list_b_%s.txt' % (self.path_save, dtype)
		un = '%s/list_u_%s.txt' % (self.path_save, dtype)

		fn_list = ['%s/%s/%s' % (self.path_output, d, f)\
				for d in os.listdir(self.path_output)   if re.search('JU',    d)\
				for f in os.listdir(self.path_output+d) if re.search('band_', f)]
		fn_list = [fn_list[i] for i in self.GroundOnly_(fn_list)]

		if   dtype == 'M': fn_list = [fn for fn in fn_list if self.ReadFn_(fn)['m']   > 1e-1]
		elif dtype == 'G': fn_list = [fn for fn in fn_list if self.ReadFn_(fn)['gap'] > 1e-1]

		with open(pn, 'w') as f:
			for p in self.params: f.write('%16s' % p)
			f.write('\n')
			for fn in fn_list:
				for p in self.params:
					f.write('%16.10f' % self.ReadFn_(fn)[p])
				f.write('\n')

		with open(bn, 'w') as f:
			f.write('%d\n' % len(fn_list))
			f.write('\n'.join(['%s\t%.16f' % (fn, self.ReadFn_(fn)['fermi']) for fn in fn_list]))

		with open(un, 'w') as f:
			f.write('%d\n' % len(fn_list))
			f.write('\n'.join([re.sub('band', 'ufw', fn) for fn in fn_list]))

		t1 = time.time()
		print('GenFileList(%s, %s, %s) : %fs' % (pn, bn, un, t1-t0))



