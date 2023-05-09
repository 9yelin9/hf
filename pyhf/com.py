# pyhf/com.py : common functions

import re
import sys
import numpy as np
import pandas as pd

def ReadConfig(dim):
	Nkg1 = 0
	Nkg  = 0
	Nkb  = 0

	fn = 'input/config.txt'
	with open(fn, 'r') as f:
		for line in f:
			if   re.search('Nkg1', line): Nkg1 = int(re.sub('Nkg1', '', line))
			elif re.search('Nkb',  line): Nkb  = int(re.sub('Nkb',  '', line))
	Nkg = Nkg1 ** dim

	if not (Nkg1 or Nkg or Nkb):
		print('Wrong %s\nNkg1 = %d, Nkg = %d, Nkb = %d' % (fn, Nkg1, Nkg, Nkb))
		sys.exit()

	return Nkg1, Nkg, Nkb

def FnDict(fn):
	fn_dict = {
		'type':  re.sub('_', '', re.search('[A-Z]\d_', fn).group()),
		'JU':    float(re.sub('JU',    '', re.search('JU\d+[.]\d+',        fn).group())),
		'N':     float(re.sub('N',     '', re.search('N\d+[.]\d+',         fn).group())),
		'U':     float(re.sub('_U',    '', re.search('_U\d+[.]\d+',        fn).group())),
		'n':     float(re.sub('_n',    '', re.search('_n\d+[.]\d+',        fn).group())),
		'm':     float(re.sub('_m',    '', re.search('_m[-]?\d+[.]\d+',    fn).group())),
		'e':     float(re.sub('_e',    '', re.search('_e[-]?\d+[.]\d+',    fn).group())),
		'fermi': float(re.sub('fermi', '', re.search('fermi[-]?\d+[.]\d+', fn).group())),
		'gap':   float(re.sub('gap',   '', re.search('gap[-]?\d+[.]\d+',   fn).group())),
	}

	return fn_dict

def GroundOnly(fn_list, n_list, u_list):
	params = ['type', 'N', 'U', 'e']

	grd_fn_list = []

	for n in n_list:
		for u in u_list:
			e_list = []
			grd_fn = []
			for fn in fn_list:
				if re.search('N%.1f_U%.1f' % (n, u), fn):
					e_list.append(np.around(FnDict(fn)['e']))
					grd_fn.append(fn)
			grd_fn_list.append(grd_fn[np.argmin(e_list)])

	return grd_fn_list
