# pyhf/mod.py : modules

import re
import numpy as np
import pandas as pd

def FnDict(fn):
	fn_dict = {
		'type':  re.sub('_', '', re.search('[A-Z]\d_', fn).group()),
		'JU':    float(re.sub('_JU',    '', re.search('_JU\d+[.]\d+',        fn).group())),
		'N':     float(re.sub('_N',     '', re.search('_N\d+[.]\d+',         fn).group())),
		'U':     float(re.sub('_U',     '', re.search('_U\d+[.]\d+',         fn).group())),
		'n':     float(re.sub('_n',     '', re.search('_n\d+[.]\d+',         fn).group())),
		'm':     float(re.sub('_m',     '', re.search('_m[-]?\d+[.]\d+',     fn).group())),
		'e':     float(re.sub('_e',     '', re.search('_e[-]?\d+[.]\d+',     fn).group())),
		'fermi': float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fn).group())),
		'gap':   float(re.sub('_gap',   '', re.search('_gap[-]?\d+[.]\d+',   fn).group())),
	}

	return fn_dict

def GroundOnly(fn_list):
	params = ['type', 'JU', 'N', 'U', 'e']

	data = np.zeros(len(params))
	for fn in fn_list: data = np.vstack((data, [FnDict(fn)[p] for p in params]))
	data = np.delete(data, 0, axis=0)

	df = pd.DataFrame(data, columns=params)
	df = df.sort_values(by=['JU', 'N', 'U', 'e'])
	df = df.drop_duplicates(subset=['JU', 'N', 'U', 'type'], keep='first')
	grd_idx = df.index.to_list()

	return grd_idx
