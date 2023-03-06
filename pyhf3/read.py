# pyhf3/read.py : read something

import os
import re
import numpy as np
import pandas as pd

def ReadInfo(name):
	path_input  = 'input/%s/'  % name.split('_')[0]
	path_output = 'output/%s/' % name

	f_path = open('%s/path.txt' % (path_input), 'r')
	f_cell = open('%s/cell.txt' % (path_input), 'r')

	line = f_path.readline()
	points = re.findall('\d+', line)
	labels = re.findall('[A-Z]+', line)
	info_path = [[int(point), label] for point, label in zip(points, labels)]
	f_path.close()

	info_cell = {}
	while 1:
		line = f_cell.readline()
		if not line: break
		lines = line.strip().split(' ')
		info_cell[lines[0]] = [int(lines[3]), int(lines[4])] # 0:type, 3:Ni, 4:Nc
	f_cell.close()

	return path_input, path_output, info_path, info_cell

def ReadFn(fn):
	fn_dict = {
		'type':  re.sub('band_', '', re.search('band_[a-z]+\d?', fn).group()),
		'JU':    float(re.sub('JU',          '', re.search('JU\d+[.]\d+',        fn).group())),
		'SOC':   float(re.sub('SOC',         '', re.search('SOC\d+[.]\d+',       fn).group())),
		'N':     float(re.sub('_N',          '', re.search('_N\d+[.]\d+',        fn).group())),
		'U':     float(re.sub('_U',          '', re.search('_U\d+[.]\d+',        fn).group())),
		'n':     float(re.sub('_n',          '', re.search('_n\d+[.]\d+',        fn).group())),
		'm':     float(re.sub('_m',          '', re.search('_m[-]?\d+[.]\d+',    fn).group())),
		'e':     float(re.sub('_e',          '', re.search('_e[-]?\d+[.]\d+',    fn).group())),
		'fermi': float(re.sub('fermi',       '', re.search('fermi[-]?\d+[.]\d+', fn).group())), 
		'dntop': float(re.sub('dntop',       '', re.search('dntop[-]?\d+[.]\d+', fn).group())), 
		'gap':   float(re.sub('gap',         '', re.search('gap[-]?\d+[.]\d+',   fn).group())), 
	}
	
	return fn_dict

def GenGroundIdx(fn_list, exclude_f=False):
	params = ['type', 'JU', 'N', 'U', 'n', 'e']

	data = np.zeros(len(params))
	for fn in fn_list:
		fn_dict = ReadFn(fn)
		fn_dict['type'] = ord(fn_dict['type'][0])
		data = np.vstack((data, [fn_dict[p] for p in params]))
	data = np.delete(data, 0, axis=0)

	df = pd.DataFrame(data, columns=params)
	if exclude_f: df = df[df['type'] != ord('f')]
	#df = df[abs(df['N'] - df['n']) < 1e-2] # delete unreliable data

	# drop higher energy
	df = df.sort_values(by=['JU', 'N', 'U', 'e'])
	df = df.drop_duplicates(subset=['JU', 'N', 'U', 'type'], keep='first')
	idx_list = df.index.to_list()
	
	return idx_list
