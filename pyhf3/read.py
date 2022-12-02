# pyhf3/read.py : read something

import re
import os
import numpy as np
import pandas as pd

def ReadInfo(input_path):
	f_path = open('%s/path.txt' % (input_path), 'r')
	f_cell = open('%s/cell.txt' % (input_path), 'r')

	path = f_path.readline()
	path_point = re.findall('\d+', path)
	path_label = re.findall('[A-Z]+', path)
	info_path = [[int(point), label] for point, label in zip(path_point, path_label)]
	f_path.close()

	info_cell = {}
	while 1:
		line = f_cell.readline()
		if not line: break
		lines = line.strip().split(' ')
		info_cell[lines[0]] = [int(lines[3]), int(lines[4])]
	f_cell.close()

	return info_path, info_cell

def ReadFn(fn, dtype='band'):
	JU    = re.sub('JU',           '', re.search('JU\d+[.]\d+',               fn).group())	
	SOC   = re.sub('_SOC',         '', re.search('_SOC\d+[.]\d+',             fn).group())	
	type  = re.sub('%s_' % dtype,  '', re.search('%s_[a-z]+\d?' % dtype,      fn).group())	
	N     = re.sub('_N',           '', re.search('_N\d+[.]\d+',               fn).group())	
	U     = re.sub('_U',           '', re.search('_U\d+[.]\d+',               fn).group())	
	n     = re.sub('_n',           '', re.search('_n\d+[.]\d+',               fn).group())	
	m     = re.sub('_m',           '', re.search('_m[-]?\d+[.]\d+',           fn).group())	
	e     = re.sub('_e',           '', re.search('_e[-]?\d+[.]\d+',           fn).group())	
	fermi = re.sub('_fermi',       '', re.search('_fermi[-]?\d+[.]\d+',       fn).group())	
	dntop = re.sub('_dntop',       '', re.search('_dntop[-]?\d+[.]\d+',       fn).group())	
	gap   = re.sub('_gap',         '', re.search('_gap\d+[.]\d+',             fn).group())	

	info_dict = {'JU': float(JU), 'SOC': float(SOC), 'type': type, 'N': float(N), 'U': float(U),\
			'n': float(n), 'm': float(m), 'e': float(e), 'fermi': float(fermi), 'dntop': float(dntop), 'gap': float(gap)}
	
	return info_dict

def MakeGroundIdx(fn_list, dtype='band', exclude_f=False):
	df = pd.DataFrame()

	for fn in fn_list:
		fn_dict = ReadFn(fn, dtype=dtype)
		data = pd.DataFrame([[fn_dict['type'][0], fn_dict['N'], fn_dict['U'], fn_dict['n'], fn_dict['e']]], columns=['type', 'N', 'U', 'n', 'e'])
		df = pd.concat([df, data], sort=False)

	df = df.reset_index(drop=True)
	if exclude_f: df = df[df['type'] != 'f']
	df = df[abs(df['N'] - df['n']) < 1e-2] # delete unreliable data

	# drop higher energy
	df = df.sort_values(by=['N', 'U', 'e'])
	df = df.drop_duplicates(subset=['N', 'U', 'type'], keep='first')
	idx_list = df.index.to_list()
	
	return idx_list
