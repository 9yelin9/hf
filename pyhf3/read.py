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

def ReadFs(fs, dtype='band'):
	JU    = re.sub('JU',           '', re.search('JU\d+[.]\d+',               fs).group())	
	SOC   = re.sub('_SOC',         '', re.search('_SOC\d+[.]\d+',             fs).group())	
	type  = re.sub('%s_' % dtype,  '', re.search('%s_[a-z]+\d?' % dtype,      fs).group())	
	N     = re.sub('_N',           '', re.search('_N\d+[.]\d+',               fs).group())	
	U     = re.sub('_U',           '', re.search('_U\d+[.]\d+',               fs).group())	
	n     = re.sub('_n',           '', re.search('_n\d+[.]\d+',               fs).group())	
	m     = re.sub('_m',           '', re.search('_m[-]?\d+[.]\d+',           fs).group())	
	e     = re.sub('_e',           '', re.search('_e[-]?\d+[.]\d+',           fs).group())	
	fermi = re.sub('_fermi',       '', re.search('_fermi[-]?\d+[.]\d+',       fs).group())	
	dntop = re.sub('_dntop',       '', re.search('_dntop[-]?\d+[.]\d+',       fs).group())	
	gap   = re.sub('_gap',         '', re.search('_gap\d+[.]\d+',             fs).group())	

	info_dict = {'JU': float(JU), 'SOC': float(SOC), 'type': type, 'N': float(N), 'U': float(U),\
			'n': float(n), 'm': float(m), 'e': float(e), 'fermi': float(fermi), 'dntop': float(dntop), 'gap': float(gap)}
	
	return info_dict

def MakeGroundIdx(fs_list, exclude_f=False):
	df = pd.DataFrame()

	for fs in fs_list:
		fs_dict = ReadFs(fs)
		data = pd.DataFrame([[fs_dict['type'], fs_dict['N'], fs_dict['U'], fs_dict['n'], fs_dict['e']]], columns=['type', 'N', 'U', 'n', 'e'])
		df = pd.concat([df, data], sort=False)

	df = df.reset_index(drop=True)
	if exclude_f: df = df[df['type'] != 'f']
	df = df[abs(df['N'] - df['n']) < 1e-2] # delete unreliable data

	# drop higher energy
	df = df.sort_values(by=['N', 'U', 'e', 'type'])
	df = df.drop_duplicates(subset=['N', 'U'], keep='first')
	idx_list = df.index.to_list()
	
	return idx_list
