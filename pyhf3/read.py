# pyhf3/read.py : read something

import re
import os
import numpy as np

def ReadInfo(input_path):
	f = open('%s/info.txt' % (input_path), 'r')

	path_len = int(f.readline())
	path_point = f.readline().strip().split(' ')
	path_label = f.readline().strip().split(' ')
	path_info = [[int(point), label] for point, label in zip(path_point, path_label)]

	type_info = {}
	line = f.readline()
	while 1:
		line = f.readline()
		if not line: break
		lines = line.strip().split(' ')
		type_info[lines[0]] = lines[1]

	f.close()

	return path_info, type_info

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
