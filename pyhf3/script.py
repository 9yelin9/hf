# pyhf3/script.py : make job script

import re
import os
import numpy as np
import pandas as pd

def MakeScript(input_path, queue):
	f = open('%s/job/path.txt' % (input_path), 'r')
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
