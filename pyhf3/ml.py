# pyhf3/ml.py : make csv file for machine learning

import re
import os
import numpy as np
import pandas as pd
from .read import ReadFs

def MakeMLData(output_path, path_info, type_info, tol):
	f = open('%s/highsym.csv' % output_path, 'w')

	f.write('JU,SOC,type,N,U,dntop,gap')
	for path_point, path_label in path_info:
		for i in range(type_info['a']):
			f.write(',e%s%d' % (path_label, i))
		for i in range(type_info['a']):
			f.write(',w%s%d' % (path_label, i))
	f.write('\n')

	dlist = [output_path + d for d in os.listdir(output_path) if os.path.isdir(output_path + d)]	
	path_info[-1][0] -= 1

	for d in dlist:
		fs_list = [d +'/'+ fs for fs in os.listdir(d) if re.match('band_', fs)]
		df = pd.DataFrame()

		for i, fs in enumerate(fs_list):
			fs_dict = ReadFs(fs)
			data = pd.DataFrame([[fs_dict['type'][0], fs_dict['N'], fs_dict['U'], fs_dict['e']]], columns=['type', 'N', 'U', 'e'])
			df = pd.concat([df, data], sort=False)

		# drop higher energy
		df = df.reset_index(drop=True)
		df = df[df['type'] != 'f']
		df = df.sort_values(by='e')
		df = df.drop_duplicates(subset=['type', 'N', 'U'], keep='first')
		idx_list = df.index.to_list()

		for i in idx_list:
			fs = fs_list[i]
			fs_dict = ReadFs(fs)

			f.write('%.2f,%.2f,%s,%.1f,%.1f,%f,%f'\
					% (fs_dict['JU'], fs_dict['SOC'], fs_dict['type'][0], fs_dict['N'], fs_dict['U'], fs_dict['dntop'], fs_dict['gap']))
			
			fe = open(fs, 'r')
			fw = open(re.sub('band_', 'ufw_', fs), 'r')
			datae = np.genfromtxt(fe)
			dataw = np.genfromtxt(fw)
			fe.close()
			fw.close()

			for path_point, path_label in path_info:
				for e in datae[path_point][1:]:
					f.write(',%s' % (e))
				for w in dataw[path_point][1:]:
					f.write(',%s' % (w))
			f.write('\n')
	f.close()
