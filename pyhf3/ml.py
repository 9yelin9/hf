# pyhf3/ml.py : make csv file for machine learning

import re
import os
import numpy as np

class ML:
	def Run(self, material, basis, path, path_label):
		f_highsym = open('%s/output/highsym.csv' % (material), 'w')
		f_highsym.write('K,JU,SOC,type,N,U,fermi')
		for pl in path_label:
			for i in range(basis[1]): f_highsym.write(',w%s%d' % (pl, i))
			for i in range(basis[1]): f_highsym.write(',e%s%d' % (pl, i))
		f_highsym.write('\n')

		dlist = ['%s/output/%s' % (material, d) for d in os.listdir('%s/output' % (material)) if re.search('_unfold', d)]	
		
		for d in dlist:
			for fs in ['%s/%s' % (d, fs) for fs in os.listdir(d) if re.match('band_', fs)]:
				K     = re.sub('K', '', re.search('K\d+', d).group())	
				JU    = re.sub('_JU', '', re.search('_JU\d+[.]\d+', d).group())	
				SOC   = re.sub('_SOC', '', re.search('_SOC\d+[.]\d+', d).group())	
				type  = re.sub('band_', '', re.search('band_[a-z]+', fs).group())	
				N     = re.sub('_N', '', re.search('_N\d+[.]\d+', fs).group())	
				U     = re.sub('_U', '', re.search('_U\d+[.]\d+', fs).group())	
				fermi = re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group())	
				f_highsym.write('%s,%s,%s,%s,%s,%s,%s' % (K, JU, SOC, type, N, U, fermi))

				f = open(fs, 'r')
				data = np.genfromtxt(f)
				f.close()
				
				for p in path:
					for v in data[p][1:]: f_highsym.write(',%s' % (v))
				f_highsym.write('\n')
				del data

		f_highsym.close()
