# pyhf3/ml.py : make csv file for machine learning

import re
import os
import numpy as np

class ML:
	def Run(self, material, output, basis, path, path_label):
		f_highsym = open('%s/%s/highsym.csv' % (material, output), 'w')

		f_highsym.write('K,JU,SOC,type,N,U,fermi')
		for pl in path_label:
			for i in range(basis[1]): f_highsym.write(',w%s%d' % (pl, i))
			for i in range(basis[1]): f_highsym.write(',e%s%d' % (pl, i))
		f_highsym.write('\n')

		dlist = ['%s/%s/%s' % (material, output, d) for d in os.listdir('%s/%s' % (material, output))\
				if os.path.isdir('%s/%s/%s' % (material, output, d))]	
		
		for d in dlist:
			for fs in ['%s/%s' % (d, fs) for fs in os.listdir(d) if re.match('band_', fs)]:
				K     = re.sub('K', '', re.search('K\d+', d).group())	
				JU    = re.sub('_JU', '', re.search('_JU\d+[.]\d+', d).group())	
				SOC   = re.sub('_SOC', '', re.search('_SOC\d+[.]\d+', d).group())	
				type  = re.sub('band_', '', re.search('band_[a-z]+', fs).group())	
				N     = re.sub('_N', '', re.search('_N\d+[.]\d+', fs).group())	
				U     = re.sub('_U', '', re.search('_U\d+[.]\d+', fs).group())	
				fermi = re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group())	

				if type != 'f':
					f_highsym.write('%s,%s,%s,%s,%s,%s,%s' % (K, JU, SOC, type, N, U, fermi))

					f_band = open(fs, 'r')
					f_ufw  = open(re.sub('band_', 'ufw_', fs), 'r')

					data_band = np.genfromtxt(f_band)
					data_ufw  = np.genfromtxt(f_ufw)

					f_band.close()
					f_ufw.close()
					
					for p in path:
						for vb in data_band[p][1:]: f_highsym.write(',%s' % (vb))
						for vu in data_ufw[p][1:]:  f_highsym.write(',%s' % (vu))
					f_highsym.write('\n')

					del data_band
					del data_ufw

		f_highsym.close()
