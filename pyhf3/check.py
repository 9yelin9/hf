# pyhf3/check.py : check reliability of data

import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Check:
	def Run(self, material, output, JU, SOC, N4ch=6):
		JU = float(JU)
		SOC = float(SOC)
		N4ch = float(N4ch)
		d = '%s/%s/K16_JU%.2f_SOC%.2f/' % (material, output, JU, SOC)	
		fs_list = ['%s%s' % (d, fs) for fs in os.listdir(d) if re.match('band_', fs)]
		markers = ['s', 'o', '^', 'x']

		df_list = []
		type_list = []
		U_list = []
		gap_list = []
		
		for fs in fs_list:
			type  = re.sub('band_', '', re.search('band_[a-z]+', fs).group())	
			N     = float(re.sub('_N', '', re.search('_N\d+[.]\d+', fs).group()))
			U     = float(re.sub('_U', '', re.search('_U\d+[.]\d+', fs).group()))
			fermi = float(re.sub('_fermi', '', re.search('_fermi\d+[.]\d+', fs).group()))

			if type != 'f':
				f = open(fs, 'r')
				data = np.genfromtxt(f)[:, 1:]
				f.close()

				is_cdt = len(data[np.where(np.abs(data - fermi) < 1e-2)])

				if not is_cdt:
					upper_bot = np.min(data[np.where(data - fermi > 0)])
					lower_top = np.max(data[np.where(data - fermi < 0)])

					gap = upper_bot - lower_top

					if np.abs(N - N4ch) < 1e-6:
						type_list.append(type)
						U_list.append(U)
						gap_list.append(gap)

				del data

		df = pd.DataFrame({'type':type_list, 'U':U_list, 'gap':gap_list})
		for type in np.unique(type_list):
			df_list.append(df[df['type'] == type].sort_values(by='U').reset_index(drop=True))

		for i, df in enumerate(df_list):
			plt.plot(df['U'], df['gap'], '--', marker=markers[i], label=df['type'][0])
	
		plt.xlabel('Interaction (U)')
		plt.ylabel('Band gap at G')
		plt.title('%s @ N=%.1f' % (d, N4ch))
		plt.legend()
		plt.grid(True)
		plt.savefig('%s/%s' % (re.sub('output', 'diagram', d), 'check_N%.1f.png' % N4ch))
		plt.show()
