# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import re
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from .read import ReadFs, MakeGroundIdx

class ML:
	def __init__(self, output_path, info_path, info_cell, ptype):
		self.output_path = output_path
		self.path_point = [path_point for path_point, path_label in info_path]
		self.path_label = [path_label for path_point, path_label in info_path]
		self.info_cell = info_cell
		self.ptype = int(ptype)

		self.Nb = self.info_cell['a']

		os.makedirs('%s/magstr/' % output_path, exist_ok=True)

	def PointType0(self, num): # high symmetry points only
		points = []
		labels = []

		for path_point, path_label in zip(self.path_point, self.path_label):	
			points.append(path_point)
			labels.append(path_label)

		return points, labels

	def PointType1(self, num): # hsp + around hsp
		points = []
		labels = []
		hnum = num // 2

		for i, (path_point, path_label) in enumerate(zip(self.path_point, self.path_label)):	
			if i == 0:
				n = hnum + 1
				v = 0
			elif i == len(self.info_path) - 1:
				n = hnum + 1
				v = -hnum
			else:
				n = num + 1
				v = -hnum

			for j in range(n):
				points.append(path_point + v + j)
				labels.append('%s%d' % (path_label, v + j))

		return points, labels

	def PointType2(self, num): # hsp + same interval between hsp
		points = []
		labels = []
		hnum = num // 2
		itvs = [np.linspace(self.path_point[i], self.path_point[i+1], num+2, dtype=int) for i in range(len(self.path_point) - 1)]

		for i, itv in enumerate(itvs):
			if i != len(itvs)-1:
				itv = np.delete(itv, -1)
			for j in range(len(itv)):
				points.append(itv[j])

		for i, (path_point, path_label) in enumerate(zip(self.path_point, self.path_label)):	
			if i == 0:
				n = hnum + 1
				v = 0
			elif i == len(self.info_path) - 1:
				n = hnum + 1
				v = -hnum
			else:
				n = num + 1
				v = -hnum

			for j in range(n):
				labels.append('%s%d' % (path_label, v + j))

		return points, labels

	def MakeBand(self, tol, num):
		num = int(num)

		f = open('%s/magstr/band_ptype%d_num%d.csv' % (self.output_path, self.ptype, num), 'w')

		if   self.ptype == 0: points, labels = self.PointType0(num)
		elif self.ptype == 1: points, labels = self.PointType1(num)
		elif self.ptype == 2: points, labels = self.PointType2(num)
		else:
			points = []
			labels = []

		f.write('JU,SOC,type,N,U,m,fermi,dntop,gap')
		for label in labels:
			for i in range(self.Nb):
				f.write(',e%s_%d' % (label, i))
			for i in range(self.Nb):
				f.write(',w%s_%d' % (label, i))
		f.write('\n')

		dlist = [self.output_path + d for d in os.listdir(self.output_path) if not re.search('magstr', d)]	
		self.info_path[-1][0] -= 1

		for d in dlist:
			fs_list = [d +'/'+ fs for fs in os.listdir(d) if re.match('band_', fs)]
			idx_list = MakeGroundIdx(fs_list, True)

			for i in idx_list:
				fs = fs_list[i]
				fs_dict = ReadFs(fs)

				f.write('%.2f,%.2f,%s,%.1f,%.1f,%f,%f,%f,%f'\
						% (fs_dict['JU'], fs_dict['SOC'], fs_dict['type'][0], fs_dict['N'], fs_dict['U'], fs_dict['m'], fs_dict['fermi'], fs_dict['dntop'], fs_dict['gap']))
				
				fe = open(fs, 'r')
				fw = open(re.sub('band_', 'ufw_', fs), 'r')
				datae = np.genfromtxt(fe)
				dataw = np.genfromtxt(fw)
				fe.close()
				fw.close()

				for point in points:
					for e in datae[point][:]:
						f.write(',%s' % (e))
					for w in dataw[point][:]:
						f.write(',%s' % (w))
				f.write('\n')
		f.close()

	def OpenBand(self):
		df = pd.read_csv('%s/magstr/band_highsym.csv' % (self.output_path, self.ptype))

		dup_cols = [col for col in df.columns if re.search('[.]', col)]
		df = df.drop(df.loc[:, dup_cols], axis=1)
		df = df[df['U'] > 1e-1]
		df = df[df['m'] > 1e-3]
		df = df.reset_index(drop=True)

		e_label = [col for col in df.columns if re.match('e', col)]
		w_label = [col for col in df.columns if re.match('w', col)]
		p_label = [col for col in df.columns if col not in (e_label + w_label)]

		for e in e_label: df.loc[:, e] = df.loc[:, e] - df.loc[:, 'dntop']

		return df, e_label, w_label, p_label

	def CalcEnergyMinMax(self, df, e_label):
		e_list = df.loc[:, e_label].values.flatten()
		e_min = min(e_list)
		e_max = max(e_list)

		return e_min, e_max

	def MakeDOS(self, bins, eta):
		bins = int(bins)
		eta = float(eta)

		df0, e_label, w_label, p_label = self.OpenBand()
		e_min, e_max = self.CalcEnergyMinMax(df0, e_label)
		e_range = np.linspace(e_min, e_max, bins)
		x_label = ['x%d' % i for i in range(bins)]
		x_label.insert(0, 'type')

		df = pd.DataFrame()

		for i in range(len(df0)):
			dos_list = []
			for e in e_range:
				dos = 0
				for el, wl in zip(e_label, w_label):
					dos += (eta / ((e - df0.loc[i, el])**2 + eta**2)) * df0.loc[i, wl] / np.pi
				dos_list.append(dos)
			dos_list.insert(0, df0.loc[i, 'type'])
			data = pd.DataFrame([dos_list], columns=x_label)
			df = pd.concat([df, data], sort=False)

		df = df.reset_index(drop=True)
		df.to_csv('%s/magstr/dos_%s_bins%d_eta%.2f.csv' % (self.output_path, self.ptype, bins, eta), sep=',')

	def DoRandomForest(self, bins, eta, test_size=0.3):
		bins = int(bins)
		eta = float(eta)
		test_size = float(test_size)

		df = pd.read_csv('%s/magstr/dos_%s_bins%d_eta%.2f.csv' % (self.output_path, self.ptype, bins, eta))
		X = df.drop(['type'], axis=1)
		y = df['type']
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=1)

		rf = RandomForestClassifier(criterion='gini', random_state=1)
		rf.fit(X_train, y_train)
		y_pred = rf.predict(X_test)
		acc = accuracy_score(y_test, y_pred)
		
		return acc, y_test, y_pred
