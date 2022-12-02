# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import re
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from .read import ReadFn, MakeGroundIdx

class MagStr:
	def __init__(self, path_output, info_path, info_cell, ptype, pnum=0):
		self.path_output = path_output

		info_path_nodup = [] # drop duplicated hsp
		for point, label in info_path:
			if label not in [label for point, label in info_path_nodup]:
				info_path_nodup.append([point, label])

		self.info_path = info_path_nodup
		self.info_cell = info_cell
		self.ptype = ptype
		self.pnum = int(pnum)

		self.Ni = self.info_cell['a'][0]
		self.Nc = self.info_cell['a'][1]
		self.Nb = self.Ni * self.Nc * 2

		self.tol_U = 1e-1
		self.params = ['JU', 'SOC', 'type', 'N', 'U', 'm', 'fermi', 'dntop', 'gap']

		os.makedirs('%s/magstr/' % path_output, exist_ok=True)

	def Points0(self): # high symmetry points
		points = []
		labels = []

		for point, label in self.info_path:	
			points.append(point)
			labels.append(label)

		return points, labels

	def Points1(self): # hsp + around hsp
		points = []
		labels = []
		hnum = -(self.pnum // 2)

		for point, label in self.info_path:	
			for i in range(self.pnum + 1):
				points.append(point + hnum + i)
				labels.append(label + str(hnum + i))

		for _ in range(-hnum):
			points = np.delete(points,  0)
			points = np.delete(points, -1)
			labels = np.delete(labels,  0)
			labels = np.delete(labels, -1)

		return points, labels

	def Points2(self): # hsp + between hsp
		labels = []
		hnum = -(self.pnum // 2)

		ps = []
		for p, label in self.info_path:
			ps.append(p)
			for i in range(self.pnum + 1):
				labels.append(label + str(hnum + i))

		itvs = [np.linspace(ps[i], ps[i+1], self.pnum+2, dtype=int) for i in range(len(ps)-1)]

		points = [itv[:-1] for itv in itvs]
		points = list(np.ravel(points))
		points.append(ps[-1])

		for _ in range(-hnum):
			labels = np.delete(labels,  0)
			labels = np.delete(labels, -1)

		return points, labels

	def MakeBand(self):
		p_dict = {
			'0': self.Points0,
			'1': self.Points1,
			'2': self.Points2,
		}

		points, labels = p_dict[self.ptype]()

		e_label = []
		w_label = []
		for label in labels:
			e_label += ['e%s_%d' % (label, i) for i in range(self.Nb)]
			w_label += ['w%s_%d' % (label, i) for i in range(self.Nb)]
		x_label = self.params + e_label + w_label

		df = pd.DataFrame()

		dir_list = [self.path_output + dir for dir in os.listdir(self.path_output) if not re.search('magstr', dir)]	

		for dir in dir_list:
			fn_list = [dir +'/'+ fn for fn in os.listdir(dir) if re.match('band_', fn)]
			idx_list = MakeGroundIdx(fn_list, exclude_f=True)

			for i in idx_list:
				v = []

				fn = fn_list[i]
				fn_dict = ReadFn(fn)
				fn_dict['type'] = fn_dict['type'][0]
				for p in self.params: v.append(fn_dict[p])
				
				fe = open(fn, 'r')
				fw = open(re.sub('band_', 'ufw_', fn), 'r')
				de = np.genfromtxt(fe)
				dw = np.genfromtxt(fw)
				fe.close()
				fw.close()

				for point in points: v += list(de[point][:])
				for point in points: v += list(dw[point][:])

				data = pd.DataFrame([v], columns=x_label)
				df = pd.concat([df, data], sort=False)

		df = df.reset_index(drop=True)
		df.to_csv('%s/magstr/band_ptype%s_pnum%d.csv' % (self.path_output, self.ptype, self.pnum), sep=',')

	def OpenBand(self, tol, under_tol):
		tol = float(tol)

		df = pd.read_csv('%s/magstr/band_ptype%s_pnum%d.csv' % (self.path_output, self.ptype, self.pnum), index_col=0)

		df = df[df['U'] > self.tol_U]
		if under_tol != '0': df = df[df['m'] < tol]
		else: df = df[df['m'] > tol]
		df = df.reset_index(drop=True)

		e_label = [col for col in df.columns if re.match('e', col)]
		w_label = [col for col in df.columns if re.match('w', col)]

		for e in e_label: df[e] = df[e] - df['dntop']

		return df, e_label, w_label

	def CalcEnergyMinMax(self, df, e_label):
		e_list = df.loc[:, e_label].values.flatten()
		e_min = min(e_list)
		e_max = max(e_list)

		return e_min, e_max

	def MakeDOS(self, dtype, tol, bins, eta, under_tol='0', broadening='0'):
		tol = float(tol)
		bins = int(bins)
		eta = float(eta)

		d_dict = {
			'0': self.MakeDOS0,
			'1': self.MakeDOS1,
			'2': self.MakeDOS2,
		}

		fn_dos = '%s/magstr/dos_ptype%s_pnum%d_dtype%s_tol%.3f_bins%d_eta%.3f.csv' % (self.path_output, self.ptype, self.pnum, dtype, tol, bins, eta)
		if under_tol  != '0': fn_dos = re.sub('dos', 'dos_und', fn_dos)
		if broadening != '0': fn_dos = re.sub('dos', 'dos_bro', fn_dos)

		d_dict[dtype](tol, bins, eta, under_tol, broadening, fn_dos)

	def MakeDOS0(self, tol, bins, eta, under_tol, broadening, fn_dos):
		df0, e_label, w_label = self.OpenBand(tol, under_tol)
		e_min, e_max = self.CalcEnergyMinMax(df0, e_label)
		e_range = np.linspace(e_min, e_max, bins)

		f_label = ['x%d' % i for i in range(bins)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		for i in df0.index.to_list():
			dos_list = list(df0.loc[i, self.params])
			for e in e_range:
				dos = 0
				for el, wl in zip(e_label, w_label):
					dos += (eta / ((e - df0.loc[i, el])**2 + eta**2)) * df0.loc[i, wl]
				dos_list.append(dos / np.pi)
			data = pd.DataFrame([dos_list], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		df.to_csv(fn_dos, sep=',')

	def MakeDOS1(self, tol, bins, eta, under_tol, broadening, fn_dos):
		df0, e_label, w_label = self.OpenBand(tol, under_tol)
		e_min, e_max = self.CalcEnergyMinMax(df0, e_label)
		e_range = np.linspace(e_min, e_max, bins)
		n_hsp = len(e_label) // self.Nb

		f_label = []
		e_label_sp = [[] for _ in range(n_hsp)]
		w_label_sp = [[] for _ in range(n_hsp)]

		for i in range(n_hsp):
			for j in range(bins):
				f_label.append('%s_%d' % (re.sub('e', '', e_label[self.Nb*i].split('_')[0]), j))
			for j in range(self.Nb):
				e_label_sp[i].append(e_label[self.Nb*i + j])
				w_label_sp[i].append(w_label[self.Nb*i + j])
		x_label = self.params + f_label

		df = pd.DataFrame()

		for i in df0.index.to_list():
			dos_list = list(df0.loc[i, self.params])
			for hsp in range(n_hsp):
				for e in e_range:
					dos = 0
					for el, wl in zip(e_label_sp[hsp], w_label_sp[hsp]):
						dos += (eta / ((e - df0.loc[i, el])**2 + eta**2)) * df0.loc[i, wl]
					dos_list.append(dos / np.pi)
			data = pd.DataFrame([dos_list], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		df.to_csv(fn_dos, sep=',')
	
	def MakeDOS2(self, tol, bins, eta, under_tol, broadening, fn_dos):
		df0, e_label, _ = self.OpenBand(tol, under_tol)
		e_min, e_max = self.CalcEnergyMinMax(df0, e_label)
		e_range = np.linspace(e_min, e_max, bins+1)

		f_label = ['x%d' % i for i in range(bins)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		dir_list = [self.path_output + dir for dir in os.listdir(self.path_output) if not re.search('magstr', dir)]	

		for dir in dir_list:
			fn_list = [dir +'/'+ fn for fn in os.listdir(dir) if re.match('dos_', fn)]
			idx_list = MakeGroundIdx(fn_list, dtype='dos', exclude_f=True)

			for i in idx_list:
				dos_list = []

				fn = fn_list[i]
				fn_dict = ReadFn(fn, dtype='dos')
				fn_dict['type'] = fn_dict['type'][0]
				for p in self.params: dos_list.append(fn_dict[p])

				if fn_dict['U'] > self.tol_U and fn_dict['m'] > tol:
					f = open(fn, 'r')
					d = np.genfromtxt(f)
					f.close()

					d = np.transpose(d)
					d = pd.DataFrame(d[1:, :], columns=d[0, :]).sum()
					x = [e - fn_dict['dntop'] for e in d.index.to_list()]
					h, _ = np.histogram(x, bins=e_range, weights=list(d))
					dos_list += list(h)

					data = pd.DataFrame([dos_list], columns=x_label)
					df = pd.concat([df, data], sort=False)

		df = df.reset_index(drop=True)
		df.to_csv(fn_dos, sep=',')

	def DoRandomForest(self, dtype, tol, bins, eta, rsp='', test_size=0.3):
		from sklearn.ensemble import RandomForestClassifier

		tol = float(tol)
		bins = int(bins)
		eta = float(eta)
		test_size = float(test_size)

		df = pd.read_csv('%s/magstr/dos_ptype%s_pnum%d_dtype%s_tol%.3f_bins%d_eta%.3f.csv' % (self.path_output, self.ptype, self.pnum, dtype, tol, bins, eta), sep=',', index_col=0)
		X = df.drop(self.params, axis=1)
		y = df['type']

		if rsp != '':
			from imblearn.under_sampling import RandomUnderSampler
			from imblearn.under_sampling import NearMiss
			from imblearn.under_sampling import EditedNearestNeighbours
			from imblearn.under_sampling import CondensedNearestNeighbour

			if   rsp == 'rus': rsp = RandomUnderSampler(random_state=1)
			elif rsp == 'nm':  rsp = NearMiss(version=1)
			elif rsp == 'enn': rsp = EditedNearestNeighbours()
			elif rsp == 'cnn': rsp = CondensedNearestNeighbour(random_state=1)
			else: 
				print('%s is wrong resampler\n' % rsp)
				sys.exit()

			X_rsp, y_rsp = rsp.fit_resample(X, y)
			X_rsp.index = X.index[rsp.sample_indices_]
			y_rsp.index = y.index[rsp.sample_indices_]

			X = X_rsp
			y = y_rsp

		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=1)

		rf = RandomForestClassifier(criterion='gini', random_state=1)
		rf.fit(X_train, y_train)
		y_pred = rf.predict(X_test)
		y_score = rf.predict_proba(X_test)
		acc = accuracy_score(y_test, y_pred)
		
		return acc, y_test, y_pred, y_score, y.index.to_list(), rf.estimators_, rf.feature_importances_

	def DoLogisticRegression(self):
		from sklearn.preprocessing import StandardScaler
		from sklearn.linear_model import LogisticRegression

		return 0
