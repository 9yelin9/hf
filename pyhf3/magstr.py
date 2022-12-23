# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import re
import os
import time
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from .read import ReadFn, MakeGroundIdx

class MagStr:
	def __init__(self, path_output, info_path, info_cell):
		self.path_output = path_output

		info_path_nodup = [] # drop duplicated hsp
		for point, label in info_path:
			if label not in [label for point, label in info_path_nodup]:
				info_path_nodup.append([point, label])

		self.info_path = info_path_nodup
		self.info_cell = info_cell

		self.Ni = self.info_cell['a'][0]
		self.Nc = self.info_cell['a'][1]
		self.Nb = self.Ni * self.Nc * 2

		self.tol_U = 0.1
		self.params = ['JU', 'type', 'N', 'U', 'm', 'fermi', 'dntop', 'gap']

		os.makedirs('%s/magstr/' % path_output, exist_ok=True)

	def Points0(self, pnum): # high symmetry points
		points = []
		labels = []

		for point, label in self.info_path:	
			points.append(point)
			labels.append(label)

		return points, labels

	def Points1(self, pnum): # hsp + around hsp
		points = []
		labels = []
		hnum = -(pnum // 2)

		for point, label in self.info_path:	
			for i in range(pnum + 1):
				points.append(point + hnum + i)
				labels.append(label + str(hnum + i))

		for _ in range(-hnum):
			points = np.delete(points,  0)
			points = np.delete(points, -1)
			labels = np.delete(labels,  0)
			labels = np.delete(labels, -1)

		return points, labels

	def Points2(self, pnum): # hsp + between hsp
		labels = []
		hnum = -(pnum // 2)

		ps = []
		for p, label in self.info_path:
			ps.append(p)
			for i in range(pnum + 1):
				labels.append(label + str(hnum + i))

		itvs = [np.linspace(ps[i], ps[i+1], pnum+2, dtype=int) for i in range(len(ps)-1)]

		points = [itv[:-1] for itv in itvs]
		points = list(np.ravel(points))
		points.append(ps[-1])

		for _ in range(-hnum):
			labels = np.delete(labels,  0)
			labels = np.delete(labels, -1)

		return points, labels

	def MakeBand(self, ptype, pnum):
		t0 = time.time()

		p_dict = {
			'0': self.Points0,
			'1': self.Points1,
			'2': self.Points2,
		}

		bn = '%s/magstr/band_pt%s_pn%s.csv' % (self.path_output, ptype, pnum)
		points, labels = p_dict[ptype](int(pnum))

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
		df.to_csv(bn, sep=',')

		t1 = time.time()
		print('MakeBand(%s) : %fs' % (bn, t1-t0))

	def OpenBand(self, dts, bins):
		bins = int(bins)

		df = pd.read_csv('%s/magstr/band_pt%s_pn%s.csv' % (self.path_output, dts[1], dts[2]), index_col=0)
		df = df[df['U'] > self.tol_U]

		e_label = [v for v in df.columns if re.match('e', v)]
		w_label = [re.sub('e', 'w', v) for v in e_label]

		for e in e_label: df[e] = df[e] - df['dntop']
		df = df.reset_index(drop=True)

		e_list = df.loc[:, e_label].values.flatten()
		e_min = min(e_list)
		e_max = max(e_list)
		e_range = np.linspace(e_min, e_max, bins)

		band = {
			'df': df,
			'e_label': e_label,
			'w_label': w_label,
			'e_range': e_range,
		}

		return band
	
	def MakeDOS0(self, dt, bins, eta, tol, band): # partial dos
		df0 = band['df']
		e_label = band['e_label']
		w_label = band['w_label']
		e_range = band['e_range']

		f_label = ['x%d' % i for i in range(bins)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		for i in df0.index:
			dos_list = list(df0.loc[i, self.params])
			for e in e_range:
				dos = 0
				for el, wl in zip(e_label, w_label):
					dos += (eta / ((e - df0.loc[i, el])**2 + eta**2)) * df0.loc[i, wl]
				dos_list.append(dos / np.pi)
			data = pd.DataFrame([dos_list], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		return df

	def MakeDOS1(self, dt, bins, eta, tol, band): # hsp dos
		df0 = band['df']
		e_label = band['e_label']
		w_label = band['w_label']
		e_range = band['e_range']

		# options
		if re.search('m', dt): df0 = df0[df0['m'] < tol]
		else:                  df0 = df0[df0['m'] > tol]

		if re.search('f', dt): w_fermi = [1 if e < 0 else 0 for e in e_range] 
		else:                  w_fermi = [1 for e in e_range]

		if re.search('b', dt): eta_brden = [0.1 + e/e_range.min() if e < 0 else 0.1 + e/e_range.max() for e in e_range]
		else:                  eta_brden = [eta for e in e_range]

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

		for i in df0.index:
			dos_list = list(df0.loc[i, self.params])
			for hsp in range(n_hsp):
				for e, wf, etab in zip(e_range, w_fermi, eta_brden):
					dos = 0
					for el, wl in zip(e_label_sp[hsp], w_label_sp[hsp]):
						dos += (etab / ((e - df0.loc[i, el])**2 + etab**2)) * df0.loc[i, wl] * wf
					dos_list.append(dos / np.pi)
			data = pd.DataFrame([dos_list], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		return df
	
	def MakeDOS2(self, dt, bins, eta, tol, band): # full dos
		df0 = band['df']
		e_label = band['e_label']
		w_label = band['w_label']
		e_range = band['e_range']

		f_label = ['x%d' % i for i in range(bins-1)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		dir_list = [self.path_output + dir for dir in os.listdir(self.path_output) if not re.search('magstr', dir)]	

		for dir in dir_list:
			fn_list = [dir +'/'+ fn for fn in os.listdir(dir) if re.match('dos', fn)]
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
					x = [e - fn_dict['dntop'] for e in d.index]
					h, _ = np.histogram(x, bins=e_range, weights=list(d))
					dos_list += list(h)

					data = pd.DataFrame([dos_list], columns=x_label)
					df = pd.concat([df, data], sort=False)

		return df

	def MakeDOS(self, dtype, bins, eta, tol):
		t0 = time.time()

		dts = dtype.split(sep='-') # dtype = 'dtype-ptype-pnum'
		tol = float(tol)
		bins = int(bins)
		eta = float(eta)

		dn = '%s/magstr/dos_dt%s_bins%d_eta%.3f_tol%.3f.csv' % (self.path_output, dtype, bins, eta, tol)

		d_dict = {
			'0': self.MakeDOS0,
			'1': self.MakeDOS1,
			'2': self.MakeDOS2,
		}

		band = self.OpenBand(dts, bins)
		df = d_dict[dts[0][0]](dts[0][1:], bins, eta, tol, band)
		df = df.reset_index(drop=True)
		df.to_csv(dn, sep=',')

		t1 = time.time()
		print('MakeDOS(%s) : %fs' % (dn, t1-t0))

	def Resample(self, rspn, X, y):
		from imblearn.under_sampling import RandomUnderSampler
		from imblearn.under_sampling import NearMiss
		from imblearn.under_sampling import EditedNearestNeighbours
		from imblearn.under_sampling import CondensedNearestNeighbour

		rsp_dict = {
			'rus': RandomUnderSampler(random_state=1),
			'nm':  NearMiss(version=1),
			'enn': EditedNearestNeighbours(),
			'cnn': CondensedNearestNeighbour(random_state=1),
		}

		rsp = rsp_dict[rspn]

		X_rsp, y_rsp = rsp.fit_resample(X, y)
		X_rsp.index = X.index[rsp.sample_indices_]
		y_rsp.index = y.index[rsp.sample_indices_]

		return X_rsp, y_rsp

	def LearnOrTune(self, dtype, bins, eta, tol, mcn, rspn):
		from sklearn.ensemble import RandomForestClassifier
		from sklearn.linear_model import LogisticRegression
		from xgboost import XGBClassifier
		from lightgbm import LGBMClassifier
		from catboost import CatBoostClassifier	
		from sklearn.preprocessing import OneHotEncoder
		from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV

		t0 = time.time()

		tol = float(tol)
		bins = int(bins)
		eta = float(eta)
		mcn = mcn.split(sep='_')

		mc_dict = {
			'rf':  RandomForestClassifier(random_state=1),
			'lr': LogisticRegression(random_state=1, multi_class='multinomial', max_iter=10000),
			'xgb': XGBClassifier(random_state=1),
			'lgb': LGBMClassifier(random_state=1),
			'cat': CatBoostClassifier(random_state=1, silent=True),
		}
		hp_dict = {
			'rf': {
				'random_state': 1,
			}
		}
		mc = mc_dict[mcn[0]]

		# prepare data
		dn = '%s/magstr/dos_dt%s_bins%d_eta%.3f_tol%.3f.csv' % (self.path_output, dtype, bins, eta, tol)
		df = pd.read_csv(dn, sep=',', index_col=0)

		X = df.drop(self.params, axis=1)
		y = df['type']
		if mcn[0] == 'xgb': y = pd.get_dummies(y)
		if rspn: X, y = self.Resample(rspn, X, y)
		idxs = y.index.to_list()

		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=1)
		
		# Learn
		if len(mcn) < 2:
			mc.fit(X_train, y_train)

			t1 = time.time()
			print('Learn(%s-%s-%s) : %fs' % (mcn[0], rspn, dn, t1-t0))

			return mc, idxs, X_test, y_test

		# Tune
		else:
			cv = StratifiedKFold(random_state=1)
			tn = RandomizedSearchCV(mc, hp_dict[mcn[0]], cv, random_state=1)
			res = tn.fit(X_train, y_train)

			t1 = time.time()
			print('Tune(%s-%s-%s) : %fs' % (mcn[0], rspn, dn, t1-t0))

			return res
		
