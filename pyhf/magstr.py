# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import os
import re
import sys
import h5py
import time
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
from itertools import product 

from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score, confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV
from sklearn.preprocessing import label_binarize
from sklearn import tree
from sklearn.decomposition import PCA

# classifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier	
from sklearn.svm import SVC

# regressor
from sklearn.linear_model import LinearRegression

# resampler
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NearMiss
from imblearn.under_sampling import EditedNearestNeighbours
from imblearn.under_sampling import CondensedNearestNeighbour

from .read import ReadInfo, ReadFn, GenGroundIdx

class MagStr:
	def __init__(self, name, num_thread=1):
		self.path_input, self.path_output, self.info_path, _ = ReadInfo(name);
		self.path_save = self.path_output + '/magstr/'
		os.makedirs(self.path_save, exist_ok=True)

		info_path_nodup = [] # drop duplicated hsp
		for point, label in self.info_path:
			if label not in [label for point, label in info_path_nodup]:
				info_path_nodup.append([point, label])
		self.points = [point for point, _ in info_path_nodup]
		self.labels = [label for _, label in info_path_nodup]
		#self.Nkb = 1024
		#self.points = range(self.Nkb)
		#self.labels = ['X%d_' for i in self.points]

		self.points_g = [point for point, _ in self.info_path]
		self.labels_g = [label.replace('G', r'$\Gamma$') for _, label in self.info_path]

		#self.info_cell = info_cell
		#self.Ni = self.info_cell['a'][0]
		#self.Nc = self.info_cell['a'][1]
		#self.Nb = self.Ni * self.Nc * 2

		self.Nbs = 12
		self.Nhsp = 4
		self.BINS_MAX = 1024
		self.PEAK_MAX = 2
		self.emin = -8
		self.emax =  8

		self.type_dict   = {'a':1, 'c':2, 'g':3}
		self.type_dict_r = {'1':'a', '2':'c', '3':'g'}

		self.params1 = ['type', 'JU', 'N', 'U']
		self.params2 = ['m', 'gap']
		self.params = self.params1 + self.params2

		self.rsp_dict = {
			'rus': RandomUnderSampler(random_state=1),
			'nm':  NearMiss(version=1),
			'enn': EditedNearestNeighbours(),
			'cnn': CondensedNearestNeighbour(random_state=1),
		}

		random_state = 1
		self.mc_dict = {
			'rf': {
				'func': RandomForestClassifier(random_state=random_state, n_jobs=num_thread),
				'onehot': False,
			},
			'xgb': {
				'func': XGBClassifier(random_state=random_state, nthread=num_thread),
				'onehot': True,
			},
			'lgb': {
				'func': LGBMClassifier(random_state=random_state, n_jobs=num_thread),
				'onehot': False,
			},
			'cat': {
				'func': CatBoostClassifier(random_state=random_state, thread_count=num_thread, silent=True),
				'onehot': False,
			},
			'lr': {
				'func': LogisticRegression(random_state=random_state, n_jobs=num_thread, multi_class='multinomial', max_iter=10000),
				'onehot': False,
			},
			'svm': {
				'func': SVC(random_state=random_state, probability=True),
				'onehot': False,
			},
		}

		self.rc_dict = {
			'lr': LinearRegression(n_jobs=num_thread),
		}

		self.figsize=(11, 5)
		plt.rcParams.update({'font.size': 35})
		plt.rcParams.update({'font.family': 'sans-serif'})
		plt.rcParams.update({'font.serif': 'Helvetica Neue'})
		#plt.rcParams.update({'mathtext.fontset': 'cm'})

	def GenEnergy(self):
		t0 = time.time()

		en = '%s/energy.h5' % self.path_save
		e = np.linspace(self.emin, self.emax, self.BINS_MAX)
		with h5py.File(en, 'w') as f: f.create_dataset('energy', data=e, dtype='d')

		t1 = time.time()
		print('GenEnergy(%s) : %fs' % (en, t1-t0))

	def GenBand(self, btype, dfermi=1e-2):
		pn = '%s/params_%s.txt' % (self.path_save, btype)
		bn = '%s/band_%s.txt' % (self.path_save, btype)
		Fp = open(pn, 'w')
		Fb = open(bn, 'w')

		t0 = time.time()

		with h5py.File('%s/energy.h5' % self.path_save, 'r') as f: energy = f['energy'][()]

		fn_list = ['%s/%s/%s' % (self.path_output, d, f)\
				for d in os.listdir(self.path_output)     if re.search('JU',    d)\
				for f in os.listdir(self.path_output + d) if re.search('band_', f)]
		fn_list = [fn_list[i] for i in GenGroundIdx(fn_list, exclude_f=True)] 
		if   btype == 'M': fn_list = [fn for fn in fn_list if ReadFn(fn)['m']   > 1e-1]
		elif btype == 'G': fn_list = [fn for fn in fn_list if ReadFn(fn)['gap'] > 1e-1]
		elif btype == 'A': pass
		else:
			print('"%s" is wrong btype.' % btype)
			sys.exit()

		#Fb.write('%d\n' % (len(fn_list) * dfermi))
		for p in self.params1: Fp.write('%10s' % p)
		for p in self.params2: Fp.write('%22s' % p)
		Fp.write('\n')

		for b in ['%s%s%d' % (t, l, n) for t in ['e', 'w'] for l in self.labels for n in range(self.Nbs)]: Fb.write('%22s' % b)
		Fb.write('\n')

		band = []
		for fn in fn_list:
			fn_dict = ReadFn(fn)
			fn_dict['type'] = self.type_dict[fn_dict['type'][0]] + float(fn_dict['type'][1]) * 0.1

			with open(fn, 'r')                        as f: db = np.genfromtxt(f)
			with open(re.sub('band', 'ufw', fn), 'r') as f: du = np.genfromtxt(f)

			if btype == 'M':
				for p in self.params1: Fp.write('%10.1f'  % fn_dict[p])
				for p in self.params2: Fp.write('%22.16f' % fn_dict[p])
				Fp.write('\n')

				for b in np.hstack((np.ravel(db[self.points] - fn_dict['fermi']), np.ravel(du[self.points]))): Fb.write('%22.16f' % b)
				Fb.write('\n')

			elif btype == 'G':
				dntop = np.max(db[:, int(fn_dict['N'])-1]) - dfermi
				upbot = np.min(db[:, int(fn_dict['N'])])   + dfermi
				Nfermi = (upbot - dntop) // dfermi
				fermi_list = np.linspace(dntop, upbot, int(Nfermi))

				for fermi in fermi_list:
					fn_dict['fermi'] = fermi
					for p in self.params1: Fp.write('%10.1f'  % fn_dict[p])
					for p in self.params2: Fp.write('%22.16f' % fn_dict[p])
					Fp.write('\n')

					for b in np.hstack((np.ravel(db[self.points] - fermi), np.ravel(du[self.points]))): Fb.write('%22.16f' % b)
					Fb.write('\n')

		Fp.close()
		Fb.close()

		t1 = time.time()
		print('GenBand(%s, %s) : %fs' % (pn, bn, t1-t0))

	def GenDOSL(self, btype, eta=0.1): # Local DOS
		eta = float(eta)

		pn = '%s/params_%s.txt' % (self.path_save, btype)
		dn = '%s/dos_%sL_eta%.2f.h5' % (self.path_save, btype, eta)
		with open(pn, 'r') as f: p = np.genfromtxt(f, skip_header=1)
		with h5py.File('%s/energy.h5' % self.path_save, 'r') as f: e = f['energy'][()]

		t0 = time.time()		

		dos = np.zeros(self.BINS_MAX)
		for i in range(p.shape[0]):
			t = p[i, 0] * 10
			dir_name = '%s/JU%.2f_SOC%.2f' % (self.path_output, p[i, 1], 0)
			dos_name = 'dos_%s_eta%.2f_N%.1f_U%.1f' % ('%s%d' % (self.type_dict_r[str(int(t//10))], t%10), eta, p[i, 2], p[i, 3])
			fn = ['%s/%s' % (dir_name, fn) for fn in os.listdir(dir_name) if re.search(dos_name, fn)][0]

			with open(fn, 'r') as f: data = np.genfromtxt(f)
			data_e   = data[:, 0]
			data_dos = np.sum(data[:, 1:], axis=1)

			itp = interpolate.interp1d(data_e, data_dos)
			dos = np.vstack((dos, itp(e)))

		dos = np.delete(dos, 0, axis=0)
		with h5py.File(dn, 'w') as f: f.create_dataset('dos', data=dos, dtype='d')

		t1 = time.time()
		print('GenDOSL(%s) : %fs' % (dn, t1-t0))

	def ReadDn_(self, Dn):
		Dn_dict = {
			'type':  int(  re.sub('AF', '', re.search('AF\d',        Dn).group())),
			'JU':    float(re.sub('_J', '', re.search('_J\d+[.]\d+', Dn).group())),
			'N':     float(re.sub('Hz', '', re.search('Hz\d+[.]\d+', Dn).group())),
			'U':     float(re.sub('UF', '', re.search('UF\d+[.]\d+', Dn).group())),
			'm':     float(re.sub('_D', '', re.search('_D\d+[.]\d+', Dn).group())),
			'fermi': 0,
			'gap':   0,
		}
		
		return Dn_dict

	def GenDOSD(self, dtype, eta=0): # DMFT DOS
		eta = float(eta)
		dn = '%s/dos_D%s_eta%.2f.h5' % (self.path_save, dtype, eta)

		if   dtype == 'K': dos = np.zeros(self.BINS_MAX * len(self.labels))
		elif dtype == 'L': dos = np.zeros(self.BINS_MAX)
		else:
			print("'%s' is wrong dtype" % dtype)
			sys.exit()

		if eta > 0: dir_name = 'vdx'
		else:       dir_name = 'matsu/ome_final_result'
		pn = '%s/params_D%s_eta%.2f.txt' % (self.path_save, dtype, eta)

		with h5py.File('%s/energy.h5' % self.path_save, 'r') as f: e = f['energy'][()]

		t0 = time.time()		

		Fp = open(pn, 'w')
		for p in self.params1: Fp.write('%10s' % p)
		for p in self.params2: Fp.write('%22s' % p)
		Fp.write('\n')

		path_spec = 'dmft_old'
		d_list = ['%s/%s/lattice/%s' % (path_spec, d, dir_name) for d in os.listdir(path_spec)\
				if (re.search('oDir', d)\
				and float(re.sub('AF', '', re.search('AF\d', d).group())) > 0
				and float(re.sub('UF', '', re.search('UF\d+[.]\d+', d).group())) > 4
				and os.path.isfile('%s/%s/mkresult.bat' % (path_spec, d)))]

		"""
		#gn_list = ['%s/lattice/vdx/%s' % (d, f) for d in d_list for f in os.listdir(d+'/lattice/vdx') if re.search('kG.*ep%.2f' % eta, f)]
		#rn_list = [re.sub('lattice\S+', 'result', gn) for gn in gn_list]

		for i, rn in enumerate(rn_list):
			rn_dict = self.ReadRn_(rn)
			with open(rn + '/u%.3f/mag.dat'     % rn_dict['U'], 'r') as fr: rn_dict['m'] = np.genfromtxt(fr)[-1, 2] * 2
			with open(rn + '/u%.3f/filling.dat' % rn_dict['U'], 'r') as fr: rn_dict['N'] = np.genfromtxt(fr)[-1, 1] / 4
		"""

		if eta > 0:
			if dtype == 'K':
				for d in d_list:
					Dn_list = ['%s/%s' % (d, f) for f in os.listdir(d) if re.search('kG.*ep%.2f' % eta, f)]
					for Dn_h_list in [[re.sub('kG', 'k%s'%l, Dn) for l in self.labels] for Dn in Dn_list]:
						data_h_list = [np.genfromtxt(Dn_h) for Dn_h in Dn_h_list]

						Dn_dict = self.ReadDn_(Dn_h_list[0])
						for p in self.params1: Fp.write('%10.1f'  % Dn_dict[p])
						for p in self.params2: Fp.write('%22.16f' % Dn_dict[p])
						Fp.write('\n')

						itp_list = [interpolate.interp1d(data_h[:, 0], data_h[:, 1] * 6, fill_value='extrapolate') for data_h in data_h_list]
						dos = np.vstack((dos, np.ravel([itp(e) for itp in itp_list])))
		else:
			a_max = 3
			if dtype == 'K':
				for d in d_list:
					Dn_list = [['%s/%s' % (d, f) for f in os.listdir(d) if re.search('optimal_spectral_functions_k%s' % l, f)] for l in self.labels]
					for Dn_h_list in product(*[Dn_list[i] for i in range(self.Nhsp)]):
						data_h_list = [np.genfromtxt(Dn_h) for Dn_h in Dn_h_list]
						a_list = [[] for _ in Dn_h_list]

						for i, data_h in enumerate(data_h_list):
							for a in range(1, a_max+1):
								if data_h[np.min(np.where(data_h[:, 0] > 0)), a] < 0.5: a_list[i].append(a)

						if not sum([len(a) > 0 for a in a_list]) < 4:
							Dn_dict = self.ReadDn_(Dn_h_list[0])
							Dn_dict['N'] = int(''.join([re.sub('_tem', '', re.search('\d_tem', Dn_h).group()) for Dn_h in Dn_h_list]))

							for a_h_list in product(*[a for a in a_list]):
								for p in self.params1:      Fp.write('%10.1f'  % Dn_dict[p])
								for p in self.params2[:-1]: Fp.write('%22.16f' % Dn_dict[p])
								Fp.write('%22s\n' % ''.join([str(a_h) for a_h in a_h_list]))

								itp_list = [interpolate.interp1d(data_h[:, 0], data_h[:, a_h] * 6, fill_value='extrapolate') for data_h, a_h in zip(data_h_list, a_h_list)]
								dos = np.vstack((dos, np.ravel([itp(e) for itp in itp_list])))

			else:
				for d in d_list:
					Dn_list = ['%s/%s' % (d, f) for f in os.listdir(d) if re.search('optimal_spectral_functions_LDOS', f)]
					for Dn in Dn_list:
						data = np.genfromtxt(Dn)
						a_list = []

						for a in range(1, a_max+1):
							if data[np.min(np.where(data_h[:, 0] > 0)), a] < 0.5: a_list.append(a)

						if not len(a_list) < 1:
							Dn_dict = self.ReadDn_(Dn)

							for a in a_list:
								for p in self.params1:      Fp.write('%10.1f'  % Dn_dict[p])
								for p in self.params2[:-1]: Fp.write('%22.16f' % Dn_dict[p])
								Fp.write('%22s\n' % str(a))

								itp = interpolate.interp1d(data[:, 0], data[:, a] * 6)
								dos = np.vstack((dos, itp(e)))

		Fp.close()
		dos = np.delete(dos, 0, axis=0)
		with h5py.File(dn, 'w') as f: f.create_dataset('dos', data=dos, dtype='d')
			
		t1 = time.time()
		print('GenDOSD(%s, %s) : %fs' % (pn, dn, t1-t0))

	def Resample_(self, X, y, resamp):
		rsp = self.rsp_dict[resamp]

		X_rsp, y_rsp = rsp.fit_resample(X, y)
		X_rsp.index = X.index[rsp.sample_indices_]
		y_rsp.index = y.index[rsp.sample_indices_]

		return X_rsp, y_rsp

	def Preprocess_(self, dn, onehot, resamp, pca=0, multi_eta=None, is_reg=None):
		kind  = dn.split('_')[0]
		dtype = re.sub('%s_'%kind, '', re.search('%s_[a-zA-Z]+'%kind, dn).group())
		eta   = re.sub('eta', '', re.search('eta\d[.]\d+',  dn).group())
		if re.search('bins', dn): bins = int(re.sub('bins', '', re.search('bins\d+', dn).group()))
		else                    : bins = self.BINS_MAX

		# feats
		Nfts_dict = {
			'dos':  bins,
			'peak': self.PEAK_MAX,
		}
		if re.search('L', dn): feats = self.params + ['x%d' % i for i in range(Nfts_dict[kind])] 
		else:                  feats = self.params + ['%s%d' % (l, i) for l in self.labels for i in range(Nfts_dict[kind])]

		if   re.search('D', dtype): pn = 'params_%s_eta%s.txt' % (dtype, eta)
		else:                       pn = 'params_%s.txt' % dtype[0]
		with open(self.path_save+pn, 'r') as f: p = np.genfromtxt(f, skip_header=1)

		if multi_eta:
			with h5py.File(self.path_save+re.sub('eta\d[.]\d+', 'eta%.2f'%multi_eta.min(), dn), 'r') as f: d = f[kind][()]
			data = np.hstack((p, d))

			for eta in multi_eta:
				with h5py.File(self.path_save+re.sub('eta\d[.]\d+', 'eta%.2f'%eta, dn), 'r') as f: d = f[kind][()]
				data = np.vstack((data, np.hstack((p, d))))
		else:
			with h5py.File(self.path_save+dn, 'r') as f: d = f[kind][()]
			data = np.hstack((p, d))

		df = pd.DataFrame(data, columns=feats)
		df['type'] = df['type'].astype('int').astype('str').replace(self.type_dict_r)

		X = df.drop(self.params, axis=1)
		if is_reg: y = df[is_reg]
		else:      y = df['type']

		if onehot: y = pd.get_dummies(y)
		if resamp != 'none': X, y = self.Resample_(resamp, X, y)
		if pca: X = PCA(n_components=pca).fit_transform(X)

		return X, y, df

	def Predict_(self, X_test, y_test, df, mc, n_dict, onehot, verbose):
		t0 = time.time()

		y_pred  = mc.predict(X_test)
		y_score = mc.predict_proba(X_test)
		y_score = np.reshape(['%.2f' % sc for sc in np.ravel(y_score)], y_score.shape)

		if onehot:
			t_dict = {}
			for i, col in enumerate(y_test.columns): t_dict[str(i)] = col
			y_test = y_test.idxmax(axis=1)
			y_pred = np.array([t_dict[str(np.argmax(y))] for y in y_pred])

		acc = accuracy_score(y_test, y_pred)

		df_test = df.loc[y_test.index, self.params]
		df_test['type_p'] = y_pred
		df_test['score'] = list(map(str, y_score))
		
		t1 = time.time()
		if verbose: print('\n# Data : %s\n# Machine : %s\n# Resampler : %s\n# Accuracy : %f (%d/%d)\n# Elapsed Time : %fs\n'\
				% (n_dict['d'], n_dict['mc'], n_dict['rsp'], acc, len(y_pred)-len(df_test[df_test['type'] != df_test['type_p']]), len(y_pred), t1-t0))
		return y_test, y_pred, y_score, df_test

	def Tune(self, dn, mcn, resamp='none', verboes=True, pca=False, multi_eta=None):
		X, y, df = self.Preprocess_(dn, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp, pca=pca, multi_eta=multi_eta)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=1)
		
		cv = StratifiedKFold(random_state=1)
		tn = RandomizedSearchCV(mc, hp_dict[mcn[0]], cv, random_state=1)
		res = tn.fit(X_train, y_train)

		return res

	def Train(self, dn, mcn, resamp='none', verbose=True, pca=False, multi_eta=False, check_dataset=False, test_size=0.3):
		X, y, df = self.Preprocess_(dn, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp, pca=pca, multi_eta=multi_eta)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=1)

		if check_dataset: return X_train, X_test, y_train, y_test
		
		mc = self.mc_dict[mcn]['func']
		mc.fit(X_train, y_train)
		y_test, y_pred, y_score, df_test = self.Predict_(X_test, y_test, df, mc, {'d':dn, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'], verbose=verbose)

		return y_test, y_pred, y_score, df_test, mc

	def Validate(self, d1n, d2n, mcn, resamp='none', verbose=True, pca=False, multi_eta=False):
		onehot = self.mc_dict[mcn]['onehot']

		# train
		y_test1, _, _, df1_test, mc = self.Train(d1n, mcn, resamp=resamp, verbose=verbose, pca=pca, multi_eta=multi_eat)
		idcs = y_test1.index

		# validate
		X, y, df2 = self.Preprocess_(d2n, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp, pca=pca, multi_eta=multi_eta)

		if re.search('_D', d2n):
			X_test = X
			y_test = y
		else:
			X_test = X.loc[idcs, :]
			if onehot: y_test = y.loc[idcs]
			else:      y_test = y[idcs]

		y_test, y_pred, y_score, df2_test = self.Predict_(X_test, y_test, df2, mc, {'d':d2n, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'], verbose=verbose)

		return y_test, y_pred, y_score, df1_test, df2_test, mc

	def Regress(self, dn, rcn, resamp='none', verbose=True, multi_eta=False):
		X, y, df = self.Preprocess_(dn, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp, multi_eta=multi_eta, is_reg='m')
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1)

		rc = self.rc_dict[rcn]
		rc.fit(X_train, y_train)
		y_pred = rc.predict(X_test)
		y_score = rc.score(X_test, y_test)

		df_test = df.loc[y_test.index, self.params]
		df_test['m_p'] = y_pred

		return y_test, y_pred, y_score, df_test, rc

	def DrawDOS(self, dn, type, JU, N, U, fermi_idx=0, point=0, ax=0):
		dn = self.path_save + dn
		dtype = re.sub('dos_', '', re.search('dos_[A-Za-z]+', dn).group())
		eta = re.sub('eta', '', re.search('eta\d[.]\d+', dn).group())
		if re.search('bins', dn): bins = int(re.sub('bins', '', re.search('bins\d+', dn).group()))
		else                    : bins = self.BINS_MAX
		Nclus = self.BINS_MAX // bins

		JU = float(JU)
		N  = float(N)
		U  = float(U)
		fermi_idx = int(fermi_idx)
			
		if   re.search('DK', dtype): pn = '%s/params_DK_eta%s.txt' % (self.path_save, eta)
		elif re.search('DL', dtype): pn = '%s/params_DL_eta%s.txt' % (self.path_save, eta)
		else:                        pn = '%s/params_%s.txt' % (self.path_save, dtype[0])

		with open(pn, 'r') as f:
			for idx, line in enumerate(f):
				if re.search('%d[.]\d\s+%.1f\s+%.1f\s+%.1f\s+' % (self.type_dict[type], JU, N, U), line):
					idx += fermi_idx - 1
					break

		with h5py.File(dn, 'r') as f: d = f['dos'][()]
		with h5py.File('%s/energy.h5' % self.path_save, 'r') as f: e = f['energy'][()]

		x = d[idx][point*bins:(point+1)*bins]
		y = [np.average(e[i*Nclus:(i+1)*Nclus]) for i in range(bins)]

		if re.search('K', dtype):
			with h5py.File(re.sub('dos', 'peak', dn), 'r') as f: p = f['peak'][()]
			peak0 = p[idx][point*self.PEAK_MAX:(point+1)*self.PEAK_MAX]
			peak_idx = [np.where(abs(y - peak0[i]) < 1e-6)[0] for i in range(self.PEAK_MAX)]
			peak = [(x[i[0]], i[0]) for i in peak_idx]
		else: peak = 0

		#if not ax: fig, ax = plt.subplots(dpi=600, constrained_layout=True)
		if ax: ax.plot(x, y, label='%s %s_eta%s/JU%.1f_N%.1f_U%.1f_fermi%d' % (type, dtype, eta, JU, N, U, fermi_idx))

		return x, y, peak
	
	def DrawSpec(self, dn, type, JU, N, U, ax=0, show_yticks=True):
		JU  = float(JU)
		N   = float(N)
		U   = float(U)

		dtype = re.sub('dos_', '', re.search('dos_[A-Za-z]+', dn).group())
		eta   = float(re.sub('eta', '', re.search('eta\d+[.]\d+', dn).group()))

		path_band = '%s/JU%.2f_SOC%.2f' % (self.path_output, JU, 0)
		os.makedirs(re.sub('output', 'diagram', path_band), exist_ok=True)

		bn = ['%s/%s' % (path_band, f) for f in os.listdir(path_band) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		bn_dict = ReadFn(bn)

		with open(bn, 'r')                        as f: db = np.genfromtxt(f) - bn_dict['dntop']
		with open(re.sub('band', 'ufw', bn), 'r') as f: du = np.genfromtxt(f)
		with h5py.File('%s/energy.h5' % self.path_save, 'r') as f: energy= f['energy'][()]

		# options
		if re.search('f', dtype): weight = [1 if e < 0 else -1 for e in energy]
		else:                     weight = [1 for _ in energy]

		if re.search('b', dtype): broad = [0.1 + abs(e / np.min(energy)) * eta for e in energy]
		else:                     broad = [eta for _ in energy]

		spec = []	
		for i in range(db.shape[0]):
			for e, w, b in zip(energy, weight, broad):
				sum = 0
				for j in range(db.shape[1]):
					sum += (b / ((e - db[i][j])**2 + b**2)) * du[i][j] * w
				spec.append([i, e, sum / np.pi])
		spec = np.array(spec)

		X = np.reshape(spec[:, 0], (db.shape[0], len(energy)))
		Y = np.reshape(spec[:, 1], (db.shape[0], len(energy)))
		Z = np.reshape(spec[:, 2], (db.shape[0], len(energy)))

		if not ax: fig, ax = plt.subplots(dpi=600)

		cmax = 11
		cf = ax.contourf(X, Y, Z, levels=np.linspace(0, cmax, 100), cmap='jet')
		#cf = ax.contourf(X, Y, Z, levels=100, cmap='jet')

		#cb = plt.colorbar(ct)
		#cb.set_ticks([0, cmax])
		#cb.set_label(r'$\rho_{\mathbf{k}}$', fontdict={'fontsize':'medium'}, rotation=0)

		if show_yticks:
			ax.set_ylabel(r'$E - E_F$')
		else:
			ax.set_yticklabels([])
			ax.set_ylabel('')

		ax.set_xticks(self.points_g)
		ax.set_xticklabels(self.labels_g)
		ax.set_ylim(np.min(db)-1, 0)
		#ax.set_title(r'%s-type $J/U = %.1f$ $N = %.1f$ $U = %.1f$' % (type[0], JU, N, U), loc='left')

		#fig.tight_layout()
		#fig.savefig('%s' % re.sub('output', 'diagram', re.sub('band', 'spec_%s' % dtype, re.sub('txt', 'pdf', bn))))
		#plt.show()

		return cf

	def annotate_heatmap_(self, im, data_txt=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None, **textkw):
		data = im.get_array()

		# Normalize the threshold to the images color range.
		if threshold is not None:
			threshold = im.norm(threshold)
		else:
			threshold = im.norm(data.max())/2.

		# Set default alignment to center, but allow it to be
		# overwritten by textkw.
		kw = dict(horizontalalignment="center",
		verticalalignment="center")
		kw.update(textkw)

		# Get the formatter in case a string is supplied
		if isinstance(valfmt, str):
			valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

		# Loop over the data and create a `Text` for each "pixel".
		# Change the text's color depending on the data.
		texts = []
		for i in range(data.shape[0]):
			for j in range(data.shape[1]):
				kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
				text = im.axes.text(j, i, valfmt(data_txt[i, j], None), fontsize=32, **kw)
				texts.append(text)

		return texts

	def DrawConfMat(self, y_test, y_pred, ax):
		labels = ['a', 'c', 'g']
		mat = confusion_matrix(y_test, y_pred, labels=labels)
		acc = [mat[i, :]/mat[i, :].sum() for i in range(3)]
		norm = plt.Normalize(0, 1)

		im = ax.imshow(acc, cmap='Blues', norm=norm)
		self.annotate_heatmap_(im, data_txt=mat, valfmt='{x:d}')

		#ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', cmap='Blues', values_format='.3f', ax=ax, text_kw={'size':23}, colorbar=False)
		#ConfusionMatrixDisplay.from_predictions(y_test, y_pred, cmap='Blues', ax=ax, text_kw={'size':30}, colorbar=False)
		#xtlabels = [item.get_text().upper() for item in ax.get_xticklabels()]
		#ytlabels = [item.get_text().upper() for item in ax.get_yticklabels()]

		return im

	def DrawROC(self, y_test, y_pred, ax=0):
		y_testb = label_binarize(y_test, classes=types)
		n_classes = y_testb.shape[1]

		fpr = dict()
		tpr = dict()
		roc_auc = dict()

		if not ax: fig, ax = plt.subplots(dpi=600)

		for i in range(n_classes):
			fpr[i], tpr[i], _ = roc_curve(y_testb[:, i], y_score[:, i])
			roc_auc[i] = auc(fpr[i], tpr[i])
			ax.plot(fpr[i], tpr[i], label='%s (%.3f)' % (types[i], roc_auc[i]))

		ax.plot([0, 1], [0, 1], 'k--')
		ax.set_xlim([-0.05, 1.0])
		ax.set_ylim([0.0, 1.05])
		ax.set_xlabel('False Positive Rate')
		ax.set_ylabel('True Positive Rate')
		ax.legend(loc="lower right")
		plt.show()

	def DrawTree(self, mc, feats, ax=0):
		if not ax: fig, ax = plt.subplots(dpi=600)
		#tree.plot_tree(mc.estimators_[0], feature_names=feats, max_depth=4, filled = True, ax=ax, fontsize=15)
		tree.plot_tree(mc.estimators_[0], feature_names=feats, filled = True, ax=ax, fontsize=15)
		plt.show()

	def DrawImps(self, mc, feats, ax=0):
		importances = mc.feature_importances_
		std = np.std([tree.feature_importances_ for tree in mc.estimators_], axis=0)
		forest_importances = pd.Series(importances, index=feats)

		if not ax: fig, ax = plt.subplots(dpi=600)
		forest_importances.plot.bar(ax=ax)
		plt.show()

		"""	
		fig, ax = plt.subplots(self.Nhsp, 1, sharey=True, figsize=figsize, dpi=600)
		for i in range(self.Nhsp):
			imp_s = forest_importances[bins*i:bins*(i+1)]
			std_s = std[bins*i:bins*(i+1)]
			#imp_s.plot.bar(yerr=std_s, ax=ax[i])
			imp_s.plot.bar(ax=ax[i])

			ax[i].set_title("%s" % self.labels[i])
			#ax.set_title("Feature importances using MDI")
			#ax.set_ylabel("Mean decrease in impurity")
		fig.tight_layout()
		"""
		plt.show()

