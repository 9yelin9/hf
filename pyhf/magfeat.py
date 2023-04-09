# pyhf3/magfeat.py : select important features to predict magnetic structure

import os
import re
import sys
import h5py
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score, confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV
from sklearn.preprocessing import label_binarize
from sklearn import tree

# classifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier	
from sklearn.svm import SVC

# resampler
from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NearMiss
from imblearn.under_sampling import EditedNearestNeighbours
from imblearn.under_sampling import CondensedNearestNeighbour

# feature selector
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import GenericUnivariateSelect, chi2
from sklearn.decomposition import PCA

class MagFeat:
	def __init__(self, name='baoso3_dU1.0', num_thread=1):
		self.bins_max = 1024
		self.Nb = 12
		self.emin = -8
		self.emax = 8
		self.Nkb = 1024

		self.Nhsp = 4
		self.hsp_point   = [0, 198, 396, 1023]
		self.hsp_label   = [r'$\Gamma$', 'X', 'M', 'R']
		self.hsp_point_g = [0, 198, 396, 677, 1023]
		self.hsp_label_g = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']

		self.type_dict = {'a':1, 'c':2, 'g':3}
		self.type_dict_r = {'1':'a', '2':'c', '3':'g'}
		self.params = ['type', 'JU', 'N', 'U', 'm', 'gap']

		self.path_output = 'output/%s/' % name
		self.path_save = 'output/%s/magfeat/' % name
		os.makedirs(self.path_save, exist_ok=True)

		self.random_state = 1
		self.mc_dict = {
			'rf':  RandomForestClassifier(random_state=self.random_state, n_jobs=num_thread),
			'xgb': XGBClassifier(random_state=self.random_state, nthread=num_thread),
			'lgb': LGBMClassifier(random_state=self.random_state, n_jobs=num_thread),
			'cat': CatBoostClassifier(random_state=self.random_state, thread_count=num_thread, silent=True),
			'lr':  LogisticRegression(random_state=self.random_state, n_jobs=num_thread, multi_class='multinomial', max_iter=10000),
			'svm': SVC(random_state=self.random_state, probability=True),
		}
		self.rsp_dict = {
			'rus': RandomUnderSampler(random_state=self.random_state),
			'nm':  NearMiss(version=1),
			'enn': EditedNearestNeighbours(),
			'cnn': CondensedNearestNeighbour(random_state=self.random_state),
		}
		self.ft_dict = {
			'pca': PCA(n_components=10),
		}

	def GenEnergy(self):
		t0 = time.time()

		en = '%s/energy.h5' % self.path_save
		e = np.linspace(self.emin, self.emax, self.bins_max)
		with h5py.File(en, 'w') as f: f.create_dataset('energy', data=e, dtype='d')

		t1 = time.time()
		print('GenEnergy(%s) : %fs' % (en, t1-t0))
	
	def FnDict_(self, fn):
		fn_dict = {
			'type':  self.type_dict[re.sub('band_', '', re.search('band_[a-z]\d', fn).group())[0]],
			'JU':    float(re.sub('JU',     '', re.search('JU\d+[.]\d+',         fn).group())),
			'N':     float(re.sub('_N',     '', re.search('_N\d+[.]\d+',         fn).group())),
			'U':     float(re.sub('_U',     '', re.search('_U\d+[.]\d+',         fn).group())),
			'n':     float(re.sub('_n',     '', re.search('_n\d+[.]\d+',         fn).group())),
			'm':     float(re.sub('_m',     '', re.search('_m[-]?\d+[.]\d+',     fn).group())),
			'e':     float(re.sub('_e',     '', re.search('_e[-]?\d+[.]\d+',     fn).group())),
			'fermi': float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fn).group())),
			'gap':   float(re.sub('_gap',   '', re.search('_gap[-]?\d+[.]\d+',   fn).group())),
		}

		return fn_dict

	def FnDictD_(self, fn):
		fn_dict = {
			'type':   int(re.sub('AF', '', re.search('AF\d',        fn).group())),
			'JU':   float(re.sub('_J', '', re.search('_J\d+[.]\d+', fn).group())),
			'N':    float(re.sub('Hz', '', re.search('Hz\d+[.]\d+', fn).group())),
			'U':    float(re.sub('UF', '', re.search('UF\d+[.]\d+', fn).group())),
			'm':    float(re.sub('_D', '', re.search('_D\d+[.]\d+', fn).group())),
			'gap':  float(re.sub('th', '', re.search('\d+th',       fn).group())),
		}
		
		return fn_dict

	def GroundOnly_(self, fn_list):
		params = ['type', 'JU', 'N', 'U', 'e']

		data = np.zeros(len(params))
		for fn in fn_list:
			fn_dict = self.FnDict_(fn)
			data = np.vstack((data, [fn_dict[p] for p in params]))
		data = np.delete(data, 0, axis=0)

		df = pd.DataFrame(data, columns=params)
		df = df.sort_values(by=['JU', 'N', 'U', 'e'])
		df = df.drop_duplicates(subset=['JU', 'N', 'U', 'type'], keep='first')
		grd_list = df.index.to_list()

		return grd_list

	def GenList(self, dtype):
		t0 = time.time()

		pn = '%s/params_%s.txt' % (self.path_save, dtype)
		ln = '%s/list_%s.txt'   % (self.path_save, dtype)

		if re.search('D', dtype):
			path_dmft = 'dmft_old'
			vdx_list = ['%s/%s/lattice/vdx' % (path_dmft, d) for d in os.listdir(path_dmft)\
					if (re.search('oDir', d)\
					and os.path.isfile('%s/%s/mkresult.bat' % (path_dmft, d))
					and float(re.sub('AF', '', re.search('AF\d', d).group())) > 0
					and float(re.sub('UF', '', re.search('UF\d+[.]\d+', d).group())) > 4)]

			fn_list = ['%s/%s' % (d, f) for d in vdx_list for f in os.listdir(d) if (re.match('u', f) and re.search('ep0.02', f))]

			with open(pn, 'w') as f:
				for p in self.params: f.write('%16s' % p)
				f.write('\n')
				for fn in fn_list:
					for p in self.params:
						f.write('%16.10f' % self.FnDictD_(fn)[p])
					f.write('\n')

			with open(ln, 'w') as f:
				f.write('%d\n' % len(fn_list))
				f.write('\n'.join(['%s' % fn for fn in fn_list]))
		else:
			fn_list = ['%s/%s/%s' % (self.path_output, d, f)\
					for d in os.listdir(self.path_output)   if re.search('JU',    d)\
					for f in os.listdir(self.path_output+d) if re.search('band_', f)]
			fn_list = [fn_list[i] for i in self.GroundOnly_(fn_list)]

			if   re.search('M', dtype): fn_list = [fn for fn in fn_list if self.FnDict_(fn)['m']   > 1e-1]
			elif re.search('G', dtype): fn_list = [fn for fn in fn_list if self.FnDict_(fn)['gap'] > 1e-1]
			elif re.search('A', dtype): pass
			else:
				print('"%s" is wrong dtype.' % dtype)
				sys.exit()

			with open(pn, 'w') as f:
				for p in self.params: f.write('%16s' % p)
				f.write('\n')
				for fn in fn_list:
					for p in self.params:
						f.write('%16.10f' % self.FnDict_(fn)[p])
					f.write('\n')

			with open(ln, 'w') as f:
				f.write('%d\n' % len(fn_list))
				if   re.search('K', dtype): f.write('\n'.join(['%s\t%.16f' % (fn, self.FnDict_(fn)['fermi']) for fn in fn_list]))
				elif re.search('L', dtype): f.write('\n'.join(['%s' % re.sub('band', 'sol', re.sub('_n.+_', '_', fn)) for fn in fn_list]))

		t1 = time.time()
		print('GenList(%s, %s) : %fs' % (pn, ln, t1-t0))

	def Resample_(self, X, y, rspn):
		rsp = self.rsp_dict[rspn]

		X_rsp, y_rsp = rsp.fit_resample(X, y)
		X_rsp.index = X.index[rsp.sample_indices_]
		y_rsp.index = y.index[rsp.sample_indices_]

		return X_rsp, y_rsp

	def Preprocess_(self, dn, mcn, rspn, ftn, is_hsp=0, is_etas=0):
		kind  = dn.split('_')[0]
		dtype = re.sub('%s_' % kind, '', re.search('%s_[a-zA-Z]+' % kind, dn).group())
		eta   = re.sub('eta', '', re.search('eta\d[.]\d+', dn).group())
		bins  = int(re.sub('bins', '', re.search('bins\d+', dn).group()))

		# feature
		Nft_dict = {
			'dos':  bins,
			'peak': self.PEAK_MAX,
		}
		if re.search('L', dn): feats = self.params + ['x%d' % i for i in range(Nft_dict[kind])] 
		else:                  feats = self.params + ['x%d_%d' % (j, i) for j in range(Nft) for i in range(Nft_dict[kind])]

		if re.search('D', dtype): pn = 'params_%s_eta%s.txt' % (dtype, eta)
		else:                     pn = 'params_%s.txt' % dtype[0]
		with open(self.path_save+pn, 'r') as f: p = np.genfromtxt(f, skip_header=1)

		if is_etas:
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
		if is_hsp:         X = X.loc[:, ['x%d_%d' % (j, i) for j in self.hsp_point for i in range(Nft_dict[kind])]]
		if mcn  == 'xgb':  y = pd.get_dummies(y)
		if rspn != 'none': X, y = self.Resample_(X, y, rspn)
		if ftn  != 'none': X = self.ft_dict[ftn].fit_transform(X)

		return X, y, df

	def Predict_(self, X_test, y_test, df, dn, mcn, rspn, ftn, is_verbose=1):
		t0 = time.time()

		y_pred  = mc.predict(X_test)
		y_score = mc.predict_proba(X_test)
		y_score = np.reshape(['%.2f' % sc for sc in np.ravel(y_score)], y_score.shape)

		if mcn == 'xgb':
			t_dict = {}
			for i, col in enumerate(y_test.columns): t_dict[str(i)] = col
			y_test = y_test.idxmax(axis=1)
			y_pred = np.array([t_dict[str(np.argmax(y))] for y in y_pred])

		acc = accuracy_score(y_test, y_pred)

		df_test = df.loc[y_test.index, self.params]
		df_test['type_p'] = y_pred
		df_test['score'] = list(map(str, y_score))
		
		t1 = time.time()
		if is_verbose:
			print('\n# Data : %s\n# Machine : %s\n# Resampler : %s\n# Feature selector : %s\n# Accuracy : %f (%d/%d)\n# Elapsed Time : %fs\n'\
					% (dn, mcn, rspn, ftn, acc, len(y_pred)-len(df_test[df_test['type'] != df_test['type_p']]), len(y_pred), t1-t0))
		return y_test, y_pred, y_score, df_test

	def Train(self, dn, mcn, rspn, ftn, tune_dict=0, test_size=0.3, is_hsp=0, is_etas=0, is_checkset=0, is_verbose=1):
		X, y, df = self.Preprocess_(dn, mcn, rspn, ftn, is_hsp=is_hsp, is_etas=is_etas)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=self.random_state)

		if tune_dict:
			cv = StratifiedKFold(random_state=self.random_state)
			tn = RandomizedSearchCV(self.mc_dict[mcn], tune_dict, cv, random_state=self.random_state)
			res = tn.fit(X_train, y_train)
			return res
		if is_checkset:
			return X_train, X_test, y_train, y_test
		
		mc = self.mc_dict[mcn]
		mc.fit(X_train, y_train)
		y_test, y_pred, y_score, df_test = self.Predict_(X_test, y_test, df, dn, mcn, rspn, ftn, is_verbose=is_verbose)

		return y_test, y_pred, y_score, df_test, mc

	def Validate(self, dn, mcn, rspn, ftn, test_size=0.3, is_hsp=0, is_etas=0, is_verbose=1):
		if not type(dn) is list:
			print('Type of dn is not a list. (dn[0] : trainset name, dn[1] : testset name)')
			sys.exit()

		y_test0, _, _, df0_test, mc = self.Train(dn[0], mcn, trspn, ftn, test_size=test_size, is_hsp=is_hsp, is_etas=is_etas, is_verbose=is_verbose)
		idx = y_test0.index

		X, y, df1 = self.Preprocess_(dn[1], onehot=self.mc_dict[mcn]['onehot'], resamp=resamp, pca=pca, multi_eta=multi_eta)

		if re.search('_D', dn[1]):
			X_test = X
			y_test = y
		else:
			X_test = X.loc[idx, :]
			y_test = y.loc[idx] if mcn == 'xgb' else y[idx]

		y_test, y_pred, y_score, df1_test = self.Predict_(X_test, y_test, df1, dn, mcn, rspn, ftn, is_verbose=is_verbose)

		return y_test, y_pred, y_score, [df0_test, df1_test], mc

