# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import os
import re
import h5py
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier	

from imblearn.under_sampling import RandomUnderSampler
from imblearn.under_sampling import NearMiss
from imblearn.under_sampling import EditedNearestNeighbours
from imblearn.under_sampling import CondensedNearestNeighbour

from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc

from .read import ReadInfo, ReadFn, GenGroundIdx

class MagStrSelf:
	def __init__(self, name, num_thread=1):
		self.path_input, self.path_output, self.info_path, self.info_cell = ReadInfo(name);
		self.path_save = self.path_output + '/magstrself/'
		os.makedirs(self.path_save, exist_ok=True)

		self.Ni = self.info_cell['a'][0]
		self.Nc = self.info_cell['a'][1]
		self.Ns = self.Ni * self.Nc
		self.Nb = self.Ns * 2

		self.type_dict   = {'a':1, 'c':2, 'g':3}
		self.type_dict_r = {'1':'a', '2':'c', '3':'g'}

		self.params = ['type', 'JU', 'N', 'U', 'm']

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
		}

		self.figsize=(11, 6)

	def ReadSn_(self, sn):
		sn_dict = {
			'type': re.sub('sol_', '', re.search('sol_[a-z]\d?', sn).group()),
			'JU': float(re.sub('JU', '', re.search('JU\d+[.]\d+', sn).group())),
			'N': float(re.sub('N', '', re.search('N\d+[.]\d+', sn).group())),
			'U': float(re.sub('U', '', re.search('U\d+[.]\d+', sn).group())),
			'm': 0,
		}
		
		return sn_dict

	def GenSelfEnergy(self):
		pn = '%s/params.txt'     % (self.path_save)
		sn = '%s/selfenergy.txt' % (self.path_save)
		fp = open(pn, 'w')
		fs = open(sn, 'w')

		t0 = time.time()

		bn_list = ['%s/%s/%s' % (self.path_output, d, f)\
				for d in os.listdir(self.path_output)     if re.search('JU',   d)\
				for f in os.listdir(self.path_output + d) if re.search('band', f)]
		bn_list = [bn_list[i] for i in GenGroundIdx(bn_list, exclude_f=True)] 
		#bn_list = [bn for bn in bn_list if ReadFn(fn)['m'] > 0.1]
		sn_list = [re.sub('band', 'sol', re.sub(re.search('_n\S+_', bn).group(), '_', bn)) for bn in bn_list]

		fp.write('%d\n' % len(sn_list))
		fs.write('%d\n' % len(sn_list))

		for p in self.params: fp.write('%22s' % p)
		fp.write('\n')

		for s in ['%d_%d%d' % (i, j, k) for i in range(self.Ni) for j in range(self.Ns) for k in range(self.Ns)]: fs.write('%22s' % s)
		fs.write('\n')

		for sn in sn_list:
			sn_dict = self.ReadSn_(sn)
			sn_dict['type'] = self.type_dict[sn_dict['type'][0]]
			with open(sn, 'r') as f: oc = np.genfromtxt(f)[-1, 1:]
			for i in range(self.Ns): sn_dict['m'] += (oc[i] - oc[i + self.Ns]) * -1**(i / self.Nc)
			for p in self.params: fp.write('%22.16f' % sn_dict[p])
			fp.write('\n')

		fp.close()
		fs.close()

		t1 = time.time()
		print('GenSelfEnergy(%s, %s) : %fs' % (pn, sn, t1-t0))

	def ReadRn_(self, rn):
		rn_dict = {
			'type': int(re.sub('AF', '', re.search('AF\d', rn).group())),
			'JU': float(re.sub('Hz', '', re.search('Hz\d+[.]\d+', rn).group())),
			'U': float(re.sub('UF', '', re.search('UF\d+[.]\d+', rn).group())),
			'n': 0,
			'm': 0,
		}
		
		return rn_dict

	def GenSelfEnergyDMFT(self):
		pn = '%s/params_d.txt' % (self.path_output)
		dn = '%s/selfenergy_d.txt' % (self.path_output, eta)
		fp = open(pn, 'w')
		fs = open(sn, 'w')

		t0 = time.time()		

		path_spec = '/home/Shared/BaOsO3/dmft_spec'
		d_list = ['%s/%s' % (path_spec, d) for d in os.listdir(path_spec)\
				if (re.match('oDir', d)\
				and float(re.sub('AF', '', re.search('AF\d', d).group())) > 0)] # choose only AF types (AF0 is FM)
		gn_list = ['%s/lattice/vdx/%s' % (d, f) for d in d_list for f in os.listdir(d+'/lattice/vdx') if re.search('kG.*ep%.2f' % eta, f)]
		rn_list = ['%s/result/%s' % (d, os.listdir(d+'/result')[0]) for d in d_list]

		fp.write('%d\n' % len(rn_list))
		fs.write('%d\n' % len(rn_list))

		for p in self.params: fp.write('%22s' % p)
		f.write('\n')

		for rn in rn_list:
			rn_dict = self.ReadRn_(rn)
			with open(rn + '/mag.dat',     'r') as fr: rn_dict['m'] = np.genfromtxt(fr)[-1, 2] * 2
			with open(rn + '/filling.dat', 'r') as fr: rn_dict['n'] = np.genfromtxt(fr)[-1, 1]
			for p in self.params: fp.write('%22.16f'  % rn_dict[p])
			fp.write('\n')

		fp.close()
		fs.close()

		t1 = time.time()
		print('GenSelfEnergyDMFT(%s) : %fs' % (dn, t1-t0))

	def Resample_(self, X, y, resamp):
		rsp = self.rsp_dict[resamp]

		X_rsp, y_rsp = rsp.fit_resample(X, y)
		X_rsp.index = X.index[rsp.sample_indices_]
		y_rsp.index = y.index[rsp.sample_indices_]

		return X_rsp, y_rsp

	def Preprocess_(self, dn, onehot=False, resamp='none'):
		dn = self.path_output + dn
		if re.search('bins', dn): bins = int(re.sub('bins', '', re.search('bins\d+', dn).group()))
		else                    : bins = self.BINS_MAX
		dtype = int(re.sub('dt', '', re.search('dt\d', dn).group()))
			
		if dtype: pn = '%s/params.txt'  % self.path_output
		else:     pn = '%s/params0.txt' % self.path_output

		with open(pn, 'r') as f: p = np.genfromtxt(f, skip_header=2)
		with h5py.File(dn, 'r') as f: d = f['dos'][()]

		data = np.hstack((p, d))
		df = pd.DataFrame(data, columns=self.params + ['%s%d' % (l, b) for l in self.labels for b in range(bins)])
		df['type'] = df['type'].astype('int').astype('str').replace(self.type_dict_r)

		X = df.drop(self.params, axis=1)
		y = df['type']

		if onehot: y = pd.get_dummies(y)
		if resamp != 'none': X, y = self.Resample_(resamp, X, y)
		#idx_rsp = y.index.to_list()

		return X, y, df

	def Predict_(self, X_test, y_test, df, mc, n_dict, onehot=False, verbose=True):
		t0 = time.time()

		y_pred  = mc.predict(X_test)
		y_score = mc.predict_proba(X_test)

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
		if verbose: print('\n# Data : %s\n# Machine : %s\n# Resampler : %s\n# Max DOS : %f\n# Accuracy : %f (%d/%d)\n# Elapsed Time : %fs\n'\
				% (n_dict['d'], n_dict['mc'], n_dict['rsp'], np.max(df.drop(self.params, axis=1).values), acc, len(y_pred)-len(df_test[df_test['type'] != df_test['type_p']]), len(y_pred), t1-t0))
		return y_test, y_pred, y_score, df_test

	def Train(self, dn, mcn, resamp='none', tune=False, verbose=True):
		t0 = time.time()

		X, y, df = self.Preprocess_(dn, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=1)
		
		# tune
		if tune:
			cv = StratifiedKFold(random_state=1)
			tn = RandomizedSearchCV(mc, hp_dict[mcn[0]], cv, random_state=1)
			res = tn.fit(X_train, y_train)

			t1 = time.time()
			return res

		# train
		else:
			mc = self.mc_dict[mcn]['func']
			mc.fit(X_train, y_train)
			y_test, y_pred, y_score, df_test = self.Predict_(X_test, y_test, df, mc, {'d':dn, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'], verbose=verbose)

			t1 = time.time()
			return y_test, y_pred, y_score, df_test, mc

	def Validate(self, d1n, d2n, mcn, resamp='none', verbose=True):
		t0 = time.time()

		# train
		_, _, _, df1_test, mc = self.Train(d1n, mcn, resamp=resamp, tune=False, verbose=verbose)

		# validate
		X_test, y_test, df2 = self.Preprocess_(d2n, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp)
		y_test, y_pred, y_score, df2_test = self.Predict_(X_test, y_test, df2, mc, {'d':d2n, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'], verbose=verbose)

		t1 = time.time()
		return y_test, y_pred, y_score, df1_test, df2_test, mc

	def DrawDOS(self, dn, type, JU, N, U, ax, xmin=0, xmax=1024, alpha=1):
		dn = self.path_output + dn
		type = int(type)
		JU   = float(JU)
		N    = float(N)
		U    = float(U)

		dtype = int(re.sub('dt', '', re.search('dt\d', dn).group()))
		eta   = float(re.sub('eta', '', re.search('eta\d+[.]\d+', dn).group()))

		if dtype: pn = '%s/params.txt'  % self.path_output
		else:     pn = '%s/params0.txt' % self.path_output

		with open(pn, 'r') as f:
			for i, line in enumerate(f):
				if re.search('%.1f\s+%.1f\s+%.1f\s+%.1f\s+' % (type, JU, N, U), line):
					i -= 2
					break

		with h5py.File(dn, 'r') as f: d = f['dos'][()]
		x = range(xmin, xmax)
		y = d[i]

		ax.plot(x, y[xmin:xmax], label='%d dt%d_eta%.4f/JU%.1f_N%.1f_U%.1f' % (type, dtype, eta, JU, N, U), alpha=alpha)

	def DrawSpec(self, dn, type, JU, N, U):
		JU  = float(JU)
		N   = float(N)
		U   = float(U)

		dtype = re.sub('dt', '', re.search('dt\d', dn).group())
		eta   = float(re.sub('eta', '', re.search('eta\d+[.]\d+', dn).group()))

		path_band = '%s/JU%.2f_SOC%.2f' % (re.sub('/magstrself', '', self.path_output), JU, 0)
		os.makedirs(re.sub('output', 'diagram', path_band), exist_ok=True)

		bn = ['%s/%s' % (path_band, f) for f in os.listdir(path_band) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		bn_dict = ReadFn(bn)

		with open(bn, 'r')                        as f: db = np.genfromtxt(f) - bn_dict['dntop']
		with open(re.sub('band', 'ufw', bn), 'r') as f: du = np.genfromtxt(f)
		e_range = np.linspace(db.min(), db.max(), 128)

		# options
		if re.search('f', dtype): weight = [1 if e < 0 else 0 for e in e_range]
		else:                     weight = [1 for _ in e_range]

		if re.search('b', dtype): broad = [eta + abs(e) / 10 for e in e_range]
		else:                     broad = [eta for _ in e_range]

		spec = []	
		for i in range(db.shape[0]):
			for e, w, b in zip(e_range, weight, broad):
				sum = 0
				for j in range(db.shape[1]):
					sum += (b / ((e - db[i][j])**2 + b**2)) * du[i][j] * w
				spec.append([i, e, sum / np.pi])
		spec = np.array(spec)

		X = np.reshape(spec[:, 0], (db.shape[0], len(e_range)))
		Y = np.reshape(spec[:, 1], (db.shape[0], len(e_range)))
		Z = np.reshape(spec[:, 2], (db.shape[0], len(e_range)))

		fig, ax = plt.subplots(figsize=self.figsize)
		ct = ax.contourf(X, Y, Z, levels=1000, cmap='jet')
		cb = fig.colorbar(ct)
		#cb.set_ticks(np.arange(0, 1.1, 0.2))
		ax.set_xticks(self.points)
		ax.set_xticklabels(self.labels)
		ax.set_ylabel(r'$E - E_F$')
		ax.set_title('%s-type J/U = %.1f N = %.1f U = %.1f eta = %.4f' % (type[0], JU, N, U, eta), loc='left')

		fig.tight_layout()
		fig.savefig('%s' % re.sub('output', 'diagram', re.sub('band', 'spec_dt%s' % dtype, re.sub('txt', 'png', bn))))
		plt.show()

	def DrawConfMat(self, y_test, y_pred, title=''):
		fig, ax = plt.subplots(figsize=self.figsize)
		ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', cmap='Blues', values_format='.2f', ax=ax)
		ax.set_title(title, loc='left')
		plt.show()

	def DrawROC(self, y_test, y_pred, title=''):
		y_testb = label_binarize(y_test, classes=types)
		n_classes = y_testb.shape[1]

		fpr = dict()
		tpr = dict()
		roc_auc = dict()

		fig, ax = plt.subplots(figsize=figsize)

		for i in range(n_classes):
			fpr[i], tpr[i], _ = roc_curve(y_testb[:, i], y_score[:, i])
			roc_auc[i] = auc(fpr[i], tpr[i])

			ax.plot(fpr[i], tpr[i], label='%s (%.3f)' % (types[i], roc_auc[i]))

		ax.plot([0, 1], [0, 1], 'k--')
		ax.set_xlim([-0.05, 1.0])
		ax.set_ylim([0.0, 1.05])
		ax.set_xlabel('False Positive Rate')
		ax.set_ylabel('True Positive Rate')
		ax.set_title(title, loc='left')
		ax.legend(loc="lower right")
		plt.show()
