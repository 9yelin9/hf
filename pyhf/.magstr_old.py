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

class MagStr:
	def __init__(self, name, num_thread=1):
		self.path_input, self.path_output, self.info_path, _ = ReadInfo(name);
		os.makedirs('%s/magstr/' % self.path_output, exist_ok=True)
		self.path_output += '/magstr/'

		info_path_nodup = [] # drop duplicated hsp
		for point, label in self.info_path:
			if label not in [label for point, label in info_path_nodup]:
				info_path_nodup.append([point, label])
		self.points = [point for point, _ in info_path_nodup]
		self.labels = [label for _, label in info_path_nodup]

		#self.info_cell = info_cell
		#self.Ni = self.info_cell['a'][0]
		#self.Nc = self.info_cell['a'][1]
		#self.Nb = self.Ni * self.Nc * 2

		self.Nbs = 12
		self.Nhsp = 4
		self.BINS_MAX = 1024
		self.tol_U = 0.1

		self.type_dict   = {'a':1, 'c':2, 'g':3}
		self.type_dict_r = {'1':'a', '2':'c', '3':'g'}

		self.params1 = ['type', 'JU', 'N', 'U']
		self.params2 = ['m', 'dntop', 'gap']
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
		}

		self.figsize=(11, 6)

	def GenBand(self):
		pn = '%s/params.txt' % (self.path_output)
		bn = '%s/band.txt' % (self.path_output)
		Fp = open(pn, 'w')
		Fb = open(bn, 'w')

		t0 = time.time()

		fn_list = ['%s/%s/%s' % (re.sub('/magstr', '', self.path_output), d, f)\
				for d in os.listdir(re.sub('/magstr', '', self.path_output))   if re.search('JU',    d)\
				for f in os.listdir(re.sub('/magstr', '', self.path_output+d)) if re.search('band_', f)]
		fn_list = [fn_list[i] for i in GenGroundIdx(fn_list, exclude_f=True)] 
		fn_list = [fn for fn in fn_list if ReadFn(fn)['m'] > 0.1]

		Fp.write('%d\n' % len(fn_list))
		Fb.write('%d\n' % len(fn_list))

		for p in self.params1: Fp.write('%10s' % p)
		for p in self.params2: Fp.write('%22s' % p)
		Fp.write('\n')

		for b in ['%s%s%d' % (t, l, n) for t in ['e', 'w'] for l in self.labels for n in range(self.Nbs)]: Fb.write('%22s' % b)
		Fb.write('\n')

		band = []
		for fn in fn_list:
			fn_dict = ReadFn(fn)
			fn_dict['type'] = self.type_dict[fn_dict['type'][0]]
			for p in self.params1: Fp.write('%10.1f'  % fn_dict[p])
			for p in self.params2: Fp.write('%22.16f' % fn_dict[p])
			Fp.write('\n')
			
			with open(fn, 'r')                        as f: db = np.genfromtxt(f) - fn_dict['dntop']
			with open(re.sub('band', 'ufw', fn), 'r') as f: du = np.genfromtxt(f)
			for b in np.hstack((np.ravel(db[self.points]), np.ravel(du[self.points]))): Fb.write('%22.16f' % b)
			Fb.write('\n')

		Fp.close()
		Fb.close()

		t1 = time.time()
		print('GenBand(%s, %s) : %fs' % (pn, bn, t1-t0))

	def ReadMn(self, gn):
		mn_dict = {
			'type':  int(re.sub('AF', '', re.search('AF\d', gn).group())),
			'JU':    float(re.sub('Hz', '', re.search('Hz\d+[.]\d+', gn).group())),
			'SOC':   0,
			'N':     0,
			'U':     float(re.sub('UF', '', re.search('UF\d+[.]\d+', gn).group())),
			'n':     0,
			'm':     -1,
			'e':     0,
			'fermi': 0,
			'dntop': 0,
			'gap':   0,
		}
		
		return mn_dict

	def GenDOS0(self, eta): # k-projected DMFT DOS at high symmetry points
		eta = float(eta)

		pn = '%s/params0.txt' % (self.path_output)
		en = '%s/energy.h5' % (self.path_output)
		dn = '%s/dos_dt0_eta%.3f.h5' % (self.path_output, eta)

		t0 = time.time()		

		path_spec = '/home/Shared/BaOsO3/dmft_spec'
		d_list = ['%s/%s' % (path_spec, d) for d in os.listdir(path_spec)\
				if (re.match('oDir', d)\
				and float(re.sub('AF', '', re.search('AF\d', d).group())) > 0)] # choose only AF types (AF0 is FM)
		mn_list = ['%s/result/%s/mag.dat' % (d, os.listdir(d+'/result')[0]) for d in d_list]
		gn_list = ['%s/lattice/vdx/%s' % (d, f) for d in d_list for f in os.listdir(d+'/lattice/vdx') if re.search('kG.*ep%.2f' % eta, f)]

		if not os.path.isfile(pn):
			with open(pn, 'w') as f:
				f.write('%d\n' % len(mn_list))

				for p in self.params1: f.write('%10s' % p)
				for p in self.params2: f.write('%22s' % p)
				f.write('\n')

				for mn in mn_list:
					with open(mn, 'r') as fm: m = np.genfromtxt(fm)[-1, 2] * 2
					mn_dict = self.ReadMn(mn)
					mn_dict['m'] = m
					for p in self.params1: f.write('%10.1f'  % mn_dict[p])
					for p in self.params2: f.write('%22.16f' % mn_dict[p])
					f.write('\n')

		if not os.path.isfile(en):
			e = np.genfromtxt(gn_list[0])[:, 0]
			with h5py.File(en, 'w') as f: f.create_dataset('energy', data=e, dtype='d')

		dos = np.zeros(self.BINS_MAX * len(self.labels))
		for gn in gn_list:
			dos = np.vstack((dos, np.ravel([np.genfromtxt(re.sub('kG', 'k%s' % l, gn))[:, 1] for l in self.labels]) * 6))
		dos = np.delete(dos, 0, axis=0)

		with h5py.File(dn, 'w') as f: f.create_dataset('dos', data=dos, dtype='d')

		t1 = time.time()
		print('GenDOS0(%s) : %fs' % (dn, t1-t0))

	def Resample_(self, X, y, resamp):
		rsp = self.rsp_dict[resamp]

		X_rsp, y_rsp = rsp.fit_resample(X, y)
		X_rsp.index = X.index[rsp.sample_indices_]
		y_rsp.index = y.index[rsp.sample_indices_]

		return X_rsp, y_rsp

	def Preprocess_(self, dn, onehot=False, resamp=None):
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
		if resamp: X, y = self.Resample_(resamp, X, y)
		#idx_rsp = y.index.to_list()

		return X, y, df

	def Predict_(self, X_test, y_test, df, mc, n_dict, onehot=False):
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
		print('\n# Data : %s\n# Machine : %s\n# Resampler : %s\n# Max DOS : %f\n# Accuracy : %f (%d/%d)\n# Elapsed Time : %fs\n'\
				% (n_dict['d'], n_dict['mc'], n_dict['rsp'], np.max(df.drop(self.params, axis=1).values), acc, len(y_pred)-len(df_test[df_test['type'] != df_test['type_p']]), len(y_pred), t1-t0))
		return y_test, y_pred, y_score, df_test

	def Train(self, dn, mcn, resamp=None, tune=False):
		resamp = int(resamp)
		tune = int(tune)

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
			y_test, y_pred, y_score, df_test = self.Predict_(X_test, y_test, df, mc, {'d':dn, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'])

			t1 = time.time()
			return y_test, y_pred, y_score, df_test, mc

	def Validate(self, d1n, d2n, mcn, resamp=None):
		resamp = int(resamp)

		t0 = time.time()

		# train
		_, _, _, df1_test, mc = self.Train(d1n, mcn, resamp=resamp, tune=False)

		# validate
		X_test, y_test, df2 = self.Preprocess_(d2n, onehot=self.mc_dict[mcn]['onehot'], resamp=resamp)
		y_test, y_pred, y_score, df2_test = self.Predict_(X_test, y_test, df2, mc, {'d':d2n, 'mc':mcn, 'rsp':resamp}, onehot=self.mc_dict[mcn]['onehot'])

		t1 = time.time()
		return y_test, y_pred, y_score, df1_test, df2_test, mc

	def DrawDOS(self, dn, type, JU, N, U, ax):
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
		x = np.arange(d.shape[1] / self.Nhsp)
		y = d[i]

		ax.plot(x, y[:len(x)], label='%d dt%d_eta%.3f/JU%.1f_N%.1f_U%.1f' % (type, dtype, eta, JU, N, U))

	def DrawSpec(self, dn, type, JU, N, U):
		JU  = float(JU)
		N   = float(N)
		U   = float(U)

		dtype = re.sub('dt', '', re.search('dt\d', dn).group())
		eta   = float(re.sub('eta', '', re.search('eta\d+[.]\d+', dn).group()))

		path_band = '%s/JU%.2f_SOC%2.f' % (self.path_output, JU, 0)
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

		X = spec[:][0].reshape(db.shape[0], len(e_range))
		Y = spec[:][1].reshape(db.shape[0], len(e_range))
		Z = spec[:][2].reshape(db.shape[0], len(e_range))

		fig, ax = plt.subplots(figsize=self.figsize)
		ct = ax.contourf(X, Y, Z, levels=1000, cmap='jet')
		cb = fig.colorbar(ct)
		#cb.set_ticks(np.arange(0, 1.1, 0.2))
		ax.set_xticks([point for point, _ in self.info_path])
		ax.set_xticklabels([label for _, labels in self.info_path])
		ax.set_ylabel(r'$E - E_F$')
		ax.set_title(r'$%s$-type, $N = %.1f$ $U = %.1f$ $J/U = %.1f$' % (type[0], N, U, JU), loc='left')

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
	
	"""
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

	def GenDOS3(self, dtype, eta): # full dos
		df0, dn, e_range, e_label, w_label = self.OpenBand(dtype, eta)

		f_label = ['x%d' % i for i in range(self.max_bins)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		dir_list = [self.path_output + dir for dir in os.listdir(self.path_output) if not re.search('magstr', dir)]	

		for dir in dir_list:
			fn_list = [dir +'/'+ fn for fn in os.listdir(dir) if re.match('dos', fn)]
			gi_list = GenGroundIdx(fn_list, dtype='dos', exclude_f=True)

			for fn in [fn_list[i] for i in gi_list]:
				fn_dict = ReadFn(fn, dtype='dos')
				fn_dict['type'] = fn_dict['type'][0]

				dos = []
				for p in self.params: dos.append(fn_dict[p])

				if fn_dict['U'] > self.tol_U and fn_dict['m'] > tol:
					f = open(fn, 'r')
					d = np.genfromtxt(f)
					f.close()

					d = np.transpose(d)
					d = pd.DataFrame(d[1:, :], columns=d[0, :]).sum()
					x = [e - fn_dict['dntop'] for e in d.index]

					dos += list(np.interp(e_range, d.index, d))

					data = pd.DataFrame([dos], columns=x_label)
					df = pd.concat([df, data], sort=False)

		return df, dn

	def GenDOS(self, dtype, eta):
		t0 = time.time() 

		eta = float(eta)
		d_dict = {
			'1': self.GenDOS1,
			'2': self.GenDOS2,
			'3': self.GenDOS3,
			'4': self.GenDOS4,
			'5': self.GenDOS5,
		}

		df, dn = d_dict[dtype[0]](dtype, eta)

		df = df.reset_index(drop=True)
		df.to_csv(dn, sep=',')

		t1 = time.time()
		print('GenDOS(%s) : %fs' % (dn, t1-t0))
	"""
