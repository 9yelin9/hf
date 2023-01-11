# pyhf3/magstr.py : predict magnetic structure of Hartree-Fock approximated 3-band Hubbard model

import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

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

from .read import ReadFn, ReadFnDMFT, GenGroundIdx

class MagStr:
	def __init__(self, path_output, info_path, info_cell):
		self.path_output = path_output

		info_path_nodup = [] # drop duplicated hsp
		for point, label in info_path:
			if label not in [label for point, label in info_path_nodup]:
				info_path_nodup.append([point, label])
		self.info_path = info_path_nodup

		#self.info_cell = info_cell
		#self.Ni = self.info_cell['a'][0]
		#self.Nc = self.info_cell['a'][1]
		#self.Nb = self.Ni * self.Nc * 2

		self.Nbs = 12
		self.Nhsp = 4
		self.BINS_MAX = 8193
		self.tol_U = 0.1

		self.type_dict = {'a':'1', 'c':'2', 'g':'3'}

		self.params1 = ['type', 'JU', 'N', 'U']
		self.params2 = ['m', 'dntop', 'gap']
		self.params = self.params1 + self.params2

		self.mc_dict = {
			'rf': RandomForestClassifier(random_state=1, n_jobs=1),
			'lr': LogisticRegression(random_state=1, n_jobs=1, multi_class='multinomial', max_iter=10000),
			'xgb': XGBClassifier(random_state=1, nthread=1),
			'lgb': LGBMClassifier(random_state=1, n_jobs=1),
			'cat': CatBoostClassifier(random_state=1, silent=True),
		}

		self.figsize=(11, 6)

		os.makedirs('%s/magstr/' % path_output, exist_ok=True)

	def GenBand(self):
		pn = '%s/magstr/params.txt' % (self.path_output)
		bn = '%s/magstr/band.txt' % (self.path_output)
		Fp = open(pn, 'w')
		Fb = open(bn, 'w')

		t0 = time.time()

		fn_list = ['%s/%s/%s' % (self.path_output, d, f)\
				for d in os.listdir(self.path_output)   if re.match('JU',    d)\
				for f in os.listdir(self.path_output+d) if re.match('band_', f)]
		fn_list = [fn_list[i] for i in GenGroundIdx(fn_list, exclude_f=True)] 

		Fp.write('%d\n' % len(fn_list))
		Fb.write('%d\n' % len(fn_list))

		points = [point for point, _ in self.info_path]
		labels = [label for _, label in self.info_path]

		for p in self.params1: Fp.write('%10s' % p)
		for p in self.params2: Fp.write('%22s' % p)
		Fp.write('\n')

		for b in ['%s%s%d' % (t, l, n) for t in ['e', 'w'] for l in labels for n in range(self.Nbs)]: Fb.write('%22s' % b)
		Fb.write('\n')

		band = []
		for fn in fn_list:
			fn_dict = ReadFn(fn)
			for key, val in self.type_dict.items(): fn_dict['type'] = fn_dict['type'][0].replace(key, val)
			
			for p in self.params1: Fp.write('%10s'    % fn_dict[p])
			for p in self.params2: Fp.write('%22.16f' % fn_dict[p])
			Fp.write('\n')
			
			with open(fn, 'r')                        as fb: db = np.genfromtxt(fb) - fn_dict['dntop']
			with open(re.sub('band', 'ufw', fn), 'r') as fu: du = np.genfromtxt(fu)
			for b in np.hstack((np.ravel(db[points]), np.ravel(du[points]))): Fb.write('%22.16f' % b)
			Fb.write('\n')

		Fp.close()
		Fb.close()

		t1 = time.time()
		print('GenBand(%s, %s) : %fs' % (pn, bn, t1-t0))

	def GenDOSDMFT(self, dtype, eta):
		df0, dn, e_range, e_label, w_label = self.OpenBand(dtype, eta)

		path_spec = '/home/Shared/BaOsO3/dmft_spec/'
		dir_list = [path_spec + dir for dir in os.listdir(path_spec)\
				if (re.match('oDir', dir)\
				and float(re.sub('AF', '', re.search('AF\d', dir).group())) > 0\
				and float(re.sub('_D', '', re.search('_D\d+[.]\d+', dir).group())) < 0.1\
				and len(os.listdir(path_spec + dir)) > 1)]

		n_hsp = len(e_label) // self.Nbs

		p_label = []
		f_label = []
		e_label_sp = [[] for _ in range(n_hsp)]
		w_label_sp = [[] for _ in range(n_hsp)]

		for i in range(n_hsp):
			p_label.append(re.sub('e', '', e_label[self.Nbs*i].split('_')[0]))
			for j in range(self.max_bins):
				f_label.append('%s_%d' % (p_label[i], j))
			for j in range(self.Nbs):
				e_label_sp[i].append(e_label[self.Nbs*i + j])
				w_label_sp[i].append(w_label[self.Nbs*i + j])
		x_label = self.params + f_label

		df = pd.DataFrame()

		for dir in dir_list:
			dir_lat = dir + '/lattice/vdx/'
			G_list = [dir_lat +'/'+ fn for fn in os.listdir(dir_lat) if re.search('_kG.*_ep%.2f' % eta, fn)]

			mag_dat = open(dir + '/result/%s/mag.dat' % os.listdir(dir + '/result/')[0], 'r')
			m = np.genfromtxt(mag_dat)[-1, 2] * 2
			mag_dat.close()

			for G in G_list:
				fn_list = [re.sub('_kG', '_k%s' % p, G) for p in p_label]
				fn_dict = ReadFnDMFT(G)

				dos = []
				for p in self.params: dos.append(fn_dict[p])

				for fn in fn_list:
					f = open(fn, 'r')
					d = np.genfromtxt(f)[:, :2]
					f.close()

					itp = interpolate.interp1d(d[:, 0], d[:, 1] * 6, kind='cubic', fill_value='extrapolate')
					dos += list(itp(e_range))

				data = pd.DataFrame([dos], columns=x_label)
				data['m'] = m
				df = pd.concat([df, data], sort=False)

			df = df.replace({'0':'f', '1':'a', '2':'c', '3':'g'})

		return df, dn

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

	def Preprocess(self, df, mcn, rspn):
		X = df.drop(self.params, axis=1)
		y = df['type']

		if mcn == 'xgb': y = pd.get_dummies(y)
		if rspn != '': X, y = self.Resample(rspn, X, y)
		#idx_rsp = y.index.to_list()

		return X, y

	def Predict(self, mc, mcn, df, X_test, y_test):
		y_pred = mc.predict(X_test)
		y_score = mc.predict_proba(X_test)

		if mcn == 'xgb':
			t_dict = {}
			for i, col in enumerate(y_test.columns): t_dict[str(i)] = col
			y_test = y_test.idxmax(axis=1)
			y_pred = np.array([t_dict[str(np.argmax(y))] for y in y_pred])

		# misclassified and well-classified data
		idx_m = []
		idx_w = []	
		type_pred_m = []
		type_pred_w = []
		score_m = []
		score_w = []

		for i, y_true in enumerate(y_test):
			if(y_true != y_pred[i]):
				idx_m.append(y_test.index[i])
				type_pred_m.append(y_pred[i])
				score_m.append(list(map(str, y_score[i])))
			else:
				idx_w.append(y_test.index[i])
				type_pred_w.append(y_pred[i])
				score_w.append(list(map(str, y_score[i])))

		df_m = df.loc[idx_m, :]
		df_w = df.loc[idx_w, :]
		df_m['type_pred'] = type_pred_m
		df_w['type_pred'] = type_pred_w
		df_m['score'] = score_m
		df_w['score'] = score_w

		return y_test, y_pred, y_score, df_m, df_w

	def Train(self, dn, mcn, rspn, tuning=False):
		t0 = time.time()

		hp_dict = {
			'rf': {
				'random_state': 1,
			}
		}
		mc = mc_dict[mcn[0]]

		# prepare data
		df = pd.read_csv(dn, sep=',', index_col=0)

		X, y = self.Preprocess(df, mcn[0], rspn)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=1)
		
		# train
		if len(mcn) < 2:
			mc.fit(X_train, y_train)

			y_test, y_pred, y_score, df_m, df_w = self.Predict(mc, mcn[0], df, X_test, y_test)
			acc = accuracy_score(y_test, y_pred)
			print('acc = %d/%d = %f' % (len(y_pred)-len(df_m), len(y_pred), acc))

			t1 = time.time()
			print('Train(%s-%s-%s) : %fs' % (mcn[0], rspn, dn, t1-t0))

			return mc, df, df_m, df_w, y_test, y_pred, y_score

		# tune
		else:
			cv = StratifiedKFold(random_state=1)
			tn = RandomizedSearchCV(mc, hp_dict[mcn[0]], cv, random_state=1)
			res = tn.fit(X_train, y_train)

			t1 = time.time()
			print('Tune(%s-%s-%s) : %fs' % (mcn[0], rspn, dn, t1-t0))

			return res

	def Validate(self, d1n, d2n, mcn, rspn):
		t0 = time.time()

		mc, df1, _, _, _, _, _ = self.TrainOrTune(d1n, mcn, rspn)

		df2 = pd.read_csv(d2n, index_col=0)
		X_test, y_test = self.Preprocess(df2, mcn, rspn)

		y_test, y_pred, y_score, df2_m, df2_w = self.Predict(mc, mcn, df2, X_test, y_test)
		acc = accuracy_score(y_test, y_pred)
		print('acc2 = %d/%d = %f' % (len(y_pred)-len(df2_m), len(y_pred), acc))

		t1 = time.time()
		print('Predict(%s-%s) : %fs' % (mcn, rspn, t1-t0))

		return mc, df1, df2, df2_m, df2_w, y_test, y_pred, y_score
		
	def DrawSpec(self, dtype, eta, JU, type, N, U):
		eta = float(eta)
		JU = float(JU)
		N = float(N)
		U = float(U)

		dir = '%s/JU%.2f_SOC%.2f/' % (self.path_output, JU, 0)
		os.makedirs(re.sub('output', 'diagram', dir), exist_ok=True)
		bn = [dir + f for f in os.listdir(dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		un = re.sub('band_', 'ufw_', bn)

		b = open(bn, 'r')
		u = open(un, 'r')
		band = np.genfromtxt(b)
		ufw  = np.genfromtxt(u)
		b.close()
		u.close()

		bn_dict = ReadFn(bn)
		fermi = bn_dict['dntop']
		band = band - fermi
		e_range = np.linspace(band.min(), band.max(), 128)

		# options
		if   re.search('f', dtype): w_fermi = [1 if e < 0 else 0 for e in e_range]
		elif re.search('F', dtype): w_fermi = [1 if (e < 0 and e > -1) else 0 for e in e_range]
		else:                       w_fermi = [1 for _ in e_range]

		if re.search('b', dtype): eta_brd = [0.1 + e/e_range.min() if e < 0 else 0.1 + e/e_range.max() for e in e_range]
		else:                     eta_brd = [eta for _ in e_range]

		Z = pd.DataFrame()
		
		for i in range(band.shape[0]):
			for e, wf, etab in zip(e_range, w_fermi, eta_brd):
				spec = [i, e]
				spec0 = 0

				for j in range(band.shape[1]):
					spec0 += (etab / ((e - band[i][j])**2 + etab**2)) * ufw[i][j] * wf
				spec.append(spec0 / np.pi)

				data = pd.DataFrame([spec], columns=['path', 'e', 'spec'])
				Z = pd.concat([Z, data], sort=False)

		X = Z['path'].values.reshape(band.shape[0], len(e_range))
		Y = Z['e'].values.reshape(band.shape[0], len(e_range))
		Z = Z['sp(0, 11),iec'].values.reshape(band.shape[0], len(e_range))

		fig, ax = plt.subplots(figsize=self.figsize)
		ct = ax.contourf(X, Y, Z, levels=len(e_range), cmap='jet')
		cb = fig.colorbar(ct)
		#cb.set_ticks(np.arange(0, 1.1, 0.2))
		ax.set_xticks([0, 198, 396, 677, 1023])
		ax.set_xticklabels([r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R'])
		ax.set_ylabel(r'$E - E_F$')
		ax.set_title(r'$%s$-type, $N = %.1f$ $U = %.1f$ $J/U = %.1f$' % (type[0], N, U, JU), loc='left')

		fig.tight_layout()
		fig.savefig('%s' % re.sub('output', 'diagram', re.sub('band_', 'spec%s_' % dtype, re.sub('txt', 'png', bn))))
		plt.show()

	def DrawConfMat(self, y_test, y_pred, title=''):
		from sklearn.metrics import ConfusionMatrixDisplay

		fig, ax = plt.subplots(figsize=self.figsize)
		ConfusionMatrixDisplay.from_predictions(y_test, y_pred, normalize='true', cmap='Blues', values_format='.2f', ax=ax)
		ax.set_title(title, loc='left')
		plt.show()

	def DrawROC(self, y_test, y_pred, title=''):
		from sklearn.preprocessing import label_binarize
		from sklearn.metrics import roc_curve, auc

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

	def OpenBand(self, dtype, eta):
		dt = dtype.split(sep='-') # dtype = 'dtype-ptype-pnum'
		dn = '%s/magstr/dos_dt%s_eta%.3f.csv' % (self.path_output, dtype, eta)

		df = pd.read_csv('%s/magstr/band_pt%s_pn%s.csv' % (self.path_output, dt[1], dt[2]), index_col=0)
		df = df[df['U'] > self.tol_U]

		e_label = [v for v in df.columns if re.match('e', v)]
		w_label = [re.sub('e', 'w', v) for v in e_label]

		for e in e_label: df[e] = df[e] - df['dntop']
		df = df.reset_index(drop=True)

		e = df.loc[:, e_label].values.flatten()
		e_range = np.linspace(min(e), max(e), self.max_bins)

		return df, dn, e_range, e_label, w_label

	def GenDOS1(self, dtype, eta): # partial dos
		df0, dn, e_range, e_label, w_label = self.OpenBand(dtype, eta)

		f_label = ['x%d' % i for i in range(self.max_bins)]
		x_label = self.params + f_label

		df = pd.DataFrame()

		for i in df0.index:
			dos = list(df0.loc[i, self.params])

			for e in e_range:
				dos0 = 0
				for el, wl in zip(e_label, w_label):
					dos0 += (eta / ((e - df0.loc[i, el])**2 + eta**2)) * df0.loc[i, wl]
				dos.append(dos0 / np.pi)

			data = pd.DataFrame([dos], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		return df, dn

	def GenDOS2(self, dtype, eta): # hsp dos
		df0, dn, e_range, e_label, w_label = self.OpenBand(dtype, eta)

		# options
		if   re.search('f', dtype): w_fermi = [1 if e < 0 else 0 for e in e_range]
		elif re.search('F', dtype): w_fermi = [1 if (e < 0 and e > -1) else 0 for e in e_range]
		else:                       w_fermi = [1 for _ in e_range]

		if re.search('b', dtype): eta_brd = [0.1 + e/e_range.min() if e < 0 else 0.1 + e/e_range.max() for e in e_range]
		else:                     eta_brd = [eta for _ in e_range]

		n_hsp = len(e_label) // self.Nbs

		p_label = []
		f_label = []
		e_label_sp = [[] for _ in range(n_hsp)]
		w_label_sp = [[] for _ in range(n_hsp)]

		for i in range(n_hsp):
			p_label.append(re.sub('e', '', e_label[self.Nbs*i].split('_')[0]))
			for j in range(self.max_bins):
				f_label.append('%s_%d' % (p_label[i], j))
			for j in range(self.Nbs):
				e_label_sp[i].append(e_label[self.Nbs*i + j])
				w_label_sp[i].append(w_label[self.Nbs*i + j])
		x_label = self.params + f_label

		df = pd.DataFrame()

		for i in df0.index:
			dos = list(df0.loc[i, self.params])

			for hsp in range(n_hsp):
				for e, wf, etab in zip(e_range, w_fermi, eta_brd):
					dos0 = 0
					for el, wl in zip(e_label_sp[hsp], w_label_sp[hsp]):
						dos0 += (etab / ((e - df0.loc[i, el])**2 + etab**2)) * df0.loc[i, wl] * wf
					dos.append(dos0 / np.pi)

			data = pd.DataFrame([dos], index=[i], columns=x_label)
			df = pd.concat([df, data], sort=False)

		return df, dn
	
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

	def GenDOS4(self, dtype, eta): # DMFT dos
		df0, dn, e_range, e_label, w_label = self.OpenBand(dtype, eta)

		path_spec = '/home/Shared/BaOsO3/dmft_spec/'
		dir_list = [path_spec + dir for dir in os.listdir(path_spec)\
				if (re.match('oDir', dir)\
				and float(re.sub('AF', '', re.search('AF\d', dir).group())) > 0\
				and float(re.sub('_D', '', re.search('_D\d+[.]\d+', dir).group())) < 0.1\
				and len(os.listdir(path_spec + dir)) > 1)]

		n_hsp = len(e_label) // self.Nbs

		p_label = []
		f_label = []
		e_label_sp = [[] for _ in range(n_hsp)]
		w_label_sp = [[] for _ in range(n_hsp)]

		for i in range(n_hsp):
			p_label.append(re.sub('e', '', e_label[self.Nbs*i].split('_')[0]))
			for j in range(self.max_bins):
				f_label.append('%s_%d' % (p_label[i], j))
			for j in range(self.Nbs):
				e_label_sp[i].append(e_label[self.Nbs*i + j])
				w_label_sp[i].append(w_label[self.Nbs*i + j])
		x_label = self.params + f_label

		df = pd.DataFrame()

		for dir in dir_list:
			dir_lat = dir + '/lattice/vdx/'
			G_list = [dir_lat +'/'+ fn for fn in os.listdir(dir_lat) if re.search('_kG.*_ep%.2f' % eta, fn)]

			mag_dat = open(dir + '/result/%s/mag.dat' % os.listdir(dir + '/result/')[0], 'r')
			m = np.genfromtxt(mag_dat)[-1, 2] * 2
			mag_dat.close()

			for G in G_list:
				fn_list = [re.sub('_kG', '_k%s' % p, G) for p in p_label]
				fn_dict = ReadFnDMFT(G)

				dos = []
				for p in self.params: dos.append(fn_dict[p])

				for fn in fn_list:
					f = open(fn, 'r')
					d = np.genfromtxt(f)[:, :2]
					f.close()

					itp = interpolate.interp1d(d[:, 0], d[:, 1] * 6, kind='cubic', fill_value='extrapolate')
					dos += list(itp(e_range))

				data = pd.DataFrame([dos], columns=x_label)
				data['m'] = m
				df = pd.concat([df, data], sort=False)

			df = df.replace({'0':'f', '1':'a', '2':'c', '3':'g'})

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
