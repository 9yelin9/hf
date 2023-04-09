# hf3/dmftt.py : dmft tools, read dmft/oDir_*/*

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.patches import Circle

class DMFTT:
	def __init__(self, AF, U, J, Hz, dmft_dir='dmft'):
		self.NB = 9 # # of bath orbitals
		self.NC = 3 # # of impurity orbitals
		self.Ni = 8 # # of impurities
		self.dmft_dir = dmft_dir
		self.p_dict = {
			'U' : float(U),
			'AF': int(AF),
			'Hz': float(Hz),
			'J' : float(J),
		}
		self.title = ' '.join(['%s=%.2f' % (p, self.p_dict[p]) for p in ['AF', 'U', 'J', 'Hz']])
		self.fname = '_'.join(['%s%.2f'  % (p, self.p_dict[p]) for p in ['AF', 'U', 'J', 'Hz']])
		
		self.oDir = self.GetoDir()
		self.path_oDir = '%s/%s' % (self.dmft_dir, self.oDir)
		self.path_res = '%s/result/u%.3f' % (self.path_oDir, self.p_dict['U'])
		if not os.path.isdir(self.path_res): print('Still in progress...')

		self.type_dict = {
			'A': [0, 2, 4, 6],
			'C': [0, 3, 4, 7],
			'G': [0, 3, 5, 6]
		}

		self.markers = ['s', 'o', '>', '<', '^', 'v', '*', '.']
		plt.rcParams.update({'font.size': 30})
	
	def GetoDir(self):
		oDirs = [d for d in os.listdir(self.dmft_dir) if re.search('UF%.2f.*AF%d_Hz%.3f.*J%.3f' % tuple(self.p_dict.values()), d)]
		if len(oDirs) > 1:
			for i, oDir in enumerate(oDirs): print('%2d)'%i, oDir)
			x = input('Which oDir? : ')
			return oDirs[int(x)]
		else:
			return oDirs[0]

	def GetParams_(self):
		df = pd.read_csv('%s/parameters_u%.2f.dat' % (self.path_oDir, self.p_dict['U']), sep='\s+', skip_blank_lines=False)
		return df

	def GetMag_(self):
		df = pd.read_csv('%s/mag.dat' % self.path_res, sep='\s+',\
				names=['count', 'm_total'] + ['m_%d%s' % (i, orb) for i in range(self.Ni) for orb in ['t2g', 'xy', 'yz', 'zx']])
		return df

	def GetNext_(self):
		next_dict = {}
		for i in range(self.Ni): next_dict[str(i)] = {'ep':[], 'v':[]}

		with open('%s/next.mak' % self.path_oDir, 'r') as f:
			for _ in range(2): f.readline()
			for i in range(self.Ni):
				for _ in range(2): f.readline()
				next_dict[str(i)]['ep'] = np.fromstring(f.readline(), sep=' ')
				f.readline()
				for _ in range(2*self.NC): next_dict[str(i)]['v'].append(np.fromstring(f.readline(), sep=' '))
				for _ in range(2): f.readline()

		return next_dict

	def PrintParams(self, dtype, cnt, imp):
		cnt = int(cnt)
		df = self.GetParams_()

		for dt in dtype.split(','):
			print(df.loc[cnt:, ['%s_%d'%(dt, i) for i in imp]], end='\n\n')

	def PrintMag(self, dtype, cnt, imp):
		cnt = int(cnt)
		df = self.GetMag_()

		for dt in dtype.split(','):
			print(df.loc[cnt:, ['m_%d%s'%(i, dt) for i in imp]], end='\n\n')
		
		sign = 1 if df.loc[df.index[-1], 'm_0t2g'] > 0 else -1
		magord = np.where(df.loc[df.index[-1], ['m_%dt2g' % i for i in range(self.Ni)]]*sign > 0)[0].tolist()
		for key, val in self.type_dict.items():
			if magord == val: print('=> %s-type' % key)

	def DrawIter(self, df, cnt, imp, ax):
		x = df.loc[:, '#count']
		y = df.loc[:, ['iter_%d'%i for i in range(self.Ni)]]

		ax2 = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

		for i in imp:
			ax.plot(       x,    y.iloc[:, i], ms=10, c=cm.tab10(i), marker=self.markers[i], label=i)
			ax2.plot(x[cnt:], y.iloc[cnt:, i], ms=10, c=cm.tab10(i), marker=self.markers[i], label=i)

		ax.grid()
		ax.set_xticks(x, labels=['%d'%xi if xi in [x.min(), x.max()] else '' for xi in x])
		ax.set_xlabel('count')
		ax.set_ylabel('iter')

		ax2.grid()
		ax2.set_ylim([0, y.iloc[cnt:, :].max().max()+1])
		ax2.set_xticks(x[cnt:])
		ax2.set_yticks(np.linspace(1, y.iloc[cnt:, :].max().max(), num=5, dtype='i'))
		ax2.tick_params(axis='both', labelsize='small')

		if not cnt: ax2.set_visible(False)

	def DrawNpos(self, df, cnt, imp, ax):
		x = df.loc[:, '#count']
		y = df.loc[:, ['npos_%d'%i for i in range(self.Ni)]]

		for i in imp:
			ax.plot(x, y.iloc[:, i], ms=10, c=cm.tab10(i), marker=self.markers[i], label=i)

		ax.grid()
		ax.set_xticks(x, labels=['%d'%xi if xi in [x.min(), x.max()] else '' for xi in x])
		ax.set_xlabel('count')
		ax.set_ylabel('npos')

	def ShowParams(self, dtype, cnt, imp):
		cnt = int(cnt)
		dtype = dtype.split(',')
		df = self.GetParams_()

		draw_dict = {
			'iter': self.DrawIter,
			'npos': self.DrawNpos,
		}

		fig, ax = plt.subplots(1, len(dtype), figsize=(3+6*len(dtype), 6), constrained_layout=True)
		axs = np.ravel([ax])
		
		for dt, ax in zip(dtype, ax):
			draw_dict[dt](df, cnt, imp, ax)
		
		fig.suptitle(self.title, fontsize='medium')
		plt.legend(fontsize='x-small', title='Ni', title_fontsize='x-small', loc='upper left',\
				   labelspacing=0.02, handletextpad=0.3, bbox_to_anchor=(1.2+0.02*len(dtype), 1.035))
		path_fig = '%s/diagram/%s_%s.%s' % (self.dmft_dir, self.fname, ''.join(dtype), 'png')
		fig.savefig(path_fig)
		print('Figure saved at %s' % path_fig)
		plt.show()

	def ShowNext(self, imp):
		next_dict = self.GetNext_()
		di = 2

		fig, ax = plt.subplots(figsize=(15, 9), constrained_layout=True)

		for i in imp:
			eps = np.array(next_dict[str(i)]['ep'])
			v   = np.array(next_dict[str(i)]['v'])
			i = di * i

			ans = []
			for j in range(2*self.NC):
				for k, ep in enumerate(eps):
					an = ax.annotate('', xy=(i+v[j, 2*k], ep+v[j, 2*k+1]), xytext=(i, ep),\
									 arrowprops=dict(arrowstyle='-', fc=cm.Paired(j), ec=cm.Paired(j)))
					c = Circle((i, ep), np.linalg.norm([v[j, 2*k], v[j, 2*k+1]]), color=cm.Paired(j), alpha=0.3, clip_on=False)
					ax.add_patch(c)
				ans.append(an.arrow_patch)
			ax.scatter([i for _ in range(2*self.NB)], eps, marker='.', s=1, color='w')

		ax.axhline(y=0, ls='--', lw=2, color='dimgrey')
		ax.grid()
		ax.set_xlabel('Impurity index')
		ax.set_ylabel(r'$E_\mathrm{g}^\mathrm{bath}$')
		ax.set_xticks(range(0, di*self.Ni, di), labels=range(self.Ni))
		ax.set_xlim([-0.5*di, di*self.Ni-0.5*di])

		fig.suptitle(self.title, fontsize='medium')
		plt.legend(ans, ['%s%s' % (ob, sp) for ob in ['xy', 'yz', 'zx'] for sp in [r'$\uparrow$', r'$\downarrow$']],\
				   fontsize='x-small', title='V', title_fontsize='x-small', loc='upper left',\
				   labelspacing=0.02, handletextpad=0.3, bbox_to_anchor=(1, 1.02))
		path_fig = '%s/diagram/%s_next.%s' % (self.dmft_dir, self.fname, 'png')
		fig.savefig(path_fig)
		print('Figure saved at %s' % path_fig)
		plt.show()
