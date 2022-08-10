# pyhf3/draw.py : draw diagrams of 3 band models

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

class Draw:
	def __init__(self, material='boo', JU=0.0, SOC=0.0, type='f', basis=[], path=[], path_label=[]):
		self.material = material
		self.JU = float(JU)
		self.SOC = float(SOC)
		self.type = type
		self.bases = basis
		self.path = path
		self.path_label = path_label

		self.K = 16
		self.obt = 3
		self.basis = basis[0] if re.search('f', type) else basis[1]

		self.title = '%s K=%d J/U=%.2f SOC=%.2f\n' % (type, self.K, self.JU, self.SOC)
		self.colors=['b', 'g', 'r']
		self.labels=['xy', 'yz', 'zx']
	
	def DrawBand(self, N, U, ax, dir):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(dir) if re.search('band_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi\d+[.]\d+', fs).group()))

		f = open(dir + fs, 'r')
		arr = np.genfromtxt(f)
		f.close()

		e_min = np.min(arr[:, 1:]) - fermi - 0.1
		e_max = np.max(arr[:, 1:]) - fermi + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, self.basis+1):
			ax.plot(arr[:, 0], arr[:, i] - fermi, color='tab:blue')

		ax.grid(True)
		ax.set_xticks(self.path)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs 
	
	def DrawBandUF(self, N, U, ax, dir):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(dir) if re.search('band_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi\d+[.]\d+', fs).group()))

		f = open(dir + fs, 'r')
		arr = np.genfromtxt(f)
		f.close()

		e_min = np.min(arr[:, self.bases[1]+1:]) - fermi - 0.1
		e_max = np.max(arr[:, self.bases[1]+1:]) - fermi + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, self.basis+1):
			#points = np.array(arr[:, 0], arr[:, self.bases[1]+i]).T.reshape(-1, 1, 2)
			#segments = np.concatenate([points[:-1], points[1:]], axis=1)
			#ax.plot(arr[:, 0], arr[:, self.bases[1]+i] - fermi, color='tab:blue', lw=arr[:, i])
			lc = LineCollection(segments, linewidths=arr[:, i], color='tab:blue')
			ax.add_collection(lc)

		ax.grid(True)
		ax.set_xticks(self.path)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs 

	def DrawDos(self, N, U, ax, dir, e_min=None, e_max=None):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(dir) if re.search('dos_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi\d+[.]\d+', fs).group()))
		title = dir + fs

		f = open(dir + fs, 'r')
		arr = np.genfromtxt(f)
		f.close()

		if e_min == None: e_min = np.min(arr[:, 0]) - fermi
		if e_max == None: e_max = np.max(arr[:, 0]) - fermi

		dos_max = np.max(arr[:, 1:])

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.obt+1): 
			ax.plot(arr[:, i] / dos_max, arr[:, 0] - fermi, lw=abs(i-5), label=self.labels[i-1], color=self.colors[i-1])
			ax.plot(-arr[:, i+self.obt] / dos_max, arr[:, 0] - fermi, lw=abs(i-5), color=self.colors[i-1])

		ax.grid(True)
		ax.set_xticks([-1, 0, 1])
		ax.set_xticklabels([1, 0, 1])
		ax.set_xlabel('Density of states')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.legend()
	
	def DrawBandDos(self, N, U, is_unfold=0):
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=[16, 8])

		if not is_unfold:
			dir = '%s/output/K%d_JU%.2f_SOC%.2f/' % (self.material, self.K, self.JU, self.SOC)
			e_min, e_max, fs = self.DrawBand(N, U, ax[0], dir)
			self.DrawDos(N, U, ax[1], dir, e_min, e_max)
		else:
			dir = '%s/output/K%d_JU%.2f_SOC%.2f_unfold/' % (self.material, self.K, self.JU, self.SOC)
			e_min, e_max, fs = self.DrawBandUF(N, U, ax[0], dir)
			self.DrawDos(N, U, ax[1], dir, e_min, e_max)

		fs_s = fs.split(sep='_')
		fs_s = [s for s in fs_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		title = self.title
		for s in fs_s:
			s_name = re.search('[a-zA-Z]+', s).group()
			s_val = re.search('\d+[.]\d+', s).group()
			title += '%s=%s ' % (s_name, s_val)

		plt.suptitle(title)
		fig.tight_layout()
		fig.savefig('%s/%s' % (re.sub('output/', 'diagram/', dir), re.sub('txt', 'png', fs)))
		plt.show()

	def DrawPhase(self):
		val = [[float(re.sub('_N', '', re.search('_N\d+[.]\d+', f).group())),\
				float(re.sub('_U', '', re.search('_U\d+[.]\d+', f).group())),\
				float(re.sub('_m', '', re.search('_m[-]?\d+[.]\d+', f).group()))] for f in os.listdir(dir) if re.match('band_%s' % (self.type), f)]
		val = np.array(val)

		x = np.unique(val[:, 0])
		y = np.unique(val[:, 1])

		X, Y = np.meshgrid(x, y)
		Z = np.zeros((len(y), len(x)))

		for v in val:
			"""
			f = open(v[4], 'r')
			e = np.genfromtxt(f)[:, 1:].astype('float')
			f.close()

			e = np.abs(e - float(v[3]))
			if np.where(e < 1e-3): ins += [[v[0], v[1]]]
			"""

			for i, xi in enumerate(x):
				for j, yj in enumerate(y):
					if v[0] == xi and v[1] == yj: Z[j][i] = v[2]

		fig, ax = plt.subplots()
		#cs = ax.contour(X, Y, Z, levels=10, cmap='binary_r')
		cb = ax.contourf(X, Y, Z, levels=10)
		#plt.clabel(cs)
		plt.colorbar(cb, label='Magnetization (M)')
		#ax.grid(True, color='black')
		ax.set_xlabel('Occupation (N)')
		ax.set_ylabel('Interaction (U)')
		plt.suptitle(dir + self.type)	
		fig.tight_layout()
		fig.savefig('diagram/%s/phase_%s' % (re.sub('output/', '', dir), self.type))
		plt.show()
