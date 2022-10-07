# pyhf3/draw.py : draw diagrams of 3 band models

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from .read import ReadFs

class Draw:
	def __init__(self, input_path, output_path, path_info, type_info, JU, SOC, tol):
		self.input_path = input_path
		self.path_point = [path_point for path_point, path_label in path_info]
		self.path_label = [path_label for path_point, path_label in path_info]
		self.type_info = type_info
		self.JU = float(JU)
		self.SOC = float(SOC)
		self.tol = tol

		self.obt = 3

		self.dir = '%s/JU%.2f_SOC%.2f/' % (output_path, self.JU, self.SOC)
		self.colors=['tab:blue', 'tab:green', 'tab:red']
		self.labels=['xy', 'yz', 'zx']
		self.markers = ['s', 'o', '^']
		self.lss = ['-', '--', '-.']

		os.makedirs(re.sub('output', 'diagram', self.dir), exist_ok=True)
	
	def DrawBandTest(self, type):
		basis = self.type_info[type[0]]
		fs = [f for f in os.listdir(self.input_path) if re.search('band_%s.txt' % (type), f)][0]

		f = open(self.input_path + fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		fig, ax = plt.subplots()

		for i in range(1, basis+1):
			ax.plot(data[:, 0], data[:, i], color='tab:blue')

		ax.grid(True)
		ax.set_xticks(self.path_point)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_title(self.input_path + fs)
		plt.show()
	
	def DrawBand(self, type, N, U, ax):
		N = float(N)
		U = float(U)

		basis = self.type_info[type[0]]
		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fs_dict = ReadFs(fs)
		fermi = fs_dict['fermi']

		f = open(fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		e_min = np.min(data[:, 1:]) - fermi - 0.1
		e_max = np.max(data[:, 1:]) - fermi + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, basis+1):
			ax.plot(data[:, 0], data[:, i] - fermi, color='tab:blue')

		ax.grid(True)
		ax.set_xticks(self.path_point)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs 
	
	def DrawBandUnfold(self, type, N, U, ax):
		N = float(N)
		U = float(U)

		basis = self.type_info[type[0]]
		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fs_dict = ReadFs(fs)
		fermi = fs_dict['fermi']

		f_band = open(fs, 'r')
		f_ufw  = open(re.sub('band_', 'ufw_', fs), 'r')

		data_e = np.genfromtxt(f_band)
		data_w = np.genfromtxt(f_ufw)

		f_band.close()
		f_ufw.close()

		e_min = np.min(data_e[:, 1:]) - fermi - 0.1
		e_max = np.max(data_e[:, 1:]) - fermi + 0.1

		norm = plt.Normalize(0, 1)
		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, basis+1):
			ax.plot(data_e[:, 0], data_e[:, i] - fermi, linewidth=1, color='gray', alpha=0.5)
			ax.scatter(data_e[:, 0], data_e[:, i] - fermi, c=data_w[:, i], s=(10*(data_w[:, i]))**2, cmap='Blues', norm=norm)

		"""
		for i in range(1, basis+1):
			points = np.array([data_e[:, 0], data_e[:, i] - fermi]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)
			norm = plt.Normalize(0, 1)
			lc = LineCollection(segments, cmap='PuBu', norm=norm)
			lc.set_array(data_w[:, i])
			lc.set_linewidth(np.log(0.5+10*data_w[:, i]))
			line = ax.add_collection(lc)
		"""

		ax.grid(True)
		ax.set_xticks(self.path_point)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs

	def DrawDos(self, type, N, U, ax, e_min=None, e_max=None):
		N = float(N)
		U = float(U)

		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('dos_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fs_dict = ReadFs(fs, 'dos')
		fermi = fs_dict['fermi']

		f = open(fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		if e_min == None: e_min = np.min(data[:, 0]) - fermi
		if e_max == None: e_max = np.max(data[:, 0]) - fermi

		dos_max = np.max(data[:, 1:])

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.obt+1): 
			ax.plot(data[:, i] / dos_max,           data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors[i-1], label=self.labels[i-1])
			ax.plot(-data[:, i+self.obt] / dos_max, data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors[i-1])

		ax.grid(True)
		ax.set_xticks([-1, 0, 1])
		ax.set_xticklabels([1, 0, 1])
		ax.set_xlabel('Density of states')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)
		ax.yaxis.set_label_position("right")
		ax.yaxis.tick_right()
		ax.legend()
	
	def DrawBandDos(self, type, N, U, is_unfold=0):
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(10, 6))

		if not is_unfold: e_min, e_max, fs = self.DrawBand(type, N, U, ax[0])
		else: e_min, e_max, fs = self.DrawBandUnfold(type, N, U, ax[0])

		self.DrawDos(type, N, U, ax[1], e_min, e_max)

		plt.suptitle(re.sub('/b', '\nb', fs))
		fig.tight_layout()

		if  not is_unfold: fig.savefig('%s' % (re.sub('output', 'diagram', re.sub('txt', 'png', fs))))
		else: fig.savefig('%s' % (re.sub('output', 'diagram', re.sub('band', 'banduf', re.sub('txt', 'png', fs)))))

		plt.show()

	def DrawBandCheck(self, type):
		fig, ax = plt.subplots(5, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(5, 10))

		U = 9
		for i, N in enumerate([2, 4, 6, 8, 10]):
			e_min, e_max, fs = self.DrawBand(type, N, U, ax[i][0])
			self.DrawDos(N, U, ax[i][1], e_min, e_max)
			ax[i][0].set_title('N=%d' % (N), loc='left')
			ax[i][1].legend().remove()

		plt.suptitle(self.dir + type)
		fig.tight_layout()
		plt.show()

	def DrawPhase(self, type):
		basis = self.type_info[type[0]]
		fs_list = [self.dir + fs for fs in os.listdir(self.dir) if re.match('band_%s' % type, fs)]
		N_list = []
		U_list = []

		for fs in fs_list:
			fs_dict = ReadFs(fs)
			N_list.append(fs_dict['N'])
			U_list.append(fs_dict['U'])

		N = sorted(list(set(N_list)))
		U = sorted(list(set(U_list)))
		m = np.zeros((len(U), len(N)))

		N_ins = []
		U_ins = []

		for fs in fs_list:
			fs_dict = ReadFs(fs)
			m[U.index(fs_dict['U'])][N.index(fs_dict['N'])] = abs(fs_dict['m'])

			if fs_dict['gap'] > self.tol:
				N_ins.append(fs_dict['N'])
				U_ins.append(fs_dict['U'])

		fig, ax = plt.subplots()

		#cs = ax.contour(N, U, Z, levels=10, cmap='binary_r')
		#plt.clabel(cs)
		plt.grid(True, alpha=0.5)
		cb = ax.contourf(N, U, m, levels=np.linspace(0, basis/2, 10), cmap='Blues_r')
		plt.colorbar(cb, label='Magnetization (M)', format='%.1f')
		ax.scatter(N_ins, U_ins, marker='+', color='black', label='Insulater')

		for n_ins in np.unique(N_ins):
			u_ins = np.array(U_ins)[np.where(np.abs(N_ins - n_ins) < 1e-6)]
			ax.plot([n_ins, n_ins], [np.min(u_ins), np.max(u_ins)], color='black')
			print(n_ins, u_ins)

		ax.set_xlabel('Occupation (N)')
		ax.set_ylabel('Interaction (U)')
		ax.legend(fontsize=8, loc='lower left')
		plt.suptitle(self.dir + type)	
		fig.tight_layout()
		fig.savefig('%s/phase_%s.png' % (re.sub('output', 'diagram', self.dir), type))
		plt.show()

	def DrawPhaseCheck(self, type1, type2):
		fs1_list = [self.dir + fs for fs in os.listdir(self.dir) if re.match('band_%s' % type1, fs)]
		fs2_list = [self.dir + fs for fs in os.listdir(self.dir) if re.match('band_%s' % type2, fs)]
		N_list = []
		U_list = []

		for fs in fs1_list:
			fs_dict = ReadFs(fs)
			N_list.append(fs_dict['N'])
			U_list.append(fs_dict['U'])

		N = sorted(list(set(N_list)))
		U = sorted(list(set(U_list)))
		m1 = np.zeros((len(U), len(N)))
		m2 = np.zeros((len(U), len(N)))
		e1 = np.zeros((len(U), len(N)))
		e2 = np.zeros((len(U), len(N)))

		for fs1, fs2 in zip(fs1_list, fs2_list):
			fs1_dict = ReadFs(fs1)
			fs2_dict = ReadFs(fs2)
			m1[U.index(fs1_dict['U'])][N.index(fs1_dict['N'])] = fs1_dict['m']
			m2[U.index(fs2_dict['U'])][N.index(fs2_dict['N'])] = fs2_dict['m']
			e1[U.index(fs1_dict['U'])][N.index(fs1_dict['N'])] = fs1_dict['e']
			e2[U.index(fs2_dict['U'])][N.index(fs2_dict['N'])] = fs2_dict['e']

		s1 = []
		s2 = []
		for i, n in enumerate(N):
			for j, u in enumerate(U):
				if (m1[j][i] > 0.1) and (m2[j][i] > 0.1):
					if e1[j][i] < e2[j][i]: s1.append([n, u])	
					else: s2.append([n, u])

		fig, ax = plt.subplots()

		plt.grid(True, alpha=0.5)
		ax.scatter(np.array(s1)[:, 0], np.array(s1)[:, 1], label=type1)
		ax.scatter(np.array(s2)[:, 0], np.array(s2)[:, 1], label=type2)
		ax.legend()
		plt.suptitle(self.dir + type1 + type2)	
		fig.tight_layout()
		fig.savefig('%s/phasec_%s%s.png' % (re.sub('output', 'diagram', self.dir), type1, type2))
		plt.show()

	def DrawSolution(self, type, N, U):
		N = float(N)
		U = float(U)

		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('sol_%s_N%.1f_U%.1f' % (type, N, U), f)][0]

		f = open(fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		if re.search('f', type):
			n_max = 2.1
			m_max = 1.1
			m_min = -1.1
		else:
			n_max = 4.2
			m_max = 2.2
			m_min = -2.2

		fig, ax = plt.subplots(1, 2, figsize=(10, 6))

		for i, n in enumerate([2, 4, 6]):
			ax[0].plot(data[:, 0], data[:, n],   lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors[i], label=self.labels[i])
			ax[1].plot(data[:, 0], data[:, n+1], lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors[i], label=self.labels[i])

		ax[0].grid(True)
		ax[0].set_xlabel('Iteration')
		ax[0].set_ylabel('Occupation per orbital')
		ax[0].set_ylim(-0.1, n_max)
		ax[0].legend()

		ax[1].grid(True)
		ax[1].set_xlabel('Iteration')
		ax[1].set_ylabel('Magnetization per orbital')
		ax[1].set_ylim(m_min, m_max)
		ax[1].legend()

		fs_s = fs.split(sep='_')
		fs_s = [s for s in fs_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		plt.suptitle(re.sub('/s', '\ns', fs))
		fig.tight_layout()
		fig.savefig('%s' % (re.sub('output', 'diagram', re.sub('txt', 'png', fs))))
		plt.show()
