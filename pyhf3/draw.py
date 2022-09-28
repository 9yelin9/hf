# pyhf3/draw.py : draw diagrams of 3 band models

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class Draw:
	def __init__(self, material, output, type, JU, SOC, basis, path, path_label):
		self.material = material
		self.output = output
		self.type = type
		self.JU = float(JU)
		self.SOC = float(SOC)
		self.bases = basis
		self.path = path
		self.path_label = path_label

		self.obt = 3
		self.basis = self.bases[0] if re.search('f', type) else self.bases[1]
		self.tol = 0.01

		self.dir = '%s/%s/JU%.2f_SOC%.2f/' % (material, output, self.JU, self.SOC)
		self.title = '%s J/U=%.2f SOC=%.2f\n' % (type, self.JU, self.SOC)
		self.colors=['tab:blue', 'tab:green', 'tab:red']
		self.labels=['xy', 'yz', 'zx']
		self.markers = ['s', 'o', '^']
		self.lss = ['-', '--', '-.']
	
	def DrawBandTest(self):
		fs = [f for f in os.listdir('%s/input/' % (self.material)) if re.search('band_%s.txt' % (self.type), f)][0]

		f = open('%s/input/%s' % (self.material, fs), 'r')
		data = np.genfromtxt(f)
		f.close()

		fig, ax = plt.subplots()

		for i in range(1, self.basis+1):
			ax.plot(data[:, 0], data[:, i], color='tab:blue')

		ax.grid(True)
		ax.set_xticks(self.path)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_title('%s/input/%s' % (self.material, fs))
		plt.show()
	
	def DrawBand(self, N, U, ax):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group()))

		f = open(self.dir + fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		e_min = np.min(data[:, 1:]) - fermi - 0.1
		e_max = np.max(data[:, 1:]) - fermi + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, self.basis+1):
			ax.plot(data[:, 0], data[:, i] - fermi, color='tab:blue')

		ax.grid(True)
		ax.set_xticks(self.path)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs 
	
	def DrawBandUF(self, N, U, ax):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group()))

		f_band = open(self.dir + fs, 'r')
		f_ufw  = open(self.dir + re.sub('band_', 'ufw_', fs), 'r')

		data_e = np.genfromtxt(f_band)
		data_w = np.genfromtxt(f_ufw)

		f_band.close()
		f_ufw.close()

		e_min = np.min(data_e[:, 1:]) - fermi - 0.1
		e_max = np.max(data_e[:, 1:]) - fermi + 0.1

		norm = plt.Normalize(0, 1)
		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, self.basis+1):
			ax.plot(data_e[:, 0], data_e[:, i] - fermi, linewidth=1, color='gray', alpha=0.5)
			ax.scatter(data_e[:, 0], data_e[:, i] - fermi, c=data_w[:, i], s=(10*(data_w[:, i]))**2, cmap='Blues', norm=norm)

		"""
		for i in range(1, self.basis+1):
			points = np.array([data_e[:, 0], data_e[:, i] - fermi]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)
			norm = plt.Normalize(0, 1)
			lc = LineCollection(segments, cmap='PuBu', norm=norm)
			lc.set_array(data_w[:, i])
			lc.set_linewidth(np.log(0.5+10*data_w[:, i]))
			line = ax.add_collection(lc)
		"""

		ax.grid(True)
		ax.set_xticks(self.path)
		ax.set_xticklabels(self.path_label)
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs

	def DrawDos(self, N, U, ax, e_min=None, e_max=None):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(self.dir) if re.search('dos_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group()))
		title = self.dir + fs

		f = open(self.dir + fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		if e_min == None: e_min = np.min(data[:, 0]) - fermi
		if e_max == None: e_max = np.max(data[:, 0]) - fermi

		dos_max = np.max(data[:, 1:])
		#dos_max = 1

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
	
	def DrawBandDos(self, N, U, is_unfold=0):
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(10, 6))

		if not is_unfold: e_min, e_max, fs = self.DrawBand(N, U, ax[0])
		else: e_min, e_max, fs = self.DrawBandUF(N, U, ax[0])

		self.DrawDos(N, U, ax[1], e_min, e_max)

		fs_s = fs.split(sep='_')
		fs_s = [s for s in fs_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		title = self.title
		for s in fs_s:
			s_name = re.search('[a-zA-Z]+', s).group()
			s_val = re.search('[-]?\d+[.]\d+', s).group()
			title += '%s=%s ' % (s_name, s_val)

		plt.suptitle(title)
		fig.tight_layout()

		if  not is_unfold: fig.savefig('%s/%s' % (re.sub('output', 'diagram', self.dir), re.sub('txt', 'png', fs)))
		else: fig.savefig('%s/%s' % (re.sub('output', 'diagram', self.dir), re.sub('band', 'banduf', re.sub('txt', 'png', fs))))

		plt.show()

	def DrawPhase(self):
		fs_list    = ['%s%s' % (self.dir, fs) for fs in os.listdir(self.dir) if re.match('band_%s' % (self.type), fs)]
		N_list     = [float(re.sub('_N', '', re.search('_N\d+[.]\d+', fs).group())) for fs in fs_list]
		U_list     = [float(re.sub('_U', '', re.search('_U\d+[.]\d+', fs).group())) for fs in fs_list]
		m_list     = [float(re.sub('_m', '', re.search('_m[-]?\d+[.]\d+', fs).group())) for fs in fs_list]
		fermi_list = [float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group())) for fs in fs_list]

		N = np.unique(N_list)
		U = np.unique(list(map(int, U_list)))
		m = np.zeros((len(U), len(N)))
		gap_tol = 5

		N_ins = []
		U_ins = []

		for i, fs in enumerate(fs_list):
			f = open(fs, 'r')
			data = np.genfromtxt(f)[:, 1:]
			f.close()

			if len(np.where(np.abs(U - U_list[i]) < 1e-6)[0]): 
				m[np.where(np.abs(U - U_list[i]) < 1e-6)[0][0]][np.where(np.abs(N - N_list[i]) < 1e-6)[0][0]] = m_list[i]

			is_cdt = len(data[np.where(np.abs(data - fermi_list[i]) < self.tol)])

			if not is_cdt:
				N_ins.append(N_list[i])
				U_ins.append(U_list[i])

			del data

		fig, ax = plt.subplots()

		#cs = ax.contour(N, U, Z, levels=10, cmap='binary_r')
		#plt.clabel(cs)
		plt.grid(True, alpha=0.5)
		cb = ax.contourf(N, U, m, cmap='Blues_r')
		plt.colorbar(cb, label='Magnetization (M)')
		ax.scatter(N_ins, U_ins, marker='+', color='black', label='Insulater')

		for n_ins in np.unique(N_ins):
			u_ins = np.array(U_ins)[np.where(np.abs(N_ins - n_ins) < 1e-6)]
			ax.plot([n_ins, n_ins], [np.min(u_ins), np.max(u_ins)], color='black')
			print(n_ins, u_ins)

		ax.set_xlabel('Occupation (N)')
		ax.set_ylabel('Interaction (U)')
		ax.legend(fontsize=8, loc='lower left')
		plt.suptitle(self.dir + self.type)	
		fig.tight_layout()
		fig.savefig('%s/phase_%s.png' % (re.sub('output', 'diagram', self.dir), self.type))
		plt.show()

	def DrawSol(self, N, U):
		N = float(N)
		U = float(U)

		fs = [f for f in os.listdir(self.dir) if re.search('sol_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]

		f = open(self.dir + fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		if re.search('f', self.type):
			n_max = 2.1
			m_max = 1.1
			m_min = -1.1
		else:
			n_max = 4.2
			m_max = 2.2
			m_min = -2.2

		fig, ax = plt.subplots(1, 2, figsize=(10, 6))

		for i in [3, 5, 7]:
			ax[0].plot(data[:, 0], data[:, i],   ls=self.lss[int(i/2)-1], marker=self.markers[int(i/2)-1], ms=abs(i-10), color=self.colors[int(i/2)-1], label=self.labels[int(i/2)-1])
			ax[1].plot(data[:, 0], data[:, i+1], ls=self.lss[int(i/2)-1], marker=self.markers[int(i/2)-1], ms=abs(i-10), color=self.colors[int(i/2)-1], label=self.labels[int(i/2)-1])

		ax[0].grid(True)
		ax[0].set_xlabel('Iteration')
		ax[0].set_ylabel('Occupation per orbital')
		ax[0].set_ylim(0, n_max)
		ax[0].legend()

		ax[1].grid(True)
		ax[1].set_xlabel('Iteration')
		ax[1].set_ylabel('Magnetization per orbital')
		ax[1].set_ylim(m_min, m_max)
		ax[1].legend()

		fs_s = fs.split(sep='_')
		fs_s = [s for s in fs_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		title = self.title
		for s in fs_s:
			s_name = re.search('[a-zA-Z]+', s).group()
			s_val = re.search('[-]?\d+[.]\d+', s).group()
			title += '%s=%s ' % (s_name, s_val)

		plt.suptitle(title)
		fig.tight_layout()
		fig.savefig('%s/%s' % (re.sub('output', 'diagram', self.dir), re.sub('band', 'banduf', re.sub('txt', 'png', fs))))
		plt.show()

	def DrawM(self, N):
		fs_list    = ['%s%s' % (self.dir, fs) for fs in os.listdir(self.dir) if re.match('band_%s_N%s' % (self.type, N), fs)]
		U_list     = [float(re.sub('_U', '', re.search('_U\d+[.]\d+', fs).group())) for fs in fs_list]
		m_list     = [float(re.sub('_m', '', re.search('_m[-]?\d+[.]\d+', fs).group())) for fs in fs_list]
		fermi_list = [float(re.sub('_fermi', '', re.search('_fermi[-]?\d+[.]\d+', fs).group())) for fs in fs_list]

		U = np.unique(list(map(int, U_list)))
		m = np.zeros(len(U))
		gap_tol = 5

		U_ins = []
		m_ins = []

		for i, fs in enumerate(fs_list):
			f = open(fs, 'r')
			data = np.genfromtxt(f)[:, 1:]
			f.close()

			is_cdt = len(data[np.where(np.abs(data - fermi_list[i]) < self.tol)])

			if not is_cdt:
				U_ins.append(U_list[i])
				m_ins.append(m_list[i])

			del data

		fig, ax = plt.subplots()

		plt.grid(True, alpha=0.5)
		ax.plot(U_list, m_list, '.-', ms=10, zorder=1)
		ax.scatter(U_ins, m_ins, color='tab:orange', label='Insulater', zorder=2)

		ax.set_xlabel('Interaction (U)')
		ax.set_ylabel('Magnetization (m)')
		ax.legend(fontsize=8, loc='upper left')
		plt.suptitle(self.dir + self.type + '_N%s' % (N))	
		fig.tight_layout()
		fig.savefig('%s/m_%s_N%s.png' % (re.sub('output', 'diagram', self.dir), self.type, N))
		plt.show()
