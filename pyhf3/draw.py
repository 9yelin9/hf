# pyhf3/draw.py : draw diagrams of 3 band models

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from .read import ReadFs, MakeGroundIdx

class Draw:
	def __init__(self, input_path, output_path, info_path, info_cell, JU, SOC, tol):
		self.input_path = input_path
		self.path_point = [path_point for path_point, path_label in info_path]
		self.path_label = [path_label for path_point, path_label in info_path]
		self.info_cell = info_cell
		self.JU = float(JU)
		self.SOC = float(SOC)
		self.tol = tol

		self.obt = 3

		self.dir = '%s/JU%.2f_SOC%.2f/' % (output_path, self.JU, self.SOC)
		self.colors=['tab:blue', 'tab:green', 'tab:red']
		self.labels=[r'$d_{xy}$', r'$d_{yz}$', r'$d_{zx}$']
		self.markers = ['s', 'o', '^']
		self.lss = ['-', '--', '-.']

		os.makedirs(re.sub('output', 'diagram', self.dir), exist_ok=True)
		plt.rcParams.update({'font.size': 25})
	
	def DrawBandTB(self, type):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		fs = [f for f in os.listdir(self.input_path) if re.search('band_%s.txt' % (type), f)][0]

		f = open(self.input_path + fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		x = np.arange(len(data))

		fig, ax = plt.subplots()

		for i in range(Nb):
			ax.plot(x, data[:, i], color='tab:blue')

		ax.grid(True, axis='x')
		ax.set_xticks(self.path_point)
		ax.set_xticklabels(self.path_label)
		ax.set_ylabel(r'$E$')
		ax.set_title(r'$%s$-type' % type[0], loc='left', fontdict={'fontsize':'medium'})
		plt.show()
	
	def DrawBand(self, type, N, U, ax, is_unfold=0):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fs_dict = ReadFs(fs)
		fermi = fs_dict['fermi']

		f = open(fs, 'r')
		data = np.genfromtxt(f)
		data = data - fermi
		f.close()

		x = np.arange(len(data))
		e_min = np.min(data[:, 1:]) - 0.1
		e_max = np.max(data[:, 1:]) + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		if is_unfold:
			f = open(re.sub('band_', 'ufw_', fs), 'r')
			data_w = np.genfromtxt(f)
			f.close()

			norm = plt.Normalize(0, 1)

			for i in range(Nb):
				ax.plot(x, data[:, i], linewidth=1, color='gray', alpha=0.5)
				sc = ax.scatter(x, data[:, i], c=data_w[:, i], s=(10*(data_w[:, i]))**2, cmap='Blues', norm=norm)
		else:
			sc = 0
			for i in range(Nb):
				ax.plot(x, data[:, i], color='tab:blue')

		ax.grid(True, axis='x')
		ax.set_xticks(self.path_point)
		ax.set_xticklabels(self.path_label)
		ax.set_ylabel(r'$E - E_{F}$')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fs, sc
	
	def DrawDOS(self, type, N, U, ax, e_min=None, e_max=None):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
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

		data[:, 1:] = data[:, 1:] / np.pi
		data[:, 1:] = data[:, 1:] / (2*np.pi)**3
		dos_max = np.max(data[:, 1:])

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.obt+1): 
			ax.plot(data[:, i],           data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors[i-1], label=self.labels[i-1])
			ax.plot(-data[:, i+self.obt], data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors[i-1])

		ax.grid(True, axis='x')
		ax.set_xticks([-1, 0, 1])
		ax.set_xticklabels([1, 0, 1])
		ax.set_xlabel('PDOS')
		ax.set_ylim(e_min, e_max)
		ax.yaxis.tick_right()
		ax.yaxis.set_ticklabels([])
		ax.legend(bbox_to_anchor=(1.05, 1.03), frameon=False, labelspacing=0.1, fontsize=22)
	
	def DrawBandDOS(self, type, N, U, is_unfold=0):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(11, 6))

		e_min, e_max, fs, sc = self.DrawBand(type, N, U, ax[0], is_unfold=is_unfold)
		self.DrawDOS(type, N, U, ax[1], e_min, e_max)

		fname = re.sub('output', 'diagram', re.sub('txt', 'png', fs))

		if sc:
			cb = plt.colorbar(sc, shrink=0.6, pad=0.155, anchor=(0.00, 0.03), format='%.1f')
			cb.ax.tick_params(labelsize=20)
			fname = re.sub('output', 'diagram', re.sub('band', 'unfold', re.sub('txt', 'png', fs)))

		ax[0].set_title(r'$%s$-type $N = %.1f$ $U = %.1f$ $J/U = %.1f$' % (type[0], N, U, self.JU), loc='left', fontdict={'fontsize':'medium'})
		fig.tight_layout()
		plt.subplots_adjust(wspace=0.03)
		fig.savefig('%s' % fname)

		print(fs)
		plt.show()

	def DrawPhase(self, type):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		fs_list = [self.dir + fs for fs in os.listdir(self.dir) if re.match('band_%s' % type, fs)]
		N_list = []
		U_list = []

		idx_list = MakeGroundIdx(fs_list)

		for i in idx_list:
			fs = fs_list[i]
			fs_dict = ReadFs(fs)
			N_list.append(fs_dict['N'])
			U_list.append(fs_dict['U'])

		N = sorted(list(set(N_list)))
		U = sorted(list(set(U_list)))
		m = np.zeros((len(U), len(N)))

		type_ins = []
		N_ins = []
		U_ins = []

		for i in idx_list:
			fs = fs_list[i]
			fs_dict = ReadFs(fs)
			m[U.index(fs_dict['U'])][N.index(fs_dict['N'])] = abs(fs_dict['m'])

			if fs_dict['gap'] > self.tol:
				type_ins.append(fs_dict['type'])
				N_ins.append(fs_dict['N'])
				U_ins.append(fs_dict['U'])

		fig, ax = plt.subplots(figsize=(10, 6))

		cf = ax.contourf(N, U, m, levels=np.linspace(0, Nb/2, 101), cmap='Blues_r')
		cb = plt.colorbar(cf, shrink=0.85, anchor=(0.00, 0.03), format='%.1f')
		cb.set_ticks(np.arange(0, Nb/2+0.1, 1))
		cb.set_ticklabels(np.arange(0, Nb/2+0.1, 1))
		cb.set_label(r'$m$', rotation=0, labelpad=20)
		cb.ax.tick_params(labelsize=18)

		labels = ['_nolegend_' for _ in range(Nb)]
		labels[0] = 'Insulator'

		for i, n_ins in enumerate(np.unique(N_ins)):
			idx = np.where(np.abs(N_ins - n_ins) < 1e-6)
			t_ins = np.array(type_ins)[idx]
			u_ins = np.array(U_ins)[idx]
			ax.plot([n_ins, n_ins], [np.min(u_ins), np.max(u_ins)], lw=7, color='black', label=labels[i])
			print(n_ins, list(zip(list(t_ins), list(u_ins))))
		ax.plot([Nb], [max(U)], alpha=1) 

		ax.set_xticks(np.arange(0, Nb+1, Ni))
		ax.set_xticklabels(range(0, Nb+1, Ni))
		ax.set_yticks(np.arange(0, max(U)+1, 1))
		ax.set_ylim(min(U), max(U))
		ax.set_xlabel(r'$N$')
		ax.set_ylabel(r'$U$', rotation=0, labelpad=20)
		ax.legend(bbox_to_anchor=(1.03, 1.03), frameon=False, fontsize=17)
		ax.set_title(r'$%s$-type $J/U = %.1f$' % (type[0], self.JU), loc='left', fontdict={'fontsize':'medium'})
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

		plt.grid(True, alpha=0.5, axis='x')
		ax.scatter(np.array(s1)[:, 0], np.array(s1)[:, 1], label=type1)
		ax.scatter(np.array(s2)[:, 0], np.array(s2)[:, 1], label=type2)
		ax.legend()
		ax.set_title(r'$%s$-type $J/U = %.1f$' % (type[0], self.JU), loc='left', fontdict={'fontsize':'medium'})
		fig.tight_layout()
		fig.savefig('%s/phasec_%s%s.png' % (re.sub('output', 'diagram', self.dir), type1, type2))
		plt.show()

	def DrawSolution(self, type, N, U):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fs = [self.dir + f for f in os.listdir(self.dir) if re.search('sol_%s_N%.1f_U%.1f' % (type, N, U), f)][0]

		f = open(fs, 'r')
		data = np.genfromtxt(f)
		f.close()

		n_max = 4.2
		m_max = 2.2
		m_min = -2.2

		fig, ax = plt.subplots(1, 2, figsize=(10, 6))

		for i, n in enumerate([2, 4, 6]):
			ax[0].plot(data[:, 0], data[:, n],   lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors[i], label=self.labels[i])
			ax[1].plot(data[:, 0], data[:, n+1], lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors[i], label=self.labels[i])

		ax[0].grid(True)
		ax[0].set_xlabel('Iteration')
		ax[0].set_ylabel('Occupation')
		ax[0].set_ylim(-0.1, n_max)
		ax[0].legend()

		ax[1].grid(True)
		ax[1].set_xlabel('Iteration')
		ax[1].set_ylabel('Magnetization')
		ax[1].set_ylim(m_min, m_max)
		ax[1].legend()

		fs_s = fs.split(sep='_')
		fs_s = [s for s in fs_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		ax[0].set_title(r'$%s$-type $N = %d$ $U = %.1f$ $J/U = %.1f$' % (type[0], N , U, self.JU), loc='left', fontdict={'fontsize':'medium'})
		fig.tight_layout()
		fig.savefig('%s' % (re.sub('output', 'diagram', re.sub('txt', 'png', fs))))
		plt.show()
