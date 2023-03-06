# pyhf3/draw.py : draw diagrams of 3 band models

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .read import ReadInfo, ReadFn, GenGroundIdx

class Draw:
	def __init__(self, name, JU, SOC):
		self.path_input, self.path_output, self.info_path, self.info_cell = ReadInfo(name)
		self.points = [point for point, _ in self.info_path]
		self.labels = [label for _, label in self.info_path]
		self.labels_g = [label.replace('G', r'$\Gamma$') for label in self.labels]
		self.JU = float(JU)
		self.SOC = float(SOC)

		self.obt = 3

		self.dir = '%s/JU%.2f_SOC%.2f/' % (self.path_output, self.JU, self.SOC)
		self.colors_ob=['tab:blue', 'tab:green', 'tab:red']
		self.labels_ob=[r'$d_{xy}$', r'$d_{yz}$', r'$d_{zx}$']
		self.markers = ['s', 'o', '^']
		self.lss = ['-', '--', '-.']

		os.makedirs(re.sub('output', 'diagram', self.dir), exist_ok=True)
		plt.rcParams.update({'font.size': 35})
		plt.rcParams.update({'font.family': 'sans-serif'})
		plt.rcParams.update({'font.serif': 'Helvetica Neue'})
		#plt.rcParams.update({'mathtext.fontset': 'cm'})
	
	def DrawBandTB(self, type):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		fn = [f for f in os.listdir(self.path_input) if re.search('band_%s.txt' % (type), f)][0]

		f = open(self.path_input +'/'+ fn, 'r')
		data = np.genfromtxt(f)
		f.close()

		x = np.arange(len(data))

		fig, ax = plt.subplots(dpi=600)

		for i in range(Nb):
			ax.plot(x, data[:, i], color='tab:blue')

		ax.grid(True, axis='x')
		ax.set_xticks(self.points)
		ax.set_xticklabels(self.labels_g)
		ax.set_ylabel(r'$E$')
		ax.set_title(r'$%s$-type' % type[0], loc='left', fontdict={'fontsize':'medium'})
		plt.show()
	
	def DrawBand(self, type, N, U, ax, is_unfold=0):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fn = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fn_dict = ReadFn(fn)
		fermi = fn_dict['fermi']

		f = open(fn, 'r')
		data = np.genfromtxt(f)
		#fermi = (np.max(data[:, int(N)-1]) + np.min(data[:, int(N)])) / 2
		data = data - fermi
		f.close()

		x = np.arange(len(data))
		e_min = np.min(data) - 0.5
		e_max = np.max(data) + 0.5

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		if is_unfold:
			f = open(re.sub('band_', 'ufw_', fn), 'r')
			data_w = np.genfromtxt(f)
			f.close()

			itv = 10
			norm = plt.Normalize(-0.5, 1)

			data_tot = np.array([[i, e, w] for i in x[::itv] for e, w in zip(data[i], data_w[i])])
			data_tot = data_tot[data_tot[:, -1].argsort()[::-1]] 
			
			#sc = ax.scatter(data_tot[:, 0], data_tot[:, 1], lw=0.1, c=data_tot[:, 2], s=(10*(data_tot[:, 2]))**(2.5), cmap='Blues', norm=norm)
			sc = ax.scatter(data_tot[:, 0], data_tot[:, 1], marker='$○$', lw=0.1, c=data_tot[:, 2], s=4+(10*(data_tot[:, 2]))**(2.5), cmap='Blues', norm=norm)

			#axins = ax.inset_axes([0.92, 0.32, 0.052, 0.32])
			axins = ax.inset_axes([0.92, 0.35, 0.052, 0.37])
			#axins = ax.inset_axes([0.92, 0.32, 0.065, 0.35])
			w = np.array([(1.35**i/1.35**9) for i in range(10)])
			axins.scatter(np.zeros(10), w, marker='$○$', lw=0.1, c=w, s=4+(10*(w))**(2.5), cmap='Blues', norm=norm)
			axins.set_xticks([])
			axins.set_yticks([np.min(w), 1])
			axins.set_yticklabels([0, 1])
			#axins.yaxis.tick_right()
			axins.set_ylim([-0.05, 1.2])
			axins.tick_params(axis='both', which='major', labelsize=18)
		else:
			sc = 0
			for i in range(Nb):
				ax.plot(x, data[:, i], color='tab:blue')
		
		ax.grid(True, axis='x')
		ax.set_xticks(self.points)
		ax.set_xticklabels(self.labels_g)
		#ax.set_ylabel(r'$\omega$', rotation=0)
		ax.set_ylabel(r'$E-E_{F}$')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fn, sc
	
	def DrawDOS(self, type, N, U, ax, e_min=None, e_max=None, eta=0.10):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fn = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fn_dict = ReadFn(fn)

		f = open(re.sub('band', 'dos', re.sub('%s_'%type, '%s_eta%.2f_'%(type, eta), fn)), 'r')
		data = np.genfromtxt(f)
		f.close()

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.obt+1): 
			ax.plot(data[:, i]+data[:, i+self.obt], data[:, 0], ls=self.lss[i-1], lw=abs(i-5), color=self.colors_ob[i-1], label=self.labels_ob[i-1])
			#ax.fill_betweenx(data[:, 0], data[:, i]+data[:, i+self.obt], color=self.colors_ob[i-1], alpha=0.5)

		#ax.grid(True, axis='x')
		ax.set_xlim([0, 1.1])
		ax.set_xticks([0, 1])
		ax.set_ylim(e_min, e_max)
		ax.yaxis.tick_right()
		ax.yaxis.set_ticklabels([])
		#ax.legend(bbox_to_anchor=(1.05, 1.03), frameon=False, labelspacing=0.1, fontsize=22)
		ax.legend(fontsize=30, labelspacing=0.02, handletextpad=0.3, handlelength=1.0, borderpad=0.1, borderaxespad=0.1, frameon=False, loc='lower right')

	def DrawDOS_(self, type, N, U, ax, e_min=None, e_max=None, eta=0.10):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fn = [self.dir + f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (type, N, U), f)][0]
		fn_dict = ReadFn(fn)
		#fermi = fn_dict['fermi']

		f = open(re.sub('band', 'dos', re.sub('%s_'%type, '%s_eta%.2f_'%(type, eta), fn)), 'r')
		data = np.genfromtxt(f)
		f.close()

		#if e_min == None: e_min = np.min(data[:, 0]) - fermi
		#if e_max == None: e_max = np.max(data[:, 0]) - fermi
		if e_min == None: e_min = np.min(data[:, 0])
		if e_max == None: e_max = np.max(data[:, 0])

		dos_max = np.max(data[:, 1:])

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.obt+1): 
			#ax.plot(data[:, i],           data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors_ob[i-1], label=self.labels_ob[i-1])
			#ax.plot(-data[:, i+self.obt], data[:, 0] - fermi, ls=self.lss[i-1], lw=abs(i-5), color=self.colors_ob[i-1])
			ax.plot( data[:, i],          data[:, 0], ls=self.lss[i-1], lw=abs(i-5), color=self.colors_ob[i-1], label=self.labels_ob[i-1])
			ax.plot(-data[:, i+self.obt], data[:, 0], ls=self.lss[i-1], lw=abs(i-5), color=self.colors_ob[i-1])
			#ax.fill_betweenx(data[:, 0],  data[:, i],          color=self.colors_ob[i-1], alpha=0.5)
			#ax.fill_betweenx(data[:, 0], -data[:, i+self.obt], color=self.colors_ob[i-1], alpha=0.5)

		#ax.text( 0.8, e_max-0.5, r'$\uparrow$',   ha='center')
		#ax.text(-0.8, e_max-0.5, r'$\downarrow$', ha='center')

		ax.grid(True, axis='x')
		#ax.set_xticks([-1, 0, 1])
		#ax.set_xticklabels([1, 0, 1])
		ax.set_xlim([-0.5, 0.5])
		ax.set_xticks([-0.4, 0, 0.4])
		ax.set_xticklabels([0.4, 0, 0.4])
		#ax.set_xlabel('PDOS')
		ax.set_ylim(e_min, e_max)
		ax.yaxis.tick_right()
		ax.yaxis.set_ticklabels([])
		#ax.legend(bbox_to_anchor=(1.05, 1.03), frameon=False, labelspacing=0.1, fontsize=22)
		ax.legend(fontsize=30, labelspacing=0.02, handletextpad=0.3, handlelength=1.0, borderpad=0.02, borderaxespad=0.1, frameon=False, loc='lower right')
	
	def DrawBandDOS(self, type, N, U, is_unfold=0):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, width_ratios=[3, 1], figsize=(10, 5), dpi=600, constrained_layout=True)

		e_min, e_max, fn, sc = self.DrawBand(type, N, U, ax[0], is_unfold=is_unfold)
		self.DrawDOS(type, N, U, ax[1], e_min, e_max)

		fname = re.sub('output', 'diagram', re.sub('txt', 'eps', fn))

		#ax[0].set_title(r'%s-type $J/U = %.1f$ $N = %.1f$ $U = %.1f$' % (type[0], self.JU, N, U), loc='left', fontdict={'fontsize':'medium'})
		#fig.tight_layout()
		#plt.subplots_adjust(wspace=0.03)
		fig.savefig('%s' % fname)

		print(fn)
		plt.show()

	def DrawPhase(self, type, tol_gap=0.09, tol_m=0.1, ax=0, show_xticks=True, show_yticks=True):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		tol_gap = float(tol_gap)
		tol_m = float(tol_m)

		fn_list = [self.dir + fn for fn in os.listdir(self.dir) if re.match('band_%s' % type, fn)]
		idx_list = GenGroundIdx(fn_list)

		mag = []
		ins = []
		
		n_list = np.arange(0.2, 12, 0.2)
		u_list = np.arange(0, 8.1, 0.1)

		for u in u_list:
			mag.append([0, u, 0])

		for n in n_list:
			for u in u_list:
				fn = [fn_list[i] for i in idx_list if re.search('N%.1f_U%.1f'%(n, u), fn_list[i])][0]
				fn_dict = ReadFn(fn)
				mag.append([fn_dict['N'], fn_dict['U'], abs(fn_dict['m'])])
				if fn_dict['gap'] > tol_gap: ins.append([fn_dict['N'], fn_dict['U'], [fn_dict['type']]])

		for u in u_list:
			mag.append([12, u, 0])

		mag = np.array(mag)
		ins = pd.DataFrame(ins, columns=['N', 'U', 'type'])

		X = np.reshape(mag[:, 0], (int(Nb/0.2)+1, 81))
		Y = np.reshape(mag[:, 1], (int(Nb/0.2)+1, 81))
		Z = np.reshape(mag[:, 2], (int(Nb/0.2)+1, 81))

		if not ax: fig, ax = plt.subplots(figsize=(10, 5), dpi=600)

		#ct = ax.contour(N, U, m, levels=[tol_m], colors='w', linestyles='dotted')
		#if abs(tol_m - 0.1) < 1e-6: ax.clabel(ct, ct.levels, inline=True, fmt='%.1f', fontsize=16)

		cf = ax.contourf(X, Y, Z, levels=np.linspace(0, Nb/2, 10), cmap='Blues_r')
		#cb = plt.colorbar(cf, format='%.1f')
		#cb.set_ticks([0, Nb/2])
		#cb.set_ticklabels(['0', '3'])

		#cb = plt.colorbar(cf, shrink=0.85, anchor=(0.00, 0.03), format='%.1f')
		#cb.set_ticks(np.arange(0, Nb/2+0.1, 2))
		#cb.set_label(r'$m$', rotation=0, labelpad=20)
		#cb.ax.tick_params(labelsize=18)

		labels = ['_nolegend_' for _ in range(Nb)]
		labels[0] = 'Insulator'

		for i, ins_n in enumerate(np.unique(ins['N'])):
			idcs = np.where(np.abs(ins['N'] - ins_n) < 1e-6)
			ins_u = np.array(ins['U'])[idcs]
			ins_t = np.array(ins['type'])[idcs]
			ax.plot([ins_n, ins_n], [np.min(ins_u), np.max(ins_u)], lw=7, color='black', label=labels[i])
			#print(ins_n, list(zip(list(ins_t), list(ins_u))))
		ax.plot([np.max(X)], [np.max(Y)], alpha=1) 

		ax.text(0.5, 0.75, type[0].upper(), bbox={'boxstyle':'Square', 'facecolor':'white'})
		ax.set_ylim(0, 8)
		ax.set_xticks(np.arange(0, Nb+1, Ni))
		ax.set_yticks(range(0, 9, 2))
		ax.set_xticklabels(range(7))
		ax.set_yticklabels(range(0, 9, 2))
		ax.set_xlabel(r'$N$')
		ax.set_ylabel(r'$U$', labelpad=20)

		if not show_xticks:
			ax.set_xticklabels([])
			ax.set_xlabel('')
		if not show_yticks:
			ax.set_yticklabels([])
			ax.set_ylabel('')

		#ax.set_yticks(np.arange(0, max(U)+1, 1))
		#ax.legend(bbox_to_anchor=(1.03, 1.03), frameon=False, fontsize=17)
		#ax.set_title(r'%s-type $J/U = %.1f$' % (type[0], self.JU), loc='left', fontdict={'fontsize':'medium'})
		#ax.set_title('%s-type' % type[0], pad=20, fontdict={'fontsize':'medium'})

		if not ax:
			fig.tight_layout()
			fig.savefig('%s/phase_%s_%.3f.png' % (re.sub('output', 'diagram', self.dir), type, tol_m))
			plt.show()

		return cf

	def DrawPhaseCheck(self, type1, type2):
		fn1_list = [self.dir + fn for fn in os.listdir(self.dir) if re.match('band_%s' % type1, fn)]
		fn2_list = [self.dir + fn for fn in os.listdir(self.dir) if re.match('band_%s' % type2, fn)]
		N_list = []
		U_list = []

		for fn in fn1_list:
			fn_dict = ReadFn(fn)
			N_list.append(fn_dict['N'])
			U_list.append(fn_dict['U'])

		N = sorted(list(set(N_list)))
		U = sorted(list(set(U_list)))
		m1 = np.zeros((len(U), len(N)))
		m2 = np.zeros((len(U), len(N)))
		e1 = np.zeros((len(U), len(N)))
		e2 = np.zeros((len(U), len(N)))

		for fn1, fn2 in zip(fn1_list, fn2_list):
			fn1_dict = ReadFn(fn1)
			fn2_dict = ReadFn(fn2)
			m1[U.index(fn1_dict['U'])][N.index(fn1_dict['N'])] = fn1_dict['m']
			m2[U.index(fn2_dict['U'])][N.index(fn2_dict['N'])] = fn2_dict['m']
			e1[U.index(fn1_dict['U'])][N.index(fn1_dict['N'])] = fn1_dict['e']
			e2[U.index(fn2_dict['U'])][N.index(fn2_dict['N'])] = fn2_dict['e']

		s1 = []
		s2 = []
		for i, n in enumerate(N):
			for j, u in enumerate(U):
				if (m1[j][i] > 0.1) and (m2[j][i] > 0.1):
					if e1[j][i] < e2[j][i]: s1.append([n, u])	
					else: s2.append([n, u])

		fig, ax = plt.subplots(dpi=600)

		plt.grid(True, alpha=0.5, axis='x')
		ax.scatter(np.array(s1)[:, 0], np.array(s1)[:, 1], label=type1)
		ax.scatter(np.array(s2)[:, 0], np.array(s2)[:, 1], label=type2)
		ax.legend()
		ax.set_title(r'$%s$-type $J/U = %.1f$' % (type[0], self.JU), loc='left', fontdict={'fontsize':'medium'})
		fig.tight_layout()
		fig.savefig('%s/phasec_%s%s.pdf' % (re.sub('output', 'diagram', self.dir), type1, type2))
		plt.show()

	def DrawSolution(self, type, N, U):
		Ni = self.info_cell[type[0]][0]
		Nc = self.info_cell[type[0]][1]
		Nb = Ni * Nc * 2
		N = float(N)
		U = float(U)

		fn = [self.dir + f for f in os.listdir(self.dir) if re.search('sol_%s_N%.1f_U%.1f' % (type, N, U), f)][0]

		f = open(fn, 'r')
		data = np.genfromtxt(f)
		f.close()

		n_max = 4.2
		m_max = 2.2
		m_min = -2.2

		fig, ax = plt.subplots(1, 2, figsize=(10, 5), dpi=600)

		for i, n in enumerate([2, 4, 6]):
			ax[0].plot(data[:, 0], data[:, n],   lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors_ob[i], label=self.labels_ob[i])
			ax[1].plot(data[:, 0], data[:, n+1], lw=abs(i-4), ls=self.lss[i], marker=self.markers[i], ms=abs(n-10), color=self.colors_ob[i], label=self.labels_ob[i])

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

		fn_s = fn.split(sep='_')
		fn_s = [s for s in fn_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		ax[0].set_title(r'$%s$-type $N = %d$ $U = %.1f$ $J/U = %.1f$' % (type[0], N , U, self.JU), loc='left', fontdict={'fontsize':'medium'})
		fig.tight_layout()
		fig.savefig('%s' % (re.sub('output', 'diagram', re.sub('txt', 'pdf', fn))))
		plt.show()
