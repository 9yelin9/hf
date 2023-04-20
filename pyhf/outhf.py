# pyhf/outhf.py : class to generate output

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .mod import FnDict, GroundOnly

class OutHF:
	def __init__(self, save, type, JU, SOC):
		self.Nc = 3
		self.Nb = 6 if re.search('F', type) else 12

		self.save = save
		self.type = type
		self.JU   = float(JU)
		self.SOC  = float(SOC)

		self.path_output = 'output/%s/%s_JU%.2f_SOC%.2f/' % (save, self.type, self.JU, self.SOC)
		self.path_save = '%s/diagram/' % self.path_output
		os.makedirs(self.path_save, exist_ok=True)

		self.hsp_point = [0, 198, 396, 677, 1023]
		self.hsp_label = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']

		self.t2g_color = ['tab:blue', 'tab:green', 'tab:red']
		self.t2g_label = [r'$d_{xy}$', r'$d_{yz}$', r'$d_{zx}$']
		self.t2g_mk    = ['s', 'o', '^']
		self.t2g_ls    = ['-', '--', '-.']

		plt.rcParams.update({'font.size': 35})
		plt.rcParams.update({'font.family': 'sans-serif'})
		plt.rcParams.update({'font.serif': 'Helvetica Neue'})
	
	def ShowBand(self, N, U, ax, is_unfold):
		N = float(N)
		U = float(U)

		fn = [self.path_output+'/band/'+f for f in os.listdir(self.path_output+'/band') if re.search('N%.1f_U%.1f' % (N, U), f)][0]
		with open(fn, 'r') as f: data = np.genfromtxt(f, skip_header=1)

		x = np.arange(data.shape[0])
		e_min = np.min(data) - 0.5
		e_max = np.max(data) + 0.5

		ax.axhline(y=0, ls=':', lw=2, color='dimgrey')

		if is_unfold:
			axins = ax.inset_axes([0.92, 0.35, 0.052, 0.37])

			itv = 10
			norm = plt.Normalize(-0.5, 1)
			w = np.array([(1.35**i/1.35**9) for i in range(10)])

			data_tot = np.array([[i, e, w] for i in x[::itv] for e, w in zip(data[i, :self.Nb], data[i, self.Nb:])])
			data_tot = data_tot[data_tot[:, -1].argsort()[::-1]] 
			
			ax.scatter(data_tot[:, 0], data_tot[:, 1], marker='$○$', lw=0.1, c=data_tot[:, 2], s=4+(10*(data_tot[:, 2]))**(2.5), cmap='Blues', norm=norm)
			axins.scatter(np.zeros(10), w, marker='$○$', lw=0.1, c=w, s=4+(10*(w))**(2.5), cmap='Blues', norm=norm)

			axins.set_xticks([])
			axins.set_yticks([np.min(w), 1])
			axins.set_yticklabels([0, 1])
			#axins.yaxis.tick_right()
			axins.set_ylim([-0.05, 1.2])
			axins.tick_params(axis='both', which='major', labelsize=18)
		else:
			for i in range(self.Nb): ax.plot(x, data[:, i], color='tab:blue')
		
		ax.grid(True, axis='x')
		ax.set_xticks(self.hsp_point)
		ax.set_xticklabels(self.hsp_label)
		ax.set_ylabel(r'$E-E_{F}$')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fn
	
	def ShowDOS(self, N, U, ax, e_min, e_max, ep):
		N = float(N)
		U = float(U)

		fn = [self.path_output+'/dos_ep%.2f/' % ep+f for f in os.listdir(self.path_output+'/dos_ep%.2f' % ep) if re.search('N%.1f_U%.1f' % (N, U), f)][0]
		with open(fn, 'r') as f: data = np.genfromtxt(f, skip_header=1)

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')
		
		for i in range(1, self.Nc+1): 
			ax.plot(data[:, i]+data[:, i+self.Nc], data[:, 0], ls=self.t2g_ls[i-1], lw=abs(i-5), color=self.t2g_color[i-1], label=self.t2g_label[i-1])
			#ax.fill_betweenx(data[:, 0], data[:, i]+data[:, i+self.Nc], color=self.t2g_color[i-1], alpha=0.5)

		#ax.grid(True, axis='x')
		ax.set_xlim([0, 1.1])
		ax.set_xticks([0, 1])
		ax.set_ylim(e_min, e_max)
		ax.yaxis.tick_right()
		ax.yaxis.set_ticklabels([])
		#ax.legend(bbox_to_anchor=(1.05, 1.03), frameon=False, labelspacing=0.1, fontsize=22)
		ax.legend(fontsize=30, labelspacing=0.02, handletextpad=0.3, handlelength=1.0, borderpad=0.1, borderaxespad=0.1, frameon=False, loc='lower right')

	def ShowBandDOS(self, N, U, ep=0.02, is_unfold=0):
		N  = float(N)
		U  = float(U)
		ep = float(ep)

		fig, ax = plt.subplots(1, 2, width_ratios=[3, 1], figsize=(10, 5), constrained_layout=True)

		e_min, e_max, fn = self.ShowBand(N, U, ax[0], is_unfold=int(is_unfold))
		self.ShowDOS(N, U, ax[1], e_min, e_max, ep)

		print(fn)
		fname = re.sub('output/', '', re.sub('txt', 'png', fn))
		#fig.savefig('%s' % fname)
		plt.show()

	def ShowPhase(self):
		tol_gap = 0.1
		tol_m   = 0.1

		save_list = ['output/%s/%s' % (self.save, s) for s in os.listdir('output/%s' % self.save) if re.search('%s\d_' % self.type[0], s)]
		fn_list = ['%s/band/%s' % (s, f) for s in save_list for f in os.listdir('%s/band' % s)]
		grd_idx = GroundOnly(fn_list)

		n_list = np.arange(0.2, 12, 0.2)
		u_list = np.arange(0, 8.1, 1)

		mag = []
		ins = []
		for u in u_list:
			mag.append([0, u, 0])
			for n in n_list:
				fn = [fn_list[i] for i in grd_idx if re.search('N%.1f_U%.1f'%(n, u), fn_list[i])][0]
				mag.append([FnDict(fn)['N'], FnDict(fn)['U'], abs(FnDict(fn)['m'])])
				if FnDict(fn)['gap'] > tol_gap: ins.append([FnDict(fn)['N'], FnDict(fn)['U'], [FnDict(fn)['type']]])
			mag.append([12, u, 0])

		mag = np.array(mag)
		ins = pd.DataFrame(ins, columns=['N', 'U', 'type'])

		X = np.reshape(mag[:, 0], (9, int(self.Nb/0.2)+1))
		Y = np.reshape(mag[:, 1], (9, int(self.Nb/0.2)+1))
		Z = np.reshape(mag[:, 2], (9, int(self.Nb/0.2)+1))

		fig, ax = plt.subplots(figsize=(10, 5))

		#ct = ax.contour(N, U, m, levels=[tol_m], colors='w', linestyles='dotted')
		#if abs(tol_m - 0.1) < 1e-6: ax.clabel(ct, ct.levels, inline=True, fmt='%.1f', fontsize=16)

		cf = ax.contourf(X, Y, Z, levels=np.linspace(0, self.Nb/2, 10), cmap='Blues_r')
		cb = plt.colorbar(cf, format='%.1f')
		cb.set_ticks([0, self.Nb/2])
		cb.set_ticklabels(['0', '3'])

		#cb = plt.colorbar(cf, shrink=0.85, anchor=(0.00, 0.03), format='%.1f')
		#cb.set_ticks(np.arange(0, self.Nb/2+0.1, 2))
		#cb.set_label(r'$m$', rotation=0, labelpad=20)
		#cb.ax.tick_params(labelsize=18)

		labels = ['_nolegend_' for _ in range(self.Nb)]
		labels[0] = 'Insulator'

		for i, ins_n in enumerate(np.unique(ins['N'])):
			idcs = np.where(np.abs(ins['N'] - ins_n) < 1e-6)
			ins_u = np.array(ins['U'])[idcs]
			ins_t = np.array(ins['type'])[idcs]
			ax.plot([ins_n, ins_n], [np.min(ins_u), np.max(ins_u)], lw=7, color='black', label=labels[i])
			#print(ins_n, list(zip(list(ins_t), list(ins_u))))
		ax.plot([np.max(X)], [np.max(Y)], alpha=1) 

		ax.text(0.5, 0.75, self.type[0], bbox={'boxstyle':'Square', 'facecolor':'white'})
		ax.set_ylim(0, 8)
		ax.set_xticks(np.arange(0, self.Nb+1, 2))
		ax.set_yticks(range(0, 9, 2))
		ax.set_xticklabels(range(7))
		ax.set_yticklabels(range(0, 9, 2))
		ax.set_xlabel(r'$N$')
		ax.set_ylabel(r'$U$', labelpad=20)

		#ax.set_yticks(np.arange(0, max(U)+1, 1))
		#ax.legend(bbox_to_anchor=(1.03, 1.03), frameon=False, fontsize=17)
		#ax.set_title(r'%s-type $J/U = %.1f$' % (type[0], self.JU), loc='left', fontdict={'fontsize':'medium'})
		#ax.set_title('%s-type' % type[0], pad=20, fontdict={'fontsize':'medium'})
		
		plt.show()

		return cf

"""
	def ShowPhaseCheck(self, type1, type2):
		fn1_list = [self.path_output + fn for fn in os.listdir(self.path_output) if re.match('band_%s' % type1, fn)]
		fn2_list = [self.path_output + fn for fn in os.listdir(self.path_output) if re.match('band_%s' % type2, fn)]
		N_list = []
		U_list = []

		for fn in fn1_list:
			fn_dict = FnDict(fn)
			N_list.append(fn_dict['N'])
			U_list.append(fn_dict['U'])

		N = sorted(list(set(N_list)))
		U = sorted(list(set(U_list)))
		m1 = np.zeros((len(U), len(N)))
		m2 = np.zeros((len(U), len(N)))
		e1 = np.zeros((len(U), len(N)))
		e2 = np.zeros((len(U), len(N)))

		for fn1, fn2 in zip(fn1_list, fn2_list):
			fn1_dict = FnDict(fn1)
			fn2_dict = FnDict(fn2)
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
		fig.savefig('%s/phasec_%s%s.pdf' % (re.sub('output', 'diagram', self.path_output), type1, type2))
		plt.show()

	def ShowSolution(self, type, N, U):
		N = float(N)
		U = float(U)

		fn = [self.path_output + f for f in os.listdir(self.path_output) if re.search('sol_%s_N%.1f_U%.1f' % (type, N, U), f)][0]

		f = open(fn, 'r')
		data = np.genfromtxt(f)
		f.close()

		n_max = 4.2
		m_max = 2.2
		m_min = -2.2

		fig, ax = plt.subplots(1, 2, figsize=(10, 5), dpi=600)

		for i, n in enumerate([2, 4, 6]):
			ax[0].plot(data[:, 0], data[:, n],   lw=abs(i-4), ls=self.t2g_ls[i], marker=self.t2g_mk[i], ms=abs(n-10), color=self.t2g_color[i], label=self.t2g_label[i])
			ax[1].plot(data[:, 0], data[:, n+1], lw=abs(i-4), ls=self.t2g_ls[i], marker=self.t2g_mk[i], ms=abs(n-10), color=self.t2g_color[i], label=self.t2g_label[i])

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
"""
