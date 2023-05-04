# pyhf/outhf.py : class to generate output

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .com import ReadConfig, FnDict, GroundOnly

class OutHF:
	def __init__(self, save, strain, type, JU, SOC):
		self.dim = 3
		self.Nc  = 3

		_, _, self.Nkb = ReadConfig(self.dim)

		self.save   = save
		self.strain = strain
		self.type   = type
		self.JU     = float(JU)
		self.SOC    = float(SOC)

		self.Nb = 6 if re.search('F', self.type) else 12

		self.path_output = 'output/%s/%s_%s_JU%.2f_SOC%.2f/' % (save, self.strain, self.type, self.JU, self.SOC)
		self.path_save = '%s/diagram/' % self.path_output
		if os.path.isdir(self.path_output): os.makedirs(self.path_save, exist_ok=True)

		self.t2g_color = ['tab:blue', 'tab:green', 'tab:red']
		self.t2g_label = [r'$d_{xy}$', r'$d_{yz}$', r'$d_{zx}$']
		self.t2g_mk    = ['s', 'o', '^']
		self.t2g_ls    = ['-', '--', '-.']

		plt.rcParams.update({'font.size': 35})
		plt.rcParams.update({'font.family': 'sans-serif'})
		plt.rcParams.update({'font.serif': 'Helvetica Neue'})
	
	def ShowBand(self, N, U, ax, Nk, is_unfold):
		N = float(N)
		U = float(U)

		with open('input/%s/tb/kb_Nk%d.txt' % (self.strain, Nk), 'r') as f: hsp_info = f.readline()
		hsp_info = hsp_info.split()
		prev = object()
		hsp_point = np.cumsum([int(x) for x in hsp_info if x.isdigit()])
		hsp_point = np.insert(hsp_point, 0, 0)
		hsp_label = [prev:=label for label in [x for x in hsp_info if x.isalpha()] if prev!=label]
		hsp_label = [r'$\Gamma$' if label == 'G' else label for label in hsp_label]

		fn = [self.path_output+'/band_Nk%d/' % Nk+f for f in os.listdir(self.path_output+'/band_Nk%d' % Nk)\
				if re.search('N%.1f_U%.1f' % (N, U), f)][0]
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
			axins.set_ylim([-0.05, 1.2])
			axins.tick_params(axis='both', which='major', labelsize=18)
			#axins.yaxis.tick_right()
		else:
			for i in range(self.Nb): ax.plot(x, data[:, i], color='tab:blue')
		
		ax.grid(True, axis='x')
		ax.set_xticks(hsp_point)
		ax.set_xticklabels(hsp_label)
		ax.set_ylabel(r'$E-E_{F}$')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, fn
	
	def ShowDOS(self, N, U, ax, e_min, e_max, ep):
		N = float(N)
		U = float(U)

		fn = [self.path_output+'/dos_ep%.2f/' % ep+f for f in os.listdir(self.path_output+'/dos_ep%.2f' % ep)\
				if re.search('N%.1f_U%.1f' % (N, U), f)][0]
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
		ax.legend(fontsize=30, labelspacing=0.02, handletextpad=0.3, handlelength=1.0, borderpad=0.1, borderaxespad=0.1, frameon=False, loc='lower right')

	def ShowBandDOS(self, N, U, Nk=0, ep=0.02, is_unfold=0):
		N  = float(N)
		U  = float(U)
		ep = float(ep)

		if Nk: Nk = int(Nk)
		else:  Nk = self.Nkb

		fig, ax = plt.subplots(1, 2, width_ratios=[3, 1], figsize=(10, 5), constrained_layout=True)

		e_min, e_max, fn = self.ShowBand(N, U, ax[0], Nk, is_unfold=int(is_unfold))
		self.ShowDOS(N, U, ax[1], e_min, e_max, ep)

		fname = self.path_save + re.sub('\S+/', '', re.sub('txt', 'png', fn))
		fig.savefig(fname, dpi=300)
		print(fname)
		plt.show()

	def ShowEnergyMag(self, N, xmin=0, xmax=0, ymin=0, ymax=0):
		N = float(N)
		xmin = float(xmin)
		xmax = float(xmax)
		ymin = float(ymin)
		ymax = float(ymax)

		dU = float(re.sub('dU', '', re.search('dU\d[.]?\d*', self.save).group()))
		UF = int(re.sub('UF', '', re.search('UF\d[.]?\d*', self.save).group()))
		u_list = np.arange(0, UF+dU, dU)

		save_list = ['output/%s/%s' % (self.save, s) for s in os.listdir('output/%s' % self.save)\
				if re.search('%s\d_JU%.2f' % (self.type[0], self.JU), s)]

		fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

		e_list = []
		m_list = []
		for s, mk in zip(save_list, ['s', 'o', '^']):
			fn_list = sorted(['%s/band_Nk%d/%s' % (s, self.Nkb, f) for f in os.listdir('%s/band_Nk%d' % (s, self.Nkb))\
					if re.search('N%.1f_' % N, f)])
			e, m = np.array([(FnDict(fn)['e'], FnDict(fn)['m']) for fn in fn_list]).T

			label = re.sub('_', '', re.search('%s\d_' % self.type[0], s).group())
			ax[0].plot(u_list, e, marker=mk, ms=12, label=label)
			ax[1].plot(u_list, m, marker=mk, ms=12, label=label)

			e_list.append(e)
			m_list.append(m)

		m_list = [m_list[i][j] for i, j in zip(np.argmin(e_list, axis=0), range(len(u_list)))]
		ax[1].scatter(u_list, m_list, marker='*', color='w', zorder=2)

		for axi in ax:
			axi.grid(True)
			axi.set_xticks(u_list)
			axi.set_xlabel(r'$U$')
		ax[0].set_ylabel(r'$E$')
		ax[1].set_ylabel(r'$m$')
		ax[1].yaxis.tick_right()
		ax[1].yaxis.set_label_position('right')

		xlim = (xmin-0.2, xmax+0.2) if xmin+xmax > 1e-6 else (np.min(u_list)-0.5, np.max(u_list)+0.5)
		ylim = (ymin-0.2, ymax+0.2) if ymin+ymax > 1e-6 else (np.min(m_list)-0.5, np.max(m_list)+0.5)
		ax[0].set_xlim(xlim)
		ax[1].set_xlim(xlim)
		ax[1].set_ylim(ylim)

		ax[0].set_title(r'$N$ = %.1f' % N, loc='left', fontsize='small')
		ax[0].legend(fontsize='small', labelspacing=0.02, handletextpad=0.3, handlelength=1.0, borderpad=0.2, borderaxespad=0.1)
			
		fname = self.path_save + 'energy_N%.1f.png' % N
		fig.savefig(fname, dpi=300)
		print(fname)
		plt.show()

	def ShowPhase(self, specific_init=0, xmin=0, xmax=0, ymin=0, ymax=0):
		xmin = float(xmin)
		xmax = float(xmax)
		ymin = float(ymin)
		ymax = float(ymax)

		tol_gap = 0.1
		tol_m   = 0.1

		dN, NF = (0.1, 6) if re.search('F', self.type) else (0.2, 12)
		n_list = np.arange(dN, NF, dN)

		dU = float(re.sub('dU', '', re.search('dU\d[.]?\d*', self.save).group()))
		UF = int(re.sub('UF', '', re.search('UF\d[.]?\d*', self.save).group()))
		u_list = np.arange(0, UF+dU, dU)

		pat_type = self.type if specific_init else '%s\d' % self.type[0]
		save_list = ['output/%s/%s' % (self.save, s) for s in os.listdir('output/%s' % self.save)\
				if re.search('%s_%s_JU%.2f' % (self.strain, pat_type, self.JU), s)]
		fn_list = sorted(['%s/band_Nk%d/%s' % (s, self.Nkb, f) for s in save_list for f in os.listdir('%s/band_Nk%d' % (s, self.Nkb))])
		grd_idx = GroundOnly(fn_list)

		m = []
		ins = []
		for u in u_list:
			m.append([0, u, 0])
			for n in n_list:
				fn = [fn_list[i] for i in grd_idx if re.search('N%.1f_U%.1f' % (n, u), fn_list[i])][0]
				m.append([FnDict(fn)['N'], FnDict(fn)['U'], abs(FnDict(fn)['m'])])
				if FnDict(fn)['gap'] > tol_gap: ins.append([FnDict(fn)['N'], FnDict(fn)['U'], [FnDict(fn)['type']]])
			m.append([NF, u, 0])

		m = np.array(m)
		ins = pd.DataFrame(ins, columns=['N', 'U', 'type'])

		X = np.reshape(m[:, 0], (len(u_list), len(n_list)+2))
		Y = np.reshape(m[:, 1], (len(u_list), len(n_list)+2))
		Z = np.reshape(m[:, 2], (len(u_list), len(n_list)+2))

		fig, ax = plt.subplots(figsize=(8, 5))

		#ct = ax.contour(N, U, m, levels=[tol_m], colors='w', linestyles='dotted')
		#if abs(tol_m - 0.1) < 1e-6: ax.clabel(ct, ct.levels, inline=True, fmt='%.1f', fontsize=16)

		cf = ax.contourf(X, Y, Z, levels=np.linspace(0, self.Nb//2, 10), cmap='Blues_r')
		cb = plt.colorbar(cf, format='%.1f')
		cb.set_ticks([0, self.Nb//2])
		cb.set_ticklabels(['0', '3'])

		labels = ['_nolegend_' for _ in range(self.Nb)]
		labels[0] = 'Insulator'

		for i, ins_n in enumerate(np.unique(ins['N'])):
			idcs = np.where(np.abs(ins['N'] - ins_n) < 1e-6)
			ins_u = np.array(ins['U'])[idcs]
			ins_t = np.array(ins['type'])[idcs]
			ax.plot([ins_n, ins_n], [np.min(ins_u), np.max(ins_u)], lw=7, color='black', label=labels[i])
		ax.plot([np.max(X)], [np.max(Y)], alpha=1) 

		ax.text(0.5, 0.75, self.type[0], bbox={'boxstyle':'Square', 'facecolor':'white'})
		ax.set_xticks(range(0, NF+1, int(dN*10)), labels=range(7))
		ax.set_yticks(range(0, UF+1, 2))
		ax.set_xlabel(r'$N$')
		ax.set_ylabel(r'$U$', labelpad=20)

		xlim = (xmin-0.2, xmax+0.2) if xmin+xmax > 1e-6 else (np.min(X), np.max(X))
		ylim = (ymin-0.2, ymax+0.2) if ymin+ymax > 1e-6 else (np.min(Y), np.max(Y))
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)

		fname = self.path_save + 'phase_%s.png' % self.type if specific_init else self.path_save + 'phase.png'
		fig.savefig(fname, dpi=300)
		print(fname)
		plt.show()

		return cf

