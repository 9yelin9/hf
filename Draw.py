# Draw.py : draw diagrams

import os
import re
import sys
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt

"""
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--values', nargs='+', default=['f', 0.0, 0.0, 5.0, 0.0], help='<type> <JU> <SOC> <N> <U>')
parser.add_argument('-c', '--cvg', type=int, default=0)
parser.add_argument('-b', '--band', type=int, default=1)
parser.add_argument('-p', '--phase', type=int, default=0)
args = parser.parse_args()
"""

class Draw:
	def __init__(self, type='f', JU=0.0, SOC=0.0):
		self.type = type
		self.K = 16
		self.JU = float(JU)
		self.SOC = float(SOC)
		self.obt = 3
		self.basis = 6 if re.search('f', type) else 12
		self.dir = 'output/K%d_JU%.2f_SOC%.2f/' % (self.K, self.JU, self.SOC)
		self.title = '%s K=%d J/U=%.2f SOC=%.2f\n' % (type, self.K, self.JU, self.SOC)
		self.colors=['b', 'g', 'r']
		self.labels=['xy', 'yz', 'zx']
	
	def DrawBandTB(self):
		f_name = 'band_%s.txt' % self.type 
		f = open(f_name, 'r')
		arr = np.genfromtxt(f)
		f.close()

		fig, ax = plt.subplots()

		for i in range(1, self.basis+1):
			ax.plot(arr[:, 0], arr[:, i], color='tab:blue')

		plt.suptitle(f_name)	
		fig.tight_layout()
		plt.show()

	def DrawBand(self, N, U, ax, path, path_label):
		N = float(N)
		U = float(U)

		f_name = [f for f in os.listdir(self.dir) if re.search('band_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('fermi', '', re.search('fermi\d+[.]\d+', f_name).group()))

		f = open(self.dir + f_name, 'r')
		arr = np.genfromtxt(f)
		f.close()

		e_min = np.min(arr[:, 1:]) - fermi - 0.1
		e_max = np.max(arr[:, 1:]) - fermi + 0.1

		ax.axhline(y=0.0, ls=':', lw=2, color='dimgrey')

		for i in range(1, self.basis+1):
			ax.plot(arr[:, 0], arr[:, i] - fermi, color='tab:blue')

		ax.grid(True)
		ax.set_xticks([0, 100, 200, 340, 510])
		ax.set_xticklabels(['G', 'X', 'M', 'G', 'R'])
		ax.set_xlabel('Path')
		ax.set_ylabel('Energy')
		ax.set_ylim(e_min, e_max)

		return e_min, e_max, f_name 
	
	def DrawDos(self, N, U, ax, e_min=None, e_max=None):
		N = float(N)
		U = float(U)

		f_name = [f for f in os.listdir(self.dir) if re.search('dos_%s_N%.1f_U%.1f' % (self.type, N, U), f)][0]
		fermi = float(re.sub('fermi', '', re.search('fermi\d+[.]\d+', f_name).group()))
		title = self.dir + f_name

		f = open(self.dir + f_name, 'r')
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
	
	def DrawBandDos(self, N, U, path, path_label):
		N = float(N)
		U = float(U)

		fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=[16, 8])

		e_min, e_max, f_name = d.DrawBand(N, U, ax[0], path, path_label)
		d.DrawDos(N, U, ax[1], e_min, e_max)

		f_name_s = f_name.split(sep='_')
		f_name_s = [s for s in f_name_s if re.search('[a-zA-Z]+\d+[.]\d+', s)]

		title = self.title
		for s in f_name_s:
			s_name = re.search('[a-zA-Z]+', s).group()
			s_val = re.search('\d+[.]\d+', s).group()
			title += '%s=%s ' % (s_name, s_val)

		plt.suptitle(title)
		fig.tight_layout()
		fig.savefig('diagram/%s/%s' % (re.sub('output/', '', self.dir), re.sub('txt', 'png', f_name)))
		plt.show()
