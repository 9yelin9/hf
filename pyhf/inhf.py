# pyhf/inhf.py : class to generate input

import os
import re
import numpy as np

class InHF:
	def __init__(self, strain='none'):
		self.DIM = 3

		self.strain = strain
		self.path_input = 'input/%s' % self.strain
		self.path_save  = '%s/tb'    % self.path_input
		os.makedirs(self.path_save, exist_ok=True)
	
	def GenLat(self):
		fwn = '%s/wann/wannier90_hr.dat' % self.path_input
		fln = '%s/lat.txt'               % self.path_save

		pat_site = '[-]?\d+\s+'
		pat_obt  = '[-]?\d+\s+'
		pat_t    = '[-]?\d+[.]\d+\s+'
		pat      = self.DIM*pat_site + 2*pat_obt + 2*pat_t

		with open(fwn, 'r') as f:
			f_len = sum(1 for line in f if re.search(pat, line))

		with open(fwn, 'r') as fw, open(fln, 'w') as fl:
			fl.write('%d\n' % f_len)
			for line in fw:
				if re.search(pat, line): fl.write(line)

		print(fln)
	
	def GenLat2(self):
		fln = '%s/lat.txt'  % self.path_save
		fsn = '%s/lat2.txt' % self.path_save

		rule_dict = {
			'000': 0,
			'100': 3,
			'010': 6,
			'110': 9,
			'001': 12,
			'101': 15,
			'011': 18,
			'111': 21,
		}

		site0 = []
		obt0  = []
		t0    = []
		with open(fln, 'r') as f:
			f_len = int(f.readline())
			for line in f:
				lines = line.split()
				site0.append(lines[:self.DIM])
				obt0.append(lines[self.DIM:self.DIM+2])
				t0.append(lines[self.DIM+2:])
		site0 = np.array(site0, dtype='i')
		obt0  = np.array(obt0,  dtype='i')
		t0    = np.array(t0,    dtype='d')

		with open(fsn, 'w') as f:
			f.write('%d\n' % (f_len * len(rule_dict)))
			for key, val in rule_dict.items():
				site = np.add([int(k) for k in key], site0)
				obt  = np.add([val, 0], obt0)

				for s, o, t in zip(site, obt, t0):
					rule = ''.join(['%d' % (si % 2) for si in s])
					s = [si // 2 for si in s]
					o[1] += rule_dict[rule]

					f.write(''.join(['%5d'    % si for si in s]))
					f.write(''.join(['%5d'    % oi for oi in o]))
					f.write(''.join(['%12.6f' % ti for ti in t]))
					f.write('\n')

		print(fsn)

	def GenKB(self, type, Nk=1024):
		Nk = int(Nk) - 1

		fwn = '%s/wann/wannier90.win' % self.path_input
		fpn = '%s/POSCAR'             % self.path_input
		fkn = '%s/kb_Nk%d.txt'        % (self.path_save, Nk+1)

		exception = ['T', 'Y']
		path = []
		path_string = []
		with open(fwn, 'r') as f:
			is_path = 0
			for line in f:
				if re.match('end kpoint_path', line): break
				if is_path:
					lines = line.split()
					if re.search(''.join(exception), lines):
						print(lines)
					else:
						path.append([lines[1:self.DIM+1], lines[self.DIM+2:]])
						path_string.append([lines[0], lines[self.DIM+1]])
				# T-Y path 삭제
				if re.match('begin kpoint_path', line): is_path = 1
		path = np.array(path, dtype='d')

		"""
		lat_dict = {
			'sc' : [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
			'fcc': [[0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0.5, 0.5]],
		}

		A = np.array(lat_dict[ltype])
		V = 2*np.pi / np.dot(A[0], np.cross(A[1], A[2]))	
		B = [V * np.cross(A[i-2], A[i-1]) for i in range(self.DIM)]
		"""

		with open(fpn, 'r') as f:
			f.readline() # skip name
			c = float(f.readline())
			B = [c * np.fromstring(f.readline(), dtype='d', sep=' ') for _ in range(self.DIM)]

		dist_double = []
		for p in path:		
			ki = kf = np.zeros(self.DIM)
			for i in range(self.DIM):
				ki = np.add(ki, np.dot(p[0, i], B[i]))
				kf = np.add(kf, np.dot(p[1, i], B[i]))
			dist_double.append(np.linalg.norm(ki - kf))

		dist = [int(Nk * d / sum(dist_double)) for d in dist_double]
		if sum(dist) != Nk: dist[dist.index(max(dist))] += Nk - sum(dist)
		
		with open(fkn, 'w') as f:
			for p, d in zip(path_string, dist): f.write('%s\t%d\t' % ('\t'.join(p), d))
			f.write('\n')
					
			path = np.dot(2*np.pi, path)
			for p, d in zip(path, dist):
				for di in range(d):
					f.write('\t'.join(['%20.16f' % (p[0, i] + (p[1, i] - p[0, i]) * di / d) for i in range(self.DIM)]))
					f.write('\n')
			f.write('\t'.join(['%20.16f' % (path[-1, 0, i] + (path[-1, 1, i] - path[-1, 0, i])) for i in range(self.DIM)]))
			f.write('\n')

		print(fkn)
		self.GenUF(fkn, type)
	
	def GenKG(self, type, Nk1=32):
		Nk1 = int(Nk1)
		Nk  = Nk1 ** self.DIM

		fkn = '%s/kg_Nk%d.txt' % (self.path_save, Nk)
					
		ki = -np.pi
		kf =  np.pi
		k, w = np.polynomial.legendre.leggauss(Nk1)
		k = k * 0.5*(kf - ki) + 0.5*(ki + kf)
		w = w * 0.5*(kf - ki)

		with open(fkn, 'w') as f:
			f.write('ki\t%20.16f\tkf\t%20.16f\n' % (ki, kf))

			for i in range(Nk):	
				i0 = i // (Nk1 * Nk1)
				i1 = (i // Nk1) % Nk1
				i2 = i % Nk1

				f.write(''.join(['%20.16f' % k[ii] for ii in [i0, i1, i2]]))
				f.write('%20.16f\n' % (w[i0] * w[i1] * w[i2]))

		print(fkn)
		self.GenUF(fkn, type)
	
	def GenUF(self, fkn, type):	
		ktype = re.sub('k', '', re.search('k[a-z]', fkn).group())
		Nk = int(re.sub('Nk', '', re.search('Nk\d+', fkn).group()))

		fcn = 'input/config_%s.txt' % type
		fun = '%s/uf%s_Nk%d_%s.txt' % (self.path_save, ktype, Nk, type)

		with open(fkn, 'r') as f: k = np.genfromtxt(f, skip_header=1)[:, :self.DIM]
		with open(fcn, 'r') as f:
			r = []
			for line in f:
				if   re.match('Ni',  line): Ni = int(line.split()[1])
				elif re.match('Rho', line): r.append(line.split()[1:])
			r = np.array(r, dtype='d')

		with open(fun, 'w') as f:
			f.write('%d\n' % Nk)

			for ki in k:
				f.write(''.join(['%20.16f' % np.dot(r[i], ki) for i in range(Ni)]))
				f.write('\n')

		print(fun)
