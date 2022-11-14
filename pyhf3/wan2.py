# pyhf3/wan2.py : convert wannier90 data to something.txt

import re
import numpy as np

class Wan2:
	def __init__(self, path_input):
		self.path_input = path_input
		dim = 3

	def Wan2Lat(self):
		f_wan = open('%s/wannier90_hr.dat' % (self.path_input), 'r')
		f_lat = open('%s/lat.txt' % (self.path_input), 'w')

		pat_site = '[-]?\d+\s+'
		pat_obt = '[-]?\d+\s+'
		pat_t = '[-]?\d+[.]\d+\s+'

		pat = 3*pat_site + 2*pat_obt + 2*pat_t
		len = 0

		for line in f_wan:
			if re.search(pat, line): len += 1
		f_lat.write("%d\n" % (len))

		f_wan.seek(0, 0)
		for line in f_wan:
			if re.search(pat, line): f_lat.write(line) 

		f_wan.close()
		f_lat.close()

	def Wan2Info(self, lat_type, dim=3):
		f_wan  = open('%s/wannier90.win' % self.path_input, 'r')
		f_info = open('%s/info.txt' % self.path_input, 'w')

		lat_dict = {
			'sc' : [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
			'fcc': [[0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0.5, 0.5]]
		}

		# num_wann
		for line in f_wan:
			if re.search('num_wann', line):
				f_info.write(line + '\n')
				break

		# a
		f_info.write('begin a\n')
		for i in range(dim):
			for j in range(dim):
				f_info.write('%f ' % lat_dict[lat_type][i][j])
			f_info.write('\n')
		f_info.write('end a\n\n')

		# kpoint_path
		is_path = 0
		for line in f_wan:
			if re.search('begin kpoint_path', line): is_path = 1
			if is_path: f_info.write(line)
			if re.search('end kpoint_path', line): break

		f_wan.close()
		f_info.close()

	def Info2Path(self, max_points=1024, dim=3):
		max_points -= 1

		f_info = open('%s/info.txt' % self.path_input, 'r')
		f_path = open('%s/path.txt' % self.path_input, 'w')

		# a
		a = []
		is_a = 0
		for line in f_info:
			if re.search('end a', line): break
			if is_a: a.append(list(map(float, line.strip().split(' '))))
			if re.search('begin a', line): is_a = 1

		# b
		b = []
		V = np.dot(a[0], np.cross(a[1], a[2]))
		for i in range(dim):
			b.append(2*np.pi * np.cross(a[i-2], a[i-1]) / V)
		
		# kpoint_path
		path = []
		is_path = 0
		for line in f_info:
			if re.search('end kpoint_path', line): break
			if is_path:
				path.append(line.strip().split()[:dim+1])
				path.append(line.strip().split()[dim+1:])
			if re.search('begin kpoint_path', line): is_path = 1

		# kpoints
		pis = []
		pfs = []
		dists = []
		idxs = []
		for i in range(len(path)//2):
			ki = [0 for _ in range(dim)]
			kf = [0 for _ in range(dim)]
			
			for j in range(dim):
				ki = np.add(ki, np.dot(float(path[2*i][1+j]),   b[j]))
				kf = np.add(kf, np.dot(float(path[2*i+1][1+j]), b[j]))

			pi = np.dot(list(map(float, path[2*i][1:])),   2*np.pi)
			pf = np.dot(list(map(float, path[2*i+1][1:])), 2*np.pi)

			pis.append(list(pi))
			pfs.append(list(pf))
			dists.append(np.linalg.norm(ki - kf))
			idxs.append(path[2*i+1][0])

		dists = [int(max_points * d / sum(dists)) for d in dists]
		if sum(dists) != max_points:	
			dists[dists.index(max(dists))] += max_points - sum(dists)

		dsum = 0
		f_path.write('%s 0 ' % ('#' + path[0][0]))
		for d, p in zip(dists, idxs):
			dsum = dsum + d
			f_path.write('%s %d ' % (p, dsum))
		f_path.write('\n')

		for pi, pf, d in zip(pis, pfs, dists):
			for i in range(d):
				for j in range(dim):
					f_path.write('%.12f ' % (pi[j] + (pf[j] - pi[j]) * i / d))
				f_path.write('\n')

		for j in range(dim):
			f_path.write('%.12f ' % (pi[j] + (pf[j] - pi[j])))
		f_path.write('\n')

		f_info.close()
		f_path.close()
		
		# txt 2 bin
		f_path = open('%s/path.txt' % self.path_input, 'r')
		f_path_bin = open('%s/kb.bin' % self.path_input, 'wb')

		kpoints = np.genfromtxt(f_path, skip_header=1)[:, :]
		kpoints_bin = bytearray(kpoints)
		f_path_bin.write(kpoints_bin)

		f_path.close()
		f_path_bin.close()
