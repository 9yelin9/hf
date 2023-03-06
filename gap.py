# hf3/acc.py : machine accuracy

import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import re
import numpy as np

T_dict = {'1':'a', '2':'c', '3':'g'}
save = 'baoso3_dU1.0_gap'

with open('gap.txt', 'r') as f: data = np.genfromtxt(f)

for i in range(data.shape[0]):
	T = data[i][0] * 10
	J = data[i][1]
	N = data[i][2]
	U = data[i][3]

	os.system('./mod/hf3 baoso3 %s %.1f %.1f %s%d %.1f %.1f' % (save, J, 0, T_dict[str(int(T//10))],T%10, N, U))
