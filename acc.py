# hf3/acc.py : machine accuracy

import os
os.environ['OMP_NUM_THREADS'] = '16'
os.environ['OPENBLAS_NUM_THREADS'] = '16'

import re
import h5py
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyhf3 import draw, magstr

from sklearn.metrics import accuracy_score, confusion_matrix

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-d', '--data', type=str)
parser.add_argument('-m', '--mcn', type=str)
parser.add_argument('-e', '--eta', type=float)
parser.add_argument('-b', '--bins', type=int)
args = parser.parse_args()                                                                     

save = 'baoso3_dU1.0'
kind = 'MK'
ms = magstr.MagStr(save, num_thread=4)

print(args.bins, type(args.bins))
if args.bins != 1024: d1n = '%s_%s_eta%.2f_bins%d.h5' % (args.data, kind, args.eta, args.bins)
else:                 d1n = '%s_%s_eta%.2f.h5' % (args.data, kind, args.eta)

with open('%s/acc_%s_%s_%s_eta%.2f_bins%d.txt' % (ms.path_save, args.mcn, args.data, kind, args.eta, args.bins), 'w') as f:
	f.write('%22s%22s%22s%22s\n' % ('eta2', 'acc_a', 'acc_c', 'acc_g'))

	for eta2 in np.arange(0.1, 0.62, 0.02):
		d2n = '%s_%s_eta%.2f_bins%d.h5' % (args.data, kind, eta2, args.bins)
		y_test, y_pred, _, _, _, _ = ms.Validate(d1n, d2n, args.mcn)
		mat = confusion_matrix(y_test, y_pred)

		f.write('%22.16f' % eta2)
		for i in range(3): f.write('%22.16f' % (mat[i][i]/np.sum(mat[i])))
		f.write('\n')
