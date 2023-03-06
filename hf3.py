# hf3/hf3.py : make input and output of hf3 model

import os
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'

import re
import sys
import scipy
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyhf3.read import ReadInfo

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-n',  '--name', type=str, default='baoso3_dU1.0', help='<name of material>')
parser.add_argument('-w2', '--wan2',   nargs='+', default=None, help='l                                         : Wan2Lat\n'\
																	+'sl                                        : Lat2SLat\n'\
																	+'i  <ltype> <dim>                          : Wan2Info\n'\
																	+'p  <pmax> <dim>                           : Info2Path')
parser.add_argument('-dr', '--draw',   nargs='+', default=None, help='b  <J/U> <SOC> <type> <N> <U>             : DrawBandDOS\n'\
																	+'bu <J/U> <SOC> <type> <N> <U>             : DrawBandDOS(unfold)\n'\
																	+'bt <J/U> <SOC> <type>                     : DrawBandTB\n'\
																	+'bc <J/U> <SOC> <type>                     : DrawBandCheck\n'\
																	+'p  <J/U> <SOC> <type> <tol_gap> <tol_m>   : DrawPhase\n'\
																	+'pc <J/U> <SOC> <type1> <type2>            : DrawPhaseCheck\n'\
																	+'s  <J/U> <SOC> <type> <N> <U>             : DrawSolution\n')
parser.add_argument('-ms', '--magstr', nargs='+', default=None, help='e                                         : GenEnergy\n'\
																	+'b  <btype>                                : GenBand\n'\
																	+'dk <btype> <eta>                          : GenDOSK\n'\
																	+'dl <btype> <eta>                          : GenDOSL\n'\
																	+'dd <dtype> <eta>                          : GenDOSD\n'\
																	+'p  <dtype> <eta> <bins>                   : GenPeak\n')
parser.add_argument('-mcn', '--mcname', type=str, default=None)
args = parser.parse_args()                                                                     

# wan2
if args.wan2:
	from pyhf3 import wan2
	w2 = wan2.Wan2(args.name)

	if   args.wan2[0] == 'l':  w2.Wan2Lat()
	elif args.wan2[0] == 'sl': w2.Lat2SLat()
	elif args.wan2[0] == 'i':  w2.Wan2Info(args.wan2[1])
	elif args.wan2[0] == 'p':  w2.Info2Path()
	sys.exit()

# draw
if args.draw:
	from pyhf3 import draw
	dr = draw.Draw(args.name, args.draw[1], args.draw[2])

	if   args.draw[0] == 'b':  dr.DrawBandDOS(args.draw[3], args.draw[4], args.draw[5])
	elif args.draw[0] == 'bu': dr.DrawBandDOS(args.draw[3], args.draw[4], args.draw[5], 1)
	elif args.draw[0] == 'bt': dr.DrawBandTB(args.draw[3])
	elif args.draw[0] == 'bc': dr.DrawBandCheck(args.draw[3])
	elif args.draw[0] == 'p':  dr.DrawPhase(args.draw[3])
	elif args.draw[0] == 'pc': dr.DrawPhaseCheck(args.draw[3], args.draw[4])
	elif args.draw[0] == 's':  dr.DrawSolution(args.draw[3], args.draw[4], args.draw[5])
	sys.exit()

# magstr
if args.magstr:
	from pyhf3 import magstr
	ms = magstr.MagStr(args.name)

	if   args.magstr[0] == 'e':  ms.GenEnergy()
	elif args.magstr[0] == 'b':  ms.GenBand(args.magstr[1])
	elif args.magstr[0] == 'dk': os.system('./mod/dos %s d %sK %s' % (args.name, args.magstr[1], args.magstr[2]))
	elif args.magstr[0] == 'dl': ms.GenDOSL(args.magstr[1], args.magstr[2])
	elif args.magstr[0] == 'dd': ms.GenDOSD(args.magstr[1], args.magstr[2])
	elif args.magstr[0] == 'p':  os.system('./mod/dos %s p %s %s %s' % (args.name, args.magstr[1], args.magstr[2], args.magstr[3]))
	sys.exit()

###############################################################################################################################################

from pyhf3 import magstr
from sklearn.metrics import accuracy_score
ms = magstr.MagStr(args.name, num_thread=8)

#for mcn in ['rf', 'xgb', 'lgb', 'cat']:
mcn = args.mcname
for eta in np.arange(0.1, 0.3, 0.0001):
	d1n = 'dos_dt1_eta%.4f.h5' % eta
	d2n = 'dos_dt0_eta0.1000.h5'
	y_test, y_pred, _, _, df2_test, _ = ms.Validate(d1n, d2n, mcn, resamp='none', verbose=False)

	acc = accuracy_score(y_test, y_pred)
	if acc > 0.5:
		print('\n-------------------- %s | %s | %s --------------------' % (mcn, d1n, d2n))
		y_test, y_pred, _, _, df2_test, _ = ms.Validate(d1n, d2n, mcn, resamp='none', verbose=True)
		print(df2_test)
