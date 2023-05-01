import os
num_thread = 4
os.environ['OMP_NUM_THREADS']      = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--inhf',  nargs='+', default=None, help='l  <strain>                                                  : GenLat\n'\
																  +'l2 <strain>                                                  : GenLat222\n'\
																  +'kg <strain> <type> [Nkg1]                                    : GenKG\n'\
																  +'kb <strain> <type> [Nkb]                                     : GenKB\n')
parser.add_argument('-o', '--outhf', nargs='+', default=None, help='b  <save> <type> <JU> <SOC> <N> <U> [eta=0.02] [is_unfold=0] : ShowBandDOS\n'\
																  +'e  <save> <type> <JU> <SOC> <N>                              : ShowEnergyMag\n'\
																  +'p  <save> <type> <JU> <SOC>                                  : ShowPhase\n')
args = parser.parse_args()                                                                     

# inhf
if args.inhf:
	from pyhf import inhf
	ih = inhf.InHF(args.inhf[1])

	if   args.inhf[0] == 'l' : ih.GenLat()
	elif args.inhf[0] == 'l2': ih.GenLat2()
	elif args.inhf[0] == 'kb': ih.GenKB(*args.inhf[2:])
	elif args.inhf[0] == 'kg': ih.GenKG(*args.inhf[2:])
	sys.exit()

# outhf
if args.outhf:
	from pyhf import outhf
	oh = outhf.OutHF(*args.outhf[1:5])

	if   args.outhf[0] == 'b': oh.ShowBandDOS(*args.outhf[5:])
	elif args.outhf[0] == 'e': oh.ShowEnergyMag(*args.outhf[5:])
	elif args.outhf[0] == 'p': oh.ShowPhase()
	sys.exit()
