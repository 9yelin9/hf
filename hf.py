#!/home/9yelin9/.local/bin/python3

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
parser.add_argument('-i', '--inhf',  nargs='+', default=None, help='l  <strain>                                                   : GenLat\n'\
																  +'l2 <strain>                                                   : GenLat2\n'\
																  +'kg <strain> <type> [Nkg1]                                     : GenKG\n'\
																  +'kb <strain> <type> [Nkb]                                      : GenKB\n')
parser.add_argument('-o', '--outhf', nargs='+', default=None, help='b  <save> <strain> <type> <JU> <N> <U> [Nk] [eta] [is_unfold] : ShowBandDOS\n'\
																  +'e  <save> <strain> <type> <JU> <N>                            : ShowEnergyMag\n'\
																  +'p  <save> <strain> <type> <JU> [specific_init]                : ShowPhase\n')
parser.add_argument('--lim', nargs='+', type=float, help='[xmin] [xmax] [ymin] [ymax]')
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
	elif args.outhf[0] == 'e': oh.ShowEnergyMag(*args.outhf[5:], *args.lim)
	elif args.outhf[0] == 'p': oh.ShowPhase(*args.lim)
	sys.exit()

parser.print_help()
