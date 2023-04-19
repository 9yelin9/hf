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
parser.add_argument('-s', '--save', default=None, help='<save>')
parser.add_argument('-i', '--inhf',  nargs='+', default=None, help='l                                         : GenLat\n'\
																  +'l222                                      : GenLat222\n'\
																  +'kb <ltype> [Nk=1024]                      : GenKB\n'\
																  +'kg [Nk1=32]                               : GenKG\n'\
																  +'uf <fkn> <type>                           : GenUF\n')
parser.add_argument('-o', '--outhf', nargs='+', default=None, help='b <type> <JU> <SOC> <N> <U> [is_unfold=0] : ShowBandDOS\n'\
																  +'p                                         : ShowPhase\n')
args = parser.parse_args()                                                                     

# inhf
if args.inhf:
	from pyhf import inhf
	ih = inhf.InHF()

	if   args.inhf[0] == 'l'   : ih.GenLat()
	elif args.inhf[0] == 'l222': ih.GenLat222()
	elif args.inhf[0] == 'kb'  : ih.GenKB(*args.inhf[1:])
	elif args.inhf[0] == 'kg'  : ih.GenKG(*args.inhf[1:])
	elif args.inhf[0] == 'uf'  : ih.GenUF(*args.inhf[1:])
	sys.exit()

# outhf
if args.outhf:
	from pyhf import outhf
	oh = outhf.OutHF(args.save, *args.outhf[1:4])

	if   args.outhf[0] == 'b': oh.ShowBandDOS(*args.outhf[4:])
	#elif args.outhf[0] == 'p': oh.ShowPhase(*args.outhf[4:])
	sys.exit()
