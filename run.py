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
parser.add_argument('-i', '--inhf', nargs='+', default=None, help='l                    : GenLat\n'\
																 +'kb <ltype> [Nk=1024] : GenKB\n'\
																 +'kg [Nk1=32]          : GenKG\n'\
																 +'uf <fkn> <type>      : GenUF\n')
args = parser.parse_args()                                                                     

# inhf
if args.inhf:
	from pyhf import inhf
	ih = inhf.InHF()

	if   args.inhf[0] == 'l' : ih.GenLat()
	elif args.inhf[0] == 'kb': ih.GenKB(*args.inhf[1:])
	elif args.inhf[0] == 'kg': ih.GenKG(*args.inhf[1:])
	elif args.inhf[0] == 'uf': ih.GenUF(*args.inhf[1:])
	sys.exit()
