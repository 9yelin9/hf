# hf3/runspec.py : run pyhf3/spec.py

import os
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'

import sys
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-n', '--name', type=str, default='baoso3_dU0.2')
parser.add_argument('-s', '--spec', nargs='+', help='e                : GenEnergy_\n'\
												   +'f  <dtype>       : GenFileList_\n'\
												   +'s  <dtype> <eta> : ShowSpec(random)\n')
args = parser.parse_args()                                                                     

if args.spec:
	from pyhf3 import spec
	sp = spec.Spec(name=args.name)

	if   args.spec[0] == 'e': sp.GenEnergy_()
	elif args.spec[0] == 'f': sp.GenFileList_(args.spec[1])
	elif args.spec[0] == 's': sp.ShowSpec(dtype=args.spec[1], eta=float(args.spec[2]), is_random=True)
	sys.exit()

