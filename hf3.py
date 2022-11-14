# hf3/hf3.py : make input and output of hf3 model

import os
import re
import sys
import scipy
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyhf3.read import ReadInfo

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-n',  '--name',   nargs='+', default=None, help='<name_of_material>')
parser.add_argument('-w2', '--wan2',   nargs='+', default=None, help='l                              : Wan2Lat\n'\
																	+'i  <lat_type> <dim>            : Wan2Info\n'\
																	+'p  <max_points> <dim>          : Info2Path')
parser.add_argument('-dr', '--draw',   nargs='+', default=None, help='b  <J/U> <SOC> <type> <N> <U>  : DrawBandDOS\n'\
																	+'bu <J/U> <SOC> <type> <N> <U>  : DrawBandDOS(unfold)\n'\
																	+'bt <J/U> <SOC> <type>          : DrawBandTB\n'\
																	+'bc <J/U> <SOC> <type>          : DrawBandCheck\n'\
																	+'p  <J/U> <SOC> <type>          : DrawPhase\n'\
																	+'pc <J/U> <SOC> <type1> <type2> : DrawPhaseCheck\n'\
																	+'s  <J/U> <SOC> <type> <N> <U>  : DrawSolution\n')
parser.add_argument('-ms', '--magstr', nargs='+', default=None, help='b  <ptype> <even_num>               : MakeBand\n'\
																	+'d  <ptype> <bins> <eta>        : MakeDOS\n'\
																	+'rf <ptype>                     : DoRandomForest\n')
args = parser.parse_args()                                                                     

# set path_input and path_output
if args.name:
	path_input = 'input/%s/' % args.name[0].split('_')[0]
	path_output = 'output/%s/' % args.name[0]

# input
if args.wan2:
	from pyhf3 import wan2
	w2 = wan2.Wan2(path_input)

	if   args.wan2[0] == 'l': w2.Wan2Lat()
	elif args.wan2[0] == 'i': w2.Wan2Info(args.wan2[1])
	elif args.wan2[0] == 'p': w2.Info2Path()
	sys.exit()

# set info_path, info_cell 
info_path, info_cell = ReadInfo(path_input)

# set tol
tol = 0.1

# draw
if args.draw:
	from pyhf3 import draw
	dr = draw.Draw(path_input, path_output, info_path, info_cell, args.draw[1], args.draw[2], tol)

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
	ms = magstr.ML(path_output, info_path, info_cell, args.magstr[1])

	bins_range = np.arange(10, 210, 10)
	eta_range = np.arange(0.1, 1.1, 0.1)
	func_dict = {
		'rf': ms.DoRandomForest,
	}

	if   args.magstr[0] == 'b':  ms.MakeBand(tol, args.magstr[2])
	elif args.magstr[0] == 'd':  ms.MakeDOS(args.magstr[1], args.magstr[2])
	else:
		f = open('%s/magstr/acc_%s.csv' % (path_output, args.magstr[0]), 'w')
		f.write('bins,eta,acc\n')
		for bins in bins_range:
			for eta in eta_range:
				acc, y_test, y_pred = func_dict[args.magstr[0]](bins, eta)
				f.write('%d,%.2f,%f\n' % (bins, eta, acc))
		f.close()
	sys.exit()

################################################################################################################

dr = draw.Draw(path_input, path_output, info_path, info_cell, 0, 0, tol)
for n, u in zip([8, 6], [2, 5]):
	for t in ['a', 'a1', 'a2', 'c', 'c1', 'c2', 'g']:
		dr.DrawBandDOS(t, n, u, 1)
