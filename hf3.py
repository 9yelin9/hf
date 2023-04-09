# hf3/hf3.py : make input and output of hf3 model

import os
num_thread = 4
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyhf3.read import ReadInfo

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-n', '--name', type=str, default='baoso3_dU1.0', help='<(HF3)material/(DMFT)directory>')
parser.add_argument('-i', '--imp', nargs='+', type=int, default=range(8), help='<(DMFT)impurities>')
parser.add_argument('-w2',     '--wan2', nargs='+', default=None, help='l                                         : Wan2Lat\n'\
							 										  +'sl                                        : Lat2SLat\n'\
																	  +'i  <ltype> <dim>                          : Wan2Info\n'\
																	  +'p  <pmax> <dim>                           : Info2Path')
parser.add_argument('-dr',     '--draw', nargs='+', default=None, help='b  <J/U> <SOC> <type> <N> <U>             : DrawBandDOS\n'\
																	  +'bu <J/U> <SOC> <type> <N> <U>             : DrawBandDOS(unfold)\n'\
																	  +'bt <J/U> <SOC> <type>                     : DrawBandTB\n'\
																 	  +'bc <J/U> <SOC> <type>                     : DrawBandCheck\n'\
																 	  +'p  <J/U> <SOC> <type> <tol_gap> <tol_m>   : DrawPhase\n'\
																 	  +'pc <J/U> <SOC> <type1> <type2>            : DrawPhaseCheck\n'\
																	  +'s  <J/U> <SOC> <type> <N> <U>             : DrawSolution\n')
parser.add_argument('-mst',  '--magstr', nargs='+', default=None, help='e                                         : GenEnergy\n'\
																	  +'b  <btype>                                : GenBand\n'\
																	  +'dk <btype> <eta>                          : GenDOSK\n'\
																	  +'dl <btype> <eta>                          : GenDOSL\n'\
																	  +'dd <dtype> <eta>                          : GenDOSD\n'\
																	  +'p  <dtype> <eta> <bins>                   : GenPeak\n')
parser.add_argument('-msp', '--magspec', nargs='+', default=None, help='e                                         : GenEnergy\n'\
																	  +'f  <dtype>                                : GenFileList\n')
parser.add_argument('-mft', '--magfeat', nargs='+', default=None, help='e                                         : GenEnergy\n'\
																	  +'l  <dtype>                                : GenList\n')
parser.add_argument('-dm',    '--dmftt', nargs='+', default=None, help='pp <AF> <U> <J> <Hz> <dtype> <cnt>        : PrintParams\n'\
																	  +'pm <AF> <U> <J> <Hz> <dtype> <cnt>        : PrintMag\n'\
																	  +'sp <AF> <U> <J> <Hz> <dtype> <cnt>        : ShowParams\n'\
																	  +'sn <AF> <U> <J> <Hz> <dtype> <cnt>        : ShowNext\n')
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
	dr = draw.Draw(args.name, *args.draw[1:3])

	if   args.draw[0] == 'b':  dr.DrawBandDOS(*args.draw[3:])
	elif args.draw[0] == 'bu': dr.DrawBandDOS(*args.draw[3:], 1)
	elif args.draw[0] == 'bt': dr.DrawBandTB(args.draw[3])
	elif args.draw[0] == 'bc': dr.DrawBandCheck(args.draw[3])
	elif args.draw[0] == 'p':  dr.DrawPhase(args.draw[3])
	elif args.draw[0] == 'pc': dr.DrawPhaseCheck(*args.draw[3:])
	elif args.draw[0] == 's':  dr.DrawSolution(*args.draw[3:])
	sys.exit()

# magstr
if args.magstr:
	from pyhf3 import magstr
	ms = magstr.MagStr(args.name)

	if   args.magstr[0] == 'e':  ms.GenEnergy()
	elif args.magstr[0] == 'b':  ms.GenBand(args.magstr[1])
	elif args.magstr[0] == 'dl': ms.GenDOSL(*args.magstr[1:])
	elif args.magstr[0] == 'dd': ms.GenDOSD(*args.magstr[1:])
	sys.exit()

# magspec
if args.magspec:
	from pyhf3 import magspec
	mp = magspec.MagSpec(name=args.name)

	if   args.magspec[0] == 'e': mp.GenEnergy()
	elif args.magspec[0] == 'f': mp.GenFileList(args.magspec[1])
	sys.exit()

# magfeat
if args.magfeat:
	from pyhf3 import magfeat
	mf = magfeat.MagFeat(name=args.name)

	if   args.magfeat[0] == 'e': mf.GenEnergy()
	elif args.magfeat[0] == 'l': mf.GenList(args.magfeat[1])
	sys.exit()

# dmftt
if args.dmftt:
	from pyhf3 import dmftt
	dm = dmftt.DMFTT(*args.dmftt[1:5], dmft_dir=args.name)

	if   args.dmftt[0] == 'pp': dm.PrintParams(*args.dmftt[5:], args.imp)
	elif args.dmftt[0] == 'pm': dm.PrintMag(*args.dmftt[5:], args.imp)
	elif args.dmftt[0] == 'sp': dm.ShowParams(*args.dmftt[5:], args.imp)
	#elif args.dmftt[0] == 'sm': dm.ShowMag(*args.dmftt[5:], args.imp)
	elif args.dmftt[0] == 'sn': dm.ShowNext(args.imp)
	sys.exit()

###

