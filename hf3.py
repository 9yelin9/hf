# hf3/hf3.py : make input and output of hf3 model

import os
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyhf3 import wan2lat
from pyhf3 import ml
from pyhf3 import draw

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--material', type=str, default='boo')
parser.add_argument('-wl', '--wan2lat', type=int, default=None)
parser.add_argument('-ml', '--ml', type=int, default=None)
parser.add_argument('-v', '--values', nargs='+', default=None, help='<JU> <SOC> <type>')
#parser.add_argument('-c', '--cvg', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-b', '--band', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-bu', '--band_uf', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-p', '--phase', type=int, default=None)
args = parser.parse_args()

f_info = open('%s/input/info.txt' % (args.material), 'r')
basis      = list(map(int, f_info.readline().strip().split(',')))
path       = list(map(int, f_info.readline().strip().split(',')))
path_label = f_info.readline().strip().split(',')
f_info.close()

path[-1] -= 1

if args.wan2lat:
	w = wan2lat.Wan2Lat()
	w.Run(args.material)
	sys.exit()

if args.ml:
	m = ml.ML()
	m.Run(args.material, basis, path, path_label)
	sys.exit()

d = draw.Draw(args.material, args.values[0], args.values[1], args.values[2], basis, path, path_label)
#if args.cvg: d.DrawCvg(args.cvg[0], args.cvg[1])
if args.band: d.DrawBandDos(args.band[0], args.band[1])
if args.band_uf: d.DrawBandDos(args.band_uf[0], args.band_uf[1], 1)
if args.phase: d.DrawPhase()
