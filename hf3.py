# hf3/hf3.py : make input and output of hf3 model

import os
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyhf3 import wan2lat
from pyhf3 import ml
from pyhf3 import check
from pyhf3 import draw

parser = argparse.ArgumentParser()
parser.add_argument('-mat', '--material', type=str, default='baoso3')
parser.add_argument('-out', '--output', type=str, default='output')

parser.add_argument('-wl', '--wan2lat', type=int, default=None)
parser.add_argument('-ml', '--ml', type=int, default=None)
parser.add_argument('-ch', '--check', nargs='+', default=None, help='<JU> <SOC> <N4ch>')
parser.add_argument('-dr', '--draw', nargs='+', default=None, help='<type> <JU> <SOC>')

#parser.add_argument('-c', '--cvg', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-b', '--band', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-u', '--banduf', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-t', '--bandtest', type=int, default=None)
parser.add_argument('-p', '--phase', type=int, default=None)
parser.add_argument('-c', '--cvg', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-m', nargs='+', default=None, help='<N>')
args = parser.parse_args()

f_info = open('%s/input/info.txt' % (args.material), 'r')
path_num   = f_info.readline()
path       = list(map(int, f_info.readline().strip().split(' ')))
path_label = f_info.readline().strip().split(' ')
f_info.close()

basis = [6, 12]

if args.wan2lat:
	w = wan2lat.Wan2Lat()
	w.Run(args.material)
	sys.exit()

if args.ml:
	m = ml.ML()
	m.Run(args.material, args.output, basis, path, path_label)
	sys.exit()

if args.check:
	c = check.Check()
	c.Run(args.material, args.output, args.check[0], args.check[1], args.check[2])
	sys.exit()

if args.draw:
	d = draw.Draw(args.material, args.output, args.draw[0], args.draw[1], args.draw[2], basis, path, path_label)
	#if args.cvg: d.DrawCvg(args.cvg[0], args.cvg[1])
	if args.band: d.DrawBandDos(args.band[0], args.band[1])
	if args.banduf: d.DrawBandDos(args.banduf[0], args.banduf[1], 1)
	if args.bandtest: d.DrawBandTest()
	if args.phase: d.DrawPhase()
	if args.cvg: d.DrawCvg(args.cvg[0], args.cvg[1])
	if args.m : d.DrawM(args.m[0])
	sys.exit()

type_list = ['sa', 'sc', 'sg']
JU_list = [0.0, 0.1, 0.2, 0.3]

for type in type_list:
	d = draw.Draw(args.material, args.output, type, 0.0, 0.0, basis, path, path_label)
	d.DrawBandDos(6.0, 2.0)
	d.DrawBandDos(8.0, 2.0)
