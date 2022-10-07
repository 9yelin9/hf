# hf3/hf3.py : make input and output of hf3 model

import os
import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyhf3.read import ReadInfo
from pyhf3.wan2lat import TransWan2Lat
from pyhf3.ml import MakeMLData
from pyhf3 import draw

parser = argparse.ArgumentParser()
parser.add_argument('-w2l',  type=int,  default=None)
parser.add_argument('-ml',   type=int,  default=None)
parser.add_argument('-draw', nargs='+', default=None, help='<JU> <SOC>')

parser.add_argument('-n', '--name', type=str, default='baoso3', help='<solid_name>')

parser.add_argument('-b',  '--band',       nargs='+', default=None, help='<type> <N> <U>')
parser.add_argument('-bu', '--bandunfold', nargs='+', default=None, help='<type> <N> <U>')
parser.add_argument('-bt', '--bandtest',   nargs='+', default=None, help='<type>')
parser.add_argument('-bc', '--bandcheck',  nargs='+', default=None, help='<type>')
parser.add_argument('-p',  '--phase',      nargs='+', default=None, help='<type>')
parser.add_argument('-pc', '--phasecheck', nargs='+', default=None, help='<type1> <type2>')
parser.add_argument('-s',  '--solution',   nargs='+', default=None, help='<type> <N> <U>')
args = parser.parse_args()

# set input_path and output_path
input_path = 'input/%s/' % args.name
output_path = 'output/%s/' % args.name

# set path_info, type_info 
path_info, type_info = ReadInfo(input_path)
for key, val in type_info.items():
	if val == 's': type_info[key] = 6
	else: type_info[key] = 12

# set tol
tol = 0.1

# wan2lat
if args.w2l:
	TransWan2Lat(input_path)
	sys.exit()

# ml
if args.ml:
	MakeMLData(output_path, path_info, type_info, tol)
	sys.exit()

# draw
if args.draw:
	d = draw.Draw(input_path, output_path, path_info, type_info, args.draw[0], args.draw[1], tol)
	if args.band: d.DrawBandDos(args.band[0], args.band[1], args.band[2])
	if args.bandunfold: d.DrawBandDos(args.bandunfold[0], args.bandunfold[1], args.bandunfold[2], 1)
	if args.bandtest: d.DrawBandTest(args.bandtest[0])
	if args.bandcheck : d.DrawBandCheck(args.bandcheck[0])
	if args.phase: d.DrawPhase(args.phase[0])
	if args.phasecheck: d.DrawPhaseCheck(args.phasecheck[0], args.phasecheck[1])
	if args.solution: d.DrawSolution(args.solution[0], args.solution[1], args.solution[2])
	sys.exit()

######################################################################################################################

type = ['a1', 'a2', 'c1', 'c2', 'g']
JU = [0, 0.1, 0.2, 0.3]

for t in type:
	for ju in JU:
		d = draw.Draw(output_path, path_info, type_info, ju, 0, tol)
		d.DrawPhase(t)
