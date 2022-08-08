# hf3/hf3.py : make input and output of hf3 model

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyhf3 import wan2lat
from pyhf3 import draw

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--material', type=str, default='boo')
parser.add_argument('-w', '--wan2lat', type=int, default=None)
parser.add_argument('-v', '--values', nargs='+', default=None, help='<type> <JU> <SOC>')
#parser.add_argument('-c', '--cvg', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-b', '--band', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-bu', '--band_unfolded', nargs='+', default=None, help='<N> <U>')
parser.add_argument('-p', '--phase', type=int, default=None)
args = parser.parse_args()

f_path = open('%s/input/path.txt' % (args.material), 'r')
path = f_path.read().splitlines()
path_label = f_path.read().splitlines()
f_path.close()

w = wan2lat.Wan2Lat()
if args.wan2lat: w.runWan2Lat(args.material)

d = draw.Draw(args.values[0], args.values[1], args.values[2], args.values[3], path, path_label)
#if args.cvg: d.DrawCvg(args.cvg[0], args.cvg[1])
if args.band: d.DrawBandDos(args.band[0], args.band[1])
if args.band_unfolded: d.DrawBandDosUnfolded(args.band_unfolded[0], args.band_unfolded[1])
if args.phase: d.DrawPhase()
