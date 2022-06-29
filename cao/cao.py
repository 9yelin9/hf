# cao/cao.py : make input and output of CuAl2O4 model

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

from pyhf3 import wan2lat
from pyhf3 import draw

path = [0, 150, 350, 400, 500, 575, 800]
path_label = ['L', 'G', 'X', 'W', 'L', 'K', 'G']

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--values', nargs='+', default=['f', 0.0, 0.0, 5.0, 0.0], help='<type> <K> <JU> <SOC> <N> <U>')
parser.add_argument('-w', '--wan2lat', type=int, default=None)
parser.add_argument('-t', '--tb', type=int, default=None)
parser.add_argument('-c', '--cvg', type=int, default=None)
parser.add_argument('-b', '--band', type=int, default=None)
parser.add_argument('-p', '--phase', type=int, default=None)
args = parser.parse_args()

w = wan2lat.Wan2Lat()
if args.wan2lat: w.Extract()

d = draw.Draw(args.values[0], args.values[1], args.values[2], path, path_label)
if args.tb: d.DrawTB()
if args.cvg: pass
if args.band: d.DrawBandDos(args.values[3], args.values[4])
if args.phase: pass
