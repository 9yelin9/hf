# boo/boo.py : make input and output of BaOsO3 model

import os
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

from pyhf3 import wan2lat
from pyhf3 import draw

path = [0, 100, 200, 340, 510]
path_label = ['G', 'X', 'M', 'G', 'R']

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--values', nargs='+', default=['f', 0, 0, 0, 0, 0], help='<type> <JU> <SOC> <is_unfold> <N> <U>')
parser.add_argument('-w', '--wan2lat', type=int, default=None)
parser.add_argument('-t', '--tb', type=int, default=None)
parser.add_argument('-c', '--cvg', type=int, default=None)
parser.add_argument('-b', '--band', type=int, default=None)
parser.add_argument('-p', '--phase', type=int, default=None)
args = parser.parse_args()

w = wan2lat.Wan2Lat()
if args.wan2lat: w.SeperateLat()
if not args.wan2lat: w.Extract()

d = draw.Draw(args.values[0], args.values[1], args.values[2], args.values[3], path, path_label)
if args.tb: d.DrawTB()
if args.cvg: pass
if args.band: d.DrawBandDos(args.values[4], args.values[5])
if args.phase: d.DrawPhase()
