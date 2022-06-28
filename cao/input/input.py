# hf3/cao/input.py : make input files

import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/9yelin9/hf3')

parser = argparse.ArgumentParser()
parser.add_argument('-w', '--wan2lat', type=int, default=None)
parser.add_argument('-d', '--draw', type=str, default=None)
args = parser.parse_args()

if args.wan2lat:
	from Wan2Lat import Wan2Lat

	w = Wan2Lat()
	w.Extract()

if args.draw:
	from Draw import Draw

	d = Draw(type=args.draw)
	d.DrawBandTB()
