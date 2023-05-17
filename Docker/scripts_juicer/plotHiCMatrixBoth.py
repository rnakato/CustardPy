#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.stats as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sys.path.append("/home/git/script_rnakato/Hi-C")
from HiCmodule import *
cm = generate_cmap(['#FFFFFF', '#d10a3f'])

usage = 'Usage: \n    python {} <matrixdir> <output prefix> <start> <end> <label>'.format(__file__)
arguments = sys.argv
if len(arguments) == 1:
    print(usage)
    exit()

filename = sys.argv[1]
output = sys.argv[2]
s1 = int(sys.argv[3])
e1 = int(sys.argv[4])
s2 = int(sys.argv[5])
e2 = int(sys.argv[6])
label = sys.argv[7]

data = pd.read_csv(filename, delimiter='\t', index_col=0)
resolution = data.index[1] - data.index[0]

s1bin = int(s1 / resolution)
e1bin = int(e1 / resolution)
s2bin = int(s2 / resolution)
e2bin = int(e2 / resolution)

# Total read normalization
data = data * 1000000 / np.nansum(data)

fig = plt.figure(figsize=(8, 8))
ymax = np.sqrt(data.unstack().max())/2
plt.imshow(data.values[s1bin:e1bin,s2bin:e2bin], clim=(0, ymax), cmap=cm)
plt.title(label)
plt.savefig(output)
