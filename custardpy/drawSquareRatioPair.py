#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy import ndimage
from HiCmodule import JuicerMatrix
from generateCmap import *
from PlotModule import *
from DirectionalFreqRatio import *

#import pdb
def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end", help="end bp", type=int, default=1000000)
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=1)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=-1)
    parser.add_argument("--xsize", help="xsize for figure", type=int, default=3)
    parser.add_argument("--ysize", help="ysize (* times of samples)", type=int, default=3)

    args = parser.parse_args()
#    print(args)

    dirs = []
    labels = []
    for input in args.input:
        dirs.append(input[0])
        if (len(input) >1):
            labels.append(input[1])
        else:
            labels.append("")

    chr = args.chr
    resolution = args.resolution
    type = args.type
    figstart = args.start
    figend = args.end
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    vmax = args.vmax
    vmin = args.vmin

    print (chr)
    print (resolution)

    samples = []
    for dir in dirs:
#        observed = dir + "/Matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        samples.append(JuicerMatrix("RPM", dir, resolution))

    smooth_median_filter = 3
    mt1 = make3dmatrixRatio([samples[0], samples[1]], smooth_median_filter)
    mt2 = make3dmatrixRatio([samples[2], samples[3]], smooth_median_filter)
#    Combined = np.concatenate((mt1), axis=0)

    matrix = np.triu(mt1[0], k=1) + np.tril(mt2[0], k=-1)

    if (labels[1] != "" and labels[3] != ""):
        label = labels[1] + "-" + labels[3]
    else:
        label = ""

    plt.subplot(1, 1, 1)
    drawHeatmapSquare(plt, matrix, resolution,
                      figstart=figstart, figend=figend,
                      vmax=vmax, vmin=vmin,
                      cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                      label=label, xticks=False)
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
