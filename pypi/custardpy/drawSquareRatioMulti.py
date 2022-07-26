#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import ndimage
from HiCmodule import JuicerMatrix
from custardpy.generateCmap import *
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
#    parser.add_argument("--ysize", help="ysize (* times of samples)", type=int, default=3)

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
    figsize_y = int((figend-figstart)/2000000)

    print (chr)
    print (resolution)
    samples = []
    for dir in dirs:
        observed = dir + "/Matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        samples.append(JuicerMatrix("RPM", observed, resolution))

    ### Plot
    smooth_median_filter = 3
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)

    nsample = len(samples)
    plt.figure(figsize=(nsample*args.xsize/2, figsize_y))

    for i, sample in enumerate(EnrichMatrices):
        # Hi-C Map
        plt.subplot(2, len(samples)/2, i+2)
        drawHeatmapSquare(plt, sample, resolution,
                          figstart=figstart, figend=figend,
                          vmax=vmax, vmin=vmin,
                          cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                          label=labels[i+1], xticks=False)

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
