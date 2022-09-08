#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from custardpy.HiCmodule import JuicerMatrix
from custardpy.generateCmap import *
from custardpy.PlotModule import *
from custardpy.loadData import *

#import pdb

def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input", help="<Input direcoty>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("--heatmap", help="heatmap", action='store_true')

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

    resolution = args.resolution
    type = args.type

    samples = []
    nchr = 22
    length = 0
    plt.figure(figsize=(18, nchr))
    for chr in range(1, nchr):
        plt.subplot2grid((nchr, 4), ((chr-1), 0), rowspan=1, colspan=4)
        for i, dir in enumerate(dirs):
            eigenfile = dir + "/eigen/" + str(resolution) + "/gd_eigen."  + type + ".chr" + str(chr) + ".txt"
            print (eigenfile)
            eigen = np.loadtxt(eigenfile)
            if(length==0): length = len(eigen)

            if i==0: Matrix = eigen
            else:    Matrix = np.vstack((Matrix, eigen))

        if(args.heatmap==True):
            plt.imshow(Matrix, clim=(-0.05, 0.05),
                       cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                       aspect="auto")
            plt.colorbar()
            plt.yticks(np.arange(len(labels)), labels)
            plt.xlim(0, length)
        else:
            plt.plot(Matrix.T)
            plt.ylim(-0.05, 0.05)
            plt.xlim(0, length)
            plt.title("chr" + str(chr))
            #            plt.legend()

    plt.tight_layout()
    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
