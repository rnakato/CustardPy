#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import argparse
from custardpy.loadData import loadDenseMatrix
from custardpy.InsulationScore import MultiInsulationScore

def calcDI(mat, resolution, *, distance=1000000):
    def getDI(mat, i, len):
        A = round(np.triu(mat[i-len:i, i-len:i]).sum(), 6)
        B = round(np.triu(mat[i:i+len, i:i+len]).sum(), 6)
        E = round((A + B)/2, 6)

        if E != 0:
            DI = ((A-E)**2)*2 / E
        else:
            DI = 0

        if A > B:
            DI = -DI
        return DI

    len = int(distance / resolution)
    array = np.zeros(mat.shape[0])
    mat2 = np.nan_to_num(mat)
    for i in range(len, mat2.shape[0]-len):
        array[i] = getDI(mat2, i, len)
    return array

def getDirectionalityIndexOfMultiSample(samples, labels, *, distance=1000000):
    if not isinstance(samples, list):
        samples = [samples]

    for i, sample in enumerate(samples):
        if i==0: Matrix = sample.getDirectionalityIndex(distance=distance)
        else:    Matrix = np.vstack((Matrix,sample.getDirectionalityIndex(distance=distance)))
    Matrix = np.nan_to_num(Matrix)
    return Matrix

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", help="Input matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="Chromosome", type=str)
    parser.add_argument("resolution", help="Resolution of the input matrix", type=int)
    parser.add_argument("--num4norm", help="Read number after normalization (default: 10000000)", type=int, default=10000000)
    parser.add_argument("--distance", help="Distance of Insulation Score (default: 500000)", type=int, default=500000)

    args = parser.parse_args()
    print(args)

    matrix = loadDenseMatrix(args.matrix)
    matrix = matrix * args.num4norm / np.nansum(matrix)

    MI = MultiInsulationScore(matrix.values, 1000000, 100000, args.resolution)

    # output InsulationScore to BedGraph
    df = MI.getInsulationScore(distance=args.distance)
    df = df.replace([np.inf, -np.inf], 0)
    df.columns = ["Insulation Score"]
    df["chr"] = args.chr
    df["start"] = df.index
    df["end"] = df["start"] + args.resolution
    df = df.loc[:,["chr","start","end","Insulation Score"]]
    df.to_csv(args.output + ".bedGraph", sep="\t", header=False, index=False)
