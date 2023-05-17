#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy import ndimage
import argparse

def make3dmatrixRatio(samples, smoooth=3):
    n = len(samples)
    Ct = ndimage.median_filter(samples[0].getlog(isNonZero=False), smoooth)
    x, y = Ct.shape
    for i, sample in enumerate(samples[1:]):
        if i==0:
            data = sample.getlog(isNonZero=False)
            Matrix = ndimage.median_filter(data - Ct, smoooth)
        else:
            data = sample.getlog(isNonZero=False)
            M = ndimage.median_filter(data - Ct, smoooth)
            Matrix = np.concatenate((Matrix, M))
    Matrix = Matrix.reshape(n-1,x,y)
    return Matrix

def getDirectionalRelativeFreq(mat, resolution, strand, *,
                               startdistance=0, distance=2000000):
    if (startdistance >= distance):
        print ("getDirectionalRelativeFreq: Error: startdistance > enddistance")
        exit(1)

    arraysize = mat.shape[0]
    array = np.zeros(arraysize)
    nbin = int(distance/resolution)
    startbin = int(startdistance/resolution) +1
    for i in range(nbin, arraysize - nbin):
        if (strand == "+"):
            val = mat[i+startbin:i+nbin+1, i].mean()
        else:
            val = mat[i, i-nbin:i-startbin+1].mean()
        array[i] = val

    return array

class DirectionalRelativeFreq:
    def __init__(self, mat, resolution, *, startdistance=0, distance=2000000):
        self.arrayplus  = getDirectionalRelativeFreq(mat, resolution, "+", startdistance=startdistance, distance=distance)
        self.arrayminus = getDirectionalRelativeFreq(mat, resolution, "-", startdistance=startdistance, distance=distance)

    def getarrayplus(self):
        return self.arrayplus

    def getarrayminus(self):
        return self.arrayminus

    def getarraydiff(self):
        return self.arrayplus - self.arrayminus

def output_DRF(args):
    from custardpy.HiCmodule import JuicerMatrix
    resolution = args.resolution
    samples = []
    samples.append(JuicerMatrix("RPM", args.control, resolution))
    samples.append(JuicerMatrix("RPM", args.input, resolution))

    smooth_median_filter = 3
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)

#    import pdb; pdb.set_trace()

    drf = DirectionalRelativeFreq(EnrichMatrices[0], resolution)
    if (args.drf_right == True):
        array = drf.getarrayplus()
    elif (args.drf_left == True):
        array = drf.getarrayminus()
    else:
        array = drf.getarraydiff()

    df = pd.DataFrame(array)
    df.columns = ["DRF"]
    df["chr"] = args.chr
    df["start"] = np.arange(len(array)) * resolution
    df["end"] = df["start"] + resolution
    df = df.loc[:,["chr","start","end","DRF"]]

    df.to_csv(args.output + ".bedGraph", sep="\t", header=False, index=False)
#    np.savetxt(args.output, array, fmt="%0.6f")

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input matrix", type=str)
    parser.add_argument("control", help="Control matrix", type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("resolution", help="Resolution of the input matrix", type=int)
    parser.add_argument("--drf_right",   help="(with --drf) plot DirectionalRelativeFreq (Right)", action='store_true')
    parser.add_argument("--drf_left",   help="(with --drf) plot DirectionalRelativeFreq (Left)", action='store_true')

    args = parser.parse_args()
    print(args)

    output_DRF(args)
