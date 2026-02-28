#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import cooler
import cooltools
import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coolfile", help="Input file (.cool format)", type=str)
    parser.add_argument("chr", help="choromosome", type=str)
    parser.add_argument("outputdir", help="Output directory", type=str)
    parser.add_argument("-r", "--resolution", help="Resolution (bp, default: 100000)", type=int, default=100000)
    parser.add_argument("--raw", help="Output raw matrix (default: balanced)", action='store_true')

    args = parser.parse_args()
#    print(args)

    coolfile = args.coolfile
    chr = args.chr
    resolution = args.resolution
    odir = args.outputdir
    os.makedirs(odir, exist_ok=True)

    if args.raw:
        outputfile = f"{odir}/raw.{chr}.matrix.gz"
    else:
        outputfile = f"{odir}/balanced.{chr}.matrix.gz"

    clr = cooler.Cooler(f"{coolfile}::/resolutions/{resolution}")

    if args.raw:
        array = clr.matrix(balance=False).fetch(chr)
    else:
        array = clr.matrix(balance=True).fetch(chr)

    index = np.arange(array.shape[0]) * resolution
    df = pd.DataFrame(array, index=index, columns=index)
    df.to_csv(outputfile, sep='\t', compression='gzip')


if __name__ == '__main__':
    main()