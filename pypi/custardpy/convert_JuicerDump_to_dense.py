#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input file (dumped from Juicer)", type=str)
    parser.add_argument("output", help="output file", type=str)
    parser.add_argument("gt", help="genome_table", type=str)
    parser.add_argument("chr", help="chromosome", type=str)
    parser.add_argument("resolution", help="resolution", type=int)
    parser.add_argument("--start", help="start position (default: 0)", type=int, default=0)
    parser.add_argument("--end", help="end position (default: chromosome end)", type=int, default=-1)
    parser.add_argument("--norm", help="normalization type (NONE|VC|VC_SQRT|KR|SCALE, default: SCALE)", type=str, default="SCALE")

    args = parser.parse_args()
    print(args)

    inputfile = args.input
    outputfile = args.output
    gtfile = args.gt
    chr = args.chr
    resolution = args.resolution

    genometable = pd.read_csv(gtfile, delimiter='\t', index_col=[0], header=None)
    chrlen = genometable.loc[chr]
    binlen = int(chrlen/resolution) +1

    d = pd.read_csv(inputfile, delimiter='\t', header=None)
    d = d.set_index([0,1])
    d = d.iloc[:,0]
    d = d.unstack(fill_value=0)
    index = np.arange(binlen) * resolution
    d = d.reindex(index, columns=index, fill_value=0)
    d.index.name = None
    d.columns.name = None

    triu = np.triu(d)
    array = triu + triu.T - np.diag(np.diag(triu))
    df = pd.DataFrame(array, index=d.index, columns=d.columns)

    df.to_csv(outputfile, sep='\t', compression='gzip')
