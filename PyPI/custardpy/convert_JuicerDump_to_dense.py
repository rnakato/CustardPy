#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys
import argparse

def parse_argv():
    usage = 'Usage: \n    python {} <inputfile> <outputfile> <genometable> <chr1> <chr2> <resolution> [--help]'.format(__file__)
    arguments = sys.argv
    if len(arguments) == 1:
        print(usage)
        exit()
    # ファイル自身を指す最初の引数を除去
    arguments.pop(0)
    # - で始まるoption
    options = [option for option in arguments if option.startswith('-')]

    if '-h' in options or '--help' in options:
        print(usage)

    return arguments

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("listfile", help="Input file (output of Juicer dump)", type=str)
    parser.add_argument("output", help="Output file name", type=str)
    parser.add_argument("genometable", help="genome_table file", type=str)
    parser.add_argument("chr1", help="1st chromosome", type=str)
    parser.add_argument("chr2", help="1st chromosome", type=str)
    parser.add_argument("-r", "--resolution", help="Resolution (bp) of the input matrix (default: 1000000)", type=int, default=1000000)

    args = parser.parse_args()
#    print(args)

    inputfile = args.listfile
    outputfile = args.output
    gtfile = args.genometable
    chr1 = args.chr1
    chr2 = args.chr2
    resolution = args.resolution

    genometable = pd.read_csv(gtfile, delimiter='\t', index_col=[0], header=None)
    chrlen1 = int(genometable.loc[chr1])
    binlen1 = int(chrlen1/resolution) +1
    chrlen2 = int(genometable.loc[chr2])
    binlen2 = int(chrlen2/resolution) +1

    d = pd.read_csv(inputfile, delimiter='\t', header=None)
    d = d.set_index([0,1])
    d = d.iloc[:,0]
    d = d.unstack(fill_value=0)

    index = np.arange(binlen1) * resolution
    columns = np.arange(binlen2) * resolution
    d = d.reindex(index, columns=columns, fill_value=0)
    d.index.name = None
    d.columns.name = None

    if chr1 == chr2:
        triu = np.triu(d)
        array = triu + triu.T - np.diag(np.diag(triu))
        df = pd.DataFrame(array, index=d.index, columns=d.columns)
        df.to_csv(outputfile, sep='\t', compression='gzip')
    else:
        d.to_csv(outputfile, sep='\t', compression='gzip')

#    index = np.arange(binlen1) * resolution
#    columns = np.arange(binlen2) * resolution
#    d = d.reindex(index, columns=columns, fill_value=0)
#    d.index.name = None
#    d.columns.name = None

#    triu = np.triu(d)
#    array = triu + triu.T - np.diag(np.diag(triu))
#    df = pd.DataFrame(array, index=d.index, columns=d.columns)

#    df.to_csv(outputfile, sep='\t', compression='gzip')
