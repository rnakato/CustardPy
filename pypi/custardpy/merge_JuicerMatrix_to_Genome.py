#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys

def parse_argv():
    usage = 'Usage: \n    python {} <matrixdir> <outputfilename> <resolution> <observed/oe> <lim_pzero> <num of chr> [--help] [--includeintra] [--evenodd]'.format(__file__)
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

def getfilename(i, j):
    return dir + "/" + str(res) + "/chr" + str(i) + "-chr" + str(j) + "/" + ntype + ".matrix.gz"

def getchrlen(nchr):
    chrlen = []
    for i in range(1,nchr+1):
        d = pd.read_csv(getfilename(1, i), delimiter='\t', header=None)
        del d[d.shape[1]-1]
        chrlen.append(d.shape[1])

    #print(chrlen)
    #print(np.sum(chrlen))

    return chrlen

def getmatrix(i, j, chrlen, include_intra_read):
    def rename(d):
        d.index = map(str, d.index * res)
        d.index = d.index.map(lambda x: "chr" + str(i) + "-" + x)
        d.columns = map(str, d.columns * res)
        d.columns = d.columns.map(lambda x: "chr" + str(j) + "-" + x)

    if j < i:
        d = pd.read_csv(getfilename(j, i), delimiter='\t', header=None)
        del d[d.shape[1]-1]
        d = d.T
        rename(d)
        return d
    elif i==j and include_intra_read==False:
        d = pd.DataFrame(np.zeros((chrlen[i-1], chrlen[j-1])))
        rename(d)
        return d
    else:
        d = pd.read_csv(getfilename(i, j), delimiter='\t', header=None)
        del d[d.shape[1]-1]
        rename(d)
        return d

if __name__ == '__main__':
    arguments = parse_argv()
    dir = arguments[0]
    outputfile = arguments[1]
    res = int(arguments[2])
    ntype = arguments[3]
    lim_pzero = float(arguments[4])
    nchr = int(arguments[5])

    include_intra_read = False
    if '--includeintra' in arguments:
        include_intra_read = True

    chrlen = getchrlen(nchr)

    if '--evenodd' in arguments:
        for i in range(1,nchr+1,2):
            #            print('i={0} j=2'.format(i))
            matrix = getmatrix(i, 2, chrlen, include_intra_read)
            for j in range(4,nchr+1,2):
 #               print('i={0} j={1}'.format(i,j))
                mat = getmatrix(i, j, chrlen, include_intra_read)
                matrix = pd.concat([matrix, mat], axis=1)
            if i==1:
                A = matrix
            else:
                A = pd.concat([A, matrix])

        print("before trim: ")
        print(A.shape)

        pzero_raw = A[A>0].count(axis=0)/A.shape[0]
        index = pzero_raw[pzero_raw > lim_pzero].index
        pzero_col = A[A>0].count(axis=1)/A.shape[1]
        columns = pzero_col[pzero_col > lim_pzero].index

        A = A[index]
        A = A.loc[columns]

        print("after trim: ")
        print(A.shape)

    else:
        for i in range(1,nchr+1):
            matrix = getmatrix(i, 1, chrlen, include_intra_read)
            for j in range(2,nchr+1):
                mat = getmatrix(i, j, chrlen, include_intra_read)
                matrix = pd.concat([matrix, mat], axis=1)
            if i==1:
                A = matrix
            else:
                A = pd.concat([A, matrix])

        print("matrix size: ")
        print(A.shape)

    A.to_pickle(outputfile)
