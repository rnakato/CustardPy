#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
import argparse

def get_args():
    # 準備
    usage = 'Usage: python {} INPUT OUTPUT [--start1 <int>] [--end1 <int>] [--start2 <int>] [--end2 <int>] [--help]'\
            .format(__file__)
    parser = argparse.ArgumentParser(usage=usage)

    parser.add_argument("input",  type=str)
    parser.add_argument("output", type=str)

    parser.add_argument("--start1", default=-1, type=int)
    parser.add_argument("--end1",   default=-1, type=int)
    parser.add_argument("--start2", default=-1, type=int)
    parser.add_argument("--end2",   default=-1, type=int)

    # 結果を受ける
    args = parser.parse_args()

    return(args)

def main():
    args = get_args()

    data = pd.read_csv(args.input, delimiter='\t', index_col=0)
    if args.start1 > 0: s1 = args.start1
    else: s1 = 0
    if args.end1 > 0: e1 = args.end1
    else: e1 = data.shape[0]
    if args.start2 > 0: s2 = args.start2
    else: s2 = 0
    if args.end2 > 0: e2 = args.end2
    else: e2 = data.shape[1]

    matrix = data.values[s1:e1,s2:e2]
    np.save(args.output, matrix)

if __name__ == '__main__':
    main()
