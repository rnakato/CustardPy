#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

def fixEigendir(filename, refFlat, chr, res):
   eigen = np.loadtxt(filename)
   gene = pd.read_csv(refFlat, delimiter='\t', header=None, usecols=[2,3,4,5])
   gene = gene[gene.iloc[:,0] == chr]
   gene = gene.iloc[:,2].values
   genenum = np.zeros(len(eigen), int)
   for row in gene:
      row = int(row)
      if int(row/res)>= len(eigen): continue
      genenum[int(row/res)] += 1
   if pearsonr(genenum[~np.isnan(eigen)], eigen[~np.isnan(eigen)])[0] < 0:
      eigen = -eigen

   return eigen

if(__name__ == '__main__'):
    parser = argparse.ArgumentParser()
    parser.add_argument("input",  help="Eigen data",    type=str)
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("gene",   help="refFlat file",  type=str)
    parser.add_argument("chr",    help="chromosome",    type=str)
    parser.add_argument("resolution", help="resolution", type=int)

    args = parser.parse_args()
#    print(args)

    eigen = fixEigendir(args.input,
                        args.gene,
                        args.chr,
                        args.resolution)
    np.savetxt(args.output, eigen, fmt="%0.6f")
