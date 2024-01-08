#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

usage = 'Usage: \n    {} <QC dir>'.format(os.path.basename(__file__))
arguments = sys.argv
if len(arguments) == 1:
    print(usage)
    exit()

qcdir   = sys.argv[1]
qcfile  = qcdir + "/output/scores/qc.genomewide.txt"
repfile = qcdir + "/output/scores/reproducibility.genomewide.txt"
figdir  = qcdir + "/pdf/"

if not os.path.isdir(figdir):
    os.makedirs(figdir)

df = pd.read_csv(qcfile, sep='\t', index_col=[0], header=None, skiprows=1)
df.columns = ["QuASAR-QC"]
#print (df)
df.plot(kind='barh', figsize=(5, 10))
plt.savefig(figdir + "QuASAR-QC.pdf")

df = pd.read_csv(repfile, sep='\t', index_col=['#Sample1','Sample2'])
fig = sns.clustermap(df["GenomeDISCO"].unstack(), figsize=(12,12))
fig.savefig(figdir + "GenomeDISCO.cluster.pdf")

fig = sns.clustermap(df["HiC-Spector"].unstack(), figsize=(12,12))
fig.savefig(figdir + "HiC-Spector.cluster.pdf")

fig = sns.clustermap(df["HiCRep"].unstack(), figsize=(12,12))
fig.savefig(figdir + "HiCRep.cluster.pdf")

print("The results are in " + figdir + ".")
