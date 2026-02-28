# Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
# All rights reserved.
 
import os
import pandas as pd

def loadDenseMatrix(filename):
    print(filename)
    data = pd.read_csv(filename, delimiter='\t', index_col=0)
    data.columns = data.columns.astype('int')
    return data

def loadTADs(filename, chr, *, start=0, end=99999999999):
    if os.path.exists(filename):
#        tads = pd.read_csv(filename, delimiter='\t', usecols=['#chr1','x1','x2'])
        tads = pd.read_csv(filename, header=0, delimiter='\t', usecols=[0,1,2], names=['chr','x1','x2'])
#        tads.rename(columns={"#chr1":"chr"},inplace =True)
        tads = tads[(tads.chr == chr) & (tads.x1 < end) & (tads.x2 >= start)]
    else:
        print ("Warning: " + filename + " is not available. Skipping")
        tads = None
    return tads

def loadloops(filename, chr, *, start=0, end=99999999999):
    if os.path.exists(filename):
#        loops = pd.read_csv(filename, delimiter='\t', usecols=['#chr1','x1','x2','chr2','y1','y2'])
        loops = pd.read_csv(filename, header=0, delimiter='\t', usecols=[0,1,2,3,4,5], names=['chr1','x1','x2','chr2','y1','y2'], comment="#")
#        loops.drop([0], inplace=True)
#        loops = loops[loops.chr1 == chr]
#        loops = loops[loops.chr2 == chr]
#        loops = loops[loops.x2 < end]
#        loops = loops[loops.y1 >= start]
        loops = loops[(loops.chr1 == chr) & (loops.chr2 == chr) & (loops.x2 < end) & (loops.y1 >= start)]
    else:
        print ("Warning: " + filename + " is not available. Skipping")
        loops = None
    return loops

def readBedGraph(file, chr, *, start=-1, end=-1):
    bedgraph = pd.read_csv(file, delimiter='\s+', header=None)
    bedgraph.columns = ["chr", "start", "end", "value"]
    bedgraph = bedgraph[bedgraph["chr"] == chr]
    bedgraph = bedgraph.set_index("start")
    bedgraph = bedgraph["value"]

    if start >= 0 and end > 0:
        bedgraph = bedgraph[bedgraph.index >= start]
        bedgraph = bedgraph[bedgraph.index <= end]
    return bedgraph
