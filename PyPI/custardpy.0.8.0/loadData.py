# Copyright(c) Ryuichiro Nakato <rnakato@iqb.u-tokyo.ac.jp>
# All rights reserved.
 
import os
import numpy as np
import pandas as pd


def _load_from_hic(filename, chrom, res):
    import hicstraw
    records = hicstraw.straw("NONE", filename, chrom, chrom, "BP", res)
    if not records:
        raise ValueError(f"No data found for {chrom} at resolution {res} in {filename}")
    positions = sorted(set([r.binX for r in records] + [r.binY for r in records]))
    pos_idx = {p: i for i, p in enumerate(positions)}
    n = len(positions)
    mat = np.zeros((n, n))
    for r in records:
        i, j = pos_idx[r.binX], pos_idx[r.binY]
        mat[i, j] = r.counts
        mat[j, i] = r.counts
    return pd.DataFrame(mat, index=positions, columns=positions)


def _load_from_cool(filename, chrom, res):
    import cooler
    if filename.endswith('.mcool'):
        c = cooler.Cooler(f"{filename}::/resolutions/{res}")
    else:
        c = cooler.Cooler(filename)
    mat = c.matrix(balance=False).fetch(chrom)
    bins = c.bins().fetch(chrom)['start'].values
    return pd.DataFrame(mat, index=bins, columns=bins)


def loadDenseMatrix(filename, chrom=None, res=None):
    print(filename)
    if filename.endswith('.hic'):
        if chrom is None or res is None:
            raise ValueError("chrom and res are required for .hic files")
        return _load_from_hic(filename, chrom, res)
    elif filename.endswith('.cool') or filename.endswith('.mcool'):
        if chrom is None or res is None:
            raise ValueError("chrom and res are required for .cool/.mcool files")
        return _load_from_cool(filename, chrom, res)
    else:
        data = pd.read_csv(filename, delimiter='\t', index_col=0)
        data.columns = data.columns.astype('int')
        return data

def loadTADs(filename, chr, *, start=0, end=99999999999):
    if os.path.exists(filename):
        tads = pd.read_csv(filename, header=0, delimiter='\t', usecols=[0,1,2], names=['chr','x1','x2'])
        tads = tads[(tads.chr == chr) & (tads.x1 < end) & (tads.x2 >= start)]
    else:
        print ("Warning: " + filename + " is not available. Skipping")
        tads = None
    return tads

def loadloops(filename, chr, *, start=0, end=99999999999):
    if os.path.exists(filename):
        loops = pd.read_csv(filename, header=0, delimiter='\t', usecols=[0,1,2,3,4,5], names=['chr1','x1','x2','chr2','y1','y2'], comment="#")
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
