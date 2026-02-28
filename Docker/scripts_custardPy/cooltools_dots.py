#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import argparse
from itertools import chain

import cooler
import bioframe
import cooltools
from cooltools.lib.numutils import fill_diag

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter

# helper functions for plotting
bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    """format ticks with genomic coordinates as human readable"""
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

# create a functions that would return a series of rectangles around called dots
# in a specific region, and exposing importnat plotting parameters
def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="cyan", fc="none"):
    """
    yield a series of rectangles around called dots in a given region
    """
    # select dots from the region:
    df_reg = bioframe.select(
        bioframe.select(dots_df, region, cols=("chrom1","start1","end1")),
        region,
        cols=("chrom2","start2","end2"),
    )
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        width1 = e1 - s1
        width2 = e2 - s2
        if loc == "upper":
            yield patches.Rectangle((s2, s1), width2, width1, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((s1, s2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")

def plot_exampleregion(dots_df, clr, odir):
    # define a region to look into as an example
    start = 34_150_000
    end = start + 1_200_000
    region = ('chr17', start, end)

    # heatmap kwargs
    matshow_kwargs = dict(
        cmap='YlOrBr',
        norm=LogNorm(vmax=0.05),
        extent=(start, end, end, start)
    )

    # colorbar kwargs
    colorbar_kwargs = dict(fraction=0.046, label='corrected frequencies')

    # compute heatmap for the region
    region_matrix = clr.matrix(balance=True).fetch(region)
    for diag in [-1,0,1]:
        region_matrix = fill_diag(region_matrix, np.nan, i=diag)

    # see viz.ipynb for details of heatmap visualization
    f, ax = plt.subplots(figsize=(7,7))
    im = ax.matshow(region_matrix, **matshow_kwargs)
    format_ticks(ax, rotate=False)
    plt.colorbar(im, ax=ax, **colorbar_kwargs)

    # draw rectangular "boxes" around pixels called as dots in the "region":
    for box in rectangles_around_dots(dots_df, region, lw=1.5):
        ax.add_patch(box)
    plt.savefig(f"{odir}/exampleregion.png", dpi=300)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coolfile", help="Input file (.cool format)", type=str)
    parser.add_argument("outputdir", help="Output directory", type=str)
    parser.add_argument("gt", help="Genometable file", type=str)
    parser.add_argument("-r", "--resolution", help="Resolution (bp, default: 10000)", type=int, default=10000)
    parser.add_argument("--maxseparation", help="How far from the main diagonal to call dots (default: 10000000)", type=int, default=10000000)
    parser.add_argument("-p", "--threads", help="Number of CPUs (default: 4)", type=int, default=4)

    args = parser.parse_args()
#    print(args)

    coolfile = args.coolfile
    resolution = args.resolution
    odir = args.outputdir + "/loops/cooltools"
    gt = args.gt
    nproc = args.threads
    maxseparation = args.maxseparation
    os.makedirs(odir, exist_ok=True)

#    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
#    hg38_cens = bioframe.fetch_centromeres('hg38')
#    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
    # Select only chromosomes that are present in the cooler.
#    hg38_arms = hg38_arms.set_index("chrom").loc[clr.chromnames].reset_index()

    clr = cooler.Cooler(f"{coolfile}::/resolutions/{resolution}")

    # intra-arm expected
    expected = cooltools.expected_cis(
        clr,
#        view_df=hg38_arms,
        nproc=24,
    )

    dots_df = cooltools.dots(
        clr,
        expected=expected,
#        view_df=hg38_arms,
        # how far from the main diagonal to call dots:
        max_loci_separation=maxseparation,
        nproc=nproc,
    )

    plot_exampleregion(dots_df, clr, odir)

    dots_df.to_csv(f"{odir}/dots.{resolution}.tsv", sep="\t", index=False)

    print("\nAll done.")

if __name__ == '__main__':
    main()

