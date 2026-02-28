#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import cooler
import cooltools.lib.plotting
from cooltools import insulation
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import bioframe
import cooltools
import sys
import argparse
from pathlib import Path
from matplotlib.ticker import EngFormatter
plt.rcParams['font.size'] = 12
import pandas.api.types as ptypes

# Functions to help with plotting
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start+resolution*i for i in range(len(matrix_c)+1)]
    import itertools
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


bp_formatter = EngFormatter('b')
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)


def plot_chr(clr, insul_table, chrom, res, windows, odir):
    """Draw Hi-C + insulation for a single chromosome."""
    chrom_len = clr.chromsizes[chrom]
    region    = (chrom, 0, chrom_len)

    # ── Hi-C matrix ────────────────────────────────────────────
    mat   = clr.matrix(balance=True).fetch(region)
    vmax  = np.nanpercentile(mat, 90)      # dynamic color scale
    norm  = LogNorm(vmin=vmax*1e-3, vmax=vmax)

    fig, ax = plt.subplots(figsize=(20, 10))
    im = pcolormesh_45deg(ax, mat, start=0, resolution=res,
                          norm=norm, cmap='fall')
    ax.set_aspect(0.5)
    # display the top 10 × smallest window (edit as you like)
    ax.set_ylim(0, 10 * windows[0])
    format_ticks(ax, rotate=False)
    ax.xaxis.set_visible(False)

    # ── color-bar ──────────────────────────────────────────────
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="1%", pad=0.1, aspect=6)
    plt.colorbar(im, cax=cax)

    # ── insulation track ──────────────────────────────────────
    ins_ax = div.append_axes("bottom", size="50%", pad=0., sharex=ax)
    ins_chr = bioframe.select(insul_table, region)

    score_col     = f"log2_insulation_score_{windows[0]}"
    bstrength_col = f"boundary_strength_{windows[0]}"
    isbound_col   = f"is_boundary_{windows[0]}"

    ins_ax.plot(ins_chr[['start', 'end']].mean(axis=1),
                ins_chr[score_col],
                label=f"Window {windows[0]} bp")

    # weak / strong boundaries
    boundaries = ins_chr[~np.isnan(ins_chr[bstrength_col])]
    weak   = boundaries[~boundaries[isbound_col]]
    strong = boundaries[boundaries[isbound_col]]

    ins_ax.scatter(weak[['start', 'end']].mean(axis=1),
                   weak[score_col], label="Weak boundaries")
    ins_ax.scatter(strong[['start', 'end']].mean(axis=1),
                   strong[score_col], label="Strong boundaries")

    ins_ax.legend(bbox_to_anchor=(0., -1), loc="lower left", ncol=4)
    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_xlim(0, chrom_len)

    # ── save & clean up ───────────────────────────────────────
    out = f"{odir}/{chrom}_insulation.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_insulation_bedgraphs(ins_table, windows, resolution,
                               metric="log2_insulation_score_{}",
                               outdir="bedgraph"):
    """Save one bedGraph per window size."""
    Path(outdir).mkdir(exist_ok=True)

    for w in windows:
        if w == resolution:
            continue
        col = metric.format(w)
        if col not in ins_table.columns:
            print(f"[warn] column {col} not found – skipping")
            continue

        bg = ins_table[['chrom', 'start', 'end', col]].copy()

        if ptypes.is_numeric_dtype(bg[col]):
            bg[col] = bg[col].fillna(0)
        elif ptypes.is_categorical_dtype(bg[col]):
            bg[col] = bg[col].astype(float).fillna(0)
        else:
            bg[col] = pd.to_numeric(bg[col], errors='coerce').fillna(0)

        fname = Path(outdir) / f"InsulationScore.w{w}.bedgraph"
        with open(fname, "w") as f:
            f.write(f'track type=bedGraph name="Insulation_{w}bp"\n')
            bg.to_csv(f, sep='\t', header=False, index=False)

        print(f"✔  {fname}")



def write_boundary_bed(ins_table, window_bp, out_path,
                       score_col_tpl="boundary_strength_{}",
                       flag_col_tpl="is_boundary_{}"):
    """
    Parameters
    ----------
    ins_table : pd.DataFrame
        cooltools.insulation の出力
    window_bp : int
        例）50000 なら 'boundary_strength_50000' 列を使う
    out_path : str or Path
        生成する BED6 ファイル
    """
    score_col = score_col_tpl.format(window_bp)
    flag_col  = flag_col_tpl.format(window_bp)

    if score_col not in ins_table.columns:
        raise ValueError(f"{score_col} not found.")

    df = ins_table[['chrom', 'start', 'end', flag_col, score_col]].copy()

    # NaN → 0 （score 列のみ）
    df[score_col] = df[score_col].fillna(0)

    # ラベル列を作成：flag=True → strong, False → weak
    df['label'] = np.where(df[flag_col], 'strong', 'weak')
    bed6 = df[['chrom', 'start', 'end', 'label', score_col]].copy()
#    bed6['strand'] = '.'          # 6 列目 (Dummy strand)

    # 出力（header 行は不要／タブ区切り）
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    strong_bed6 = bed6[bed6['label'] == 'strong']
    out_strong = Path(out_path).with_suffix('').as_posix() + ".strong.bed"
    strong_bed6.to_csv(out_strong, sep='\t', header=False, index=False)

    print(f"✔  all boundaries → {out_path}")


def extract_TADs(insulation_table, is_boundary_col, max_TAD_length = 5_000_000):
    tads = bioframe.merge(insulation_table[insulation_table[is_boundary_col] == False])
    return tads[ (tads["end"] - tads["start"]) <= max_TAD_length].reset_index(drop=True)[['chrom','start','end']]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coolfile", help="Input file (.cool format)", type=str)
    parser.add_argument("outputdir", help="Output directory", type=str)
    parser.add_argument("gt", help="Genometable file", type=str)
    parser.add_argument("-r", "--resolution", help="Resolution (bp, default: 25000)", type=int, default=25000)

    args = parser.parse_args()
#    print(args)

    coolfile = args.coolfile
    resolution = args.resolution
    odir = args.outputdir + "/" + str(resolution)
    os.makedirs(odir, exist_ok=True)

    clr = cooler.Cooler(f"{coolfile}::/resolutions/{resolution}")

#    windows = [3*resolution, 5*resolution, 10*resolution, 25*resolution]
    windows = [100000, 500000, 1000000]

    print("Computing insulation …")
    ins_table = insulation(clr, windows, verbose=True)

    write_insulation_bedgraphs(ins_table, windows, resolution, outdir=odir)

    for w in windows:
        if w == resolution:
            print(f"Skipping {w} bp window (same as resolution)")
            continue
        out_bed = f"{odir}/Boundaries.w{w}.bed"
        write_boundary_bed(ins_table, w, out_bed)

        # TADs extraction
        TADs_table = extract_TADs(ins_table, f'is_boundary_{w}')
        out_tads = f"{odir}/TADs.w{w}.bed"
        TADs_table.to_csv(out_tads, sep='\t', header=True, index=False)
        print(f"✔  TADs (w={w}) → {out_tads}  (rows: {len(TADs_table):,})")



    # loop over chromosomes
    print(f"Plotting insulation score..")
    for chrom in clr.chromnames:
        if chrom in ["chrY","chrM","chrMT"]:
            continue
#        print(f"\nDrawing {chrom} …")
        plot_chr(clr, ins_table, chrom, resolution, windows, odir)

    print("\nAll done.")

if __name__ == '__main__':
    main()