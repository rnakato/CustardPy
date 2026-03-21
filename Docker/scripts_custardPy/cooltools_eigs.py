#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, subprocess
import cooler
import cooltools.lib.plotting
import bioframe
import cooltools
import sys
from functools import partial

from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cytoolz import merge

def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    groupmean = track[track.columns[3]].groupby(
        digitized_track[digitized_track.columns[3]]
    ).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)
    else:
        lo, hi = binedges[0], binedges[-1]

    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = np.asarray(saddledata, dtype=float)

    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # log表示では non-positive を mask
    if scale == "log":
        C = np.ma.masked_invalid(C)
        C = np.ma.masked_less_equal(C, 0)

        valid = np.asarray(C.compressed(), dtype=float)
        if valid.size == 0:
            raise ValueError("Saddle matrix has no positive finite values for log scale.")

        if vmin is None or (not np.isfinite(vmin)) or vmin <= 0:
            vmin = np.nanpercentile(valid, 5)
        if vmax is None or (not np.isfinite(vmax)) or vmax <= 0:
            vmax = np.nanpercentile(valid, 95)
        if vmin >= vmax:
            vmin = max(np.min(valid), vmax * 1e-3)

        norm = LogNorm(vmin=vmin, vmax=vmax)

    elif scale == "linear":
        C = np.ma.masked_invalid(C)
        if vmin is None:
            vmin = np.nanmin(C)
        if vmax is None:
            vmax = np.nanmax(C)
        norm = Normalize(vmin=vmin, vmax=vmax)

    else:
        raise ValueError("Only linear and log color scaling is supported")

    if subplot_spec is not None:
        GridSpecLocal = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    else:
        GridSpecLocal = GridSpec

    gs = GridSpecLocal(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    grid = {}
    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap=cmap, rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})

    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])
    plt.barh(
        binedges, height=1 / len(binedges), width=groupmean, align="edge", **margin_kws
    )
    plt.xlim(plt.xlim()[1], plt.xlim()[0])
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)

    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])
    plt.bar(
        binedges, width=1 / len(binedges), height=groupmean, align="edge", **margin_kws
    )
    plt.xlim(lo, hi)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})

    if scale == "linear":
        cb = plt.colorbar(img, cax=grid["ax_cbar"], **cbar_kws)
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        cb = plt.colorbar(
            img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws
        )
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid["ax_heatmap"].grid(False)

    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid


def main():
    coolfile = sys.argv[1]
    cis_eigs_file = sys.argv[2]

    odir = os.path.dirname(cis_eigs_file)

    clr = cooler.Cooler(coolfile)
    matrix_balanced = clr.matrix(balance=True)  # re-use inside loop

    view_df = pd.DataFrame({'chrom': clr.chromnames,
                            'start': 0,
                            'end': clr.chromsizes.values,
                            'name': clr.chromnames}
                        )

    cis_eigs = pd.read_csv(cis_eigs_file, sep='\t', header=0, index_col=None)
    eigenvector_track = cis_eigs[['chrom','start','end','E1']]

    # ── loop over chromosomes ──────────────────────────────────────
    for chrom in clr.chromnames:
        if chrom in ["chrY", "chrM", "chrMT"]:
            continue

        mat_chr = matrix_balanced.fetch(chrom)
        eig_chr = eigenvector_track.query("chrom == @chrom")["E1"].values

        # finite & positive only for LogNorm
        valid = mat_chr[np.isfinite(mat_chr) & (mat_chr > 0)]
        if valid.size == 0:
            print(f"[warn] {chrom}: no positive finite balanced contacts, skipping.")
            continue

        vmax = np.nanpercentile(valid, 99)
        vmin = np.nanpercentile(valid, 5)

        if not np.isfinite(vmax) or vmax <= 0:
            print(f"[warn] {chrom}: invalid vmax={vmax}, skipping.")
            continue
        if not np.isfinite(vmin) or vmin <= 0 or vmin >= vmax:
            vmin = vmax * 1e-3

        norm = LogNorm(vmin=vmin, vmax=vmax)

        mat_plot = np.ma.masked_invalid(mat_chr)
        mat_plot = np.ma.masked_less_equal(mat_plot, 0)

        fig, ax_hic = plt.subplots(figsize=(7, 6))
        im = ax_hic.matshow(mat_plot, norm=norm, cmap="fall")
        ax_hic.set_title(f"{chrom} Hi-C map", loc="left", fontsize=10)
        ax_hic.set_xlabel(f"{chrom} (bins)")
        ax_hic.set_ylabel(f"{chrom} (bins)")

        divider = make_axes_locatable(ax_hic)
        cax = divider.append_axes("right", size="3%", pad=0.1)
        plt.colorbar(im, cax=cax, label="balanced contacts")

        ax_pc1 = divider.append_axes("top", size="20%", pad=0.25, sharex=ax_hic)
        ax_pc1.plot(eig_chr, lw=0.5)
        ax_pc1.set_ylabel("PC1")
        ax_pc1.set_xticks([])

        fig.tight_layout()
        fig.savefig(f"{odir}/Contactmap_PC1.{chrom}.pdf", dpi=300)
        plt.close(fig)

    print("Done: individual PDFs written for each chromosome.")

    ### Saddle plot
    print("Generating saddle plot...")
    Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
    Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
    N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each

    cvd = cooltools.expected_cis(
            clr=clr,
            view_df=view_df,
    )

    interaction_sum, interaction_count = cooltools.saddle(
            clr,
            cvd,
            eigenvector_track,
            'cis',
            n_bins=N_GROUPS,
            qrange=(Q_LO,Q_HI),
            view_df=view_df
    )

    saddleplot(eigenvector_track,
            interaction_sum/interaction_count,
            N_GROUPS,
            qrange=(Q_LO,Q_HI),
            cbar_kws={'label':'average observed/expected contact frequency'}
            );
    plt.savefig(f"{odir}/saddleplot.pdf", dpi=300)
    plt.close() 

    print("Done: saddle plot saved.")
    
    ### Saddle strength profile
    print("Generating saddle strength profile...")
    from cooltools.api.saddle import saddle_strength

    x = np.arange(N_GROUPS + 2)

    plt.step(x, saddle_strength(interaction_sum, interaction_count), where='pre')

    plt.xlabel('extent')
    plt.ylabel('(AA + BB) / (AB + BA)')
    plt.title('saddle strength profile')
    plt.axhline(1, c='grey', ls='--', lw=1) 
    plt.xlim(0, 1);
    plt.savefig(f"{odir}/saddle_strength.pdf", dpi=300)
    print("Done: saddle strength profile saved.")


if __name__ == '__main__':
    main()