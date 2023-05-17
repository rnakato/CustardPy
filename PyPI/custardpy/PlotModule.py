import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage
from custardpy.generateCmap import *

Mega = 1000000

def pltxticks_mega_subplot2grid(start, end, figstart, figend, *, ax=None):
    xlen = (figend - figstart)
    nmem = int(xlen/Mega) + 1
    x = start + np.arange(nmem) * (end-start) * Mega/xlen
    xval = np.array(((figstart/Mega + np.arange(nmem))*10), dtype=np.int64) / 10

    if ax is None:
        plt.xticks(x, xval)
    else:
        ax.set_xticks(x, xval)

def pltxticks_subplot2grid(start, end, figstart, figend, nmem, *, ax=None):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = np.array(((figstart + np.arange(nmem+1) * val) / Mega * 10), dtype=np.int64) / 10
    
#    print (x)
#    print (xval)
    if ax is None:
        plt.xticks(x, xval)
    else:
        ax.set_xticks(x, xval)

def pltyticks_subplot2grid(start, end, figstart, figend, nmem, *, ax=None):
    mem = int((end - start)/nmem)
    x = start + np.arange(nmem+1) * mem
    val = (figend - figstart) / nmem
    xval = (figstart + np.arange(nmem+1)*val)/Mega

    if ax is None:
        plt.yticks(x, xval)
    else:
        ax.set_yticks(x, xval)

def ytickoff_ax(*, ax=None):
    if ax is None:
        plt.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left=False,        # ticks along the bottom edge are off
            right=False,       # ticks along the top edge are off
            labelleft=False
        )
    else:
        ax.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left=False,        # ticks along the bottom edge are off
            right=False,       # ticks along the top edge are off
            labelleft=False
        )

def xtickoff_ax(*, ax=None):
    if ax is None:
        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False  # labels along the bottom edge are off
        )
    else:
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False  # labels along the bottom edge are off
        )

def calculate_start_end(matrix, figstart, figend, resolution):
    s = int(figstart / resolution)
    e = int(figend / resolution)
    if e == 0:
        e = matrix.shape[0]
    return s, e

def extract_submatrix(matrix, s, e):
    if isinstance(matrix, pd.DataFrame):
        return matrix.iloc[s:e, s:e]
    else:
        return matrix[s:e, s:e]

def drawHeatmapSquare(matrix, resolution, *, tads=None, loops=None,
                      figstart=0, figend=0, distancemax=0,
                      vmin=0, vmax=50, label="", xticks=True,
                      cmap=generate_cmap(['#FFFFFF', '#d10a3f'])):

    s, e = calculate_start_end(matrix, figstart, figend, resolution)
    mat = extract_submatrix(matrix, s, e)

    plt.imshow(mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect=1)

    if not tads is None:
        for tad in tads.itertuples(name=None):
            x1 = int(tad[2])/resolution - s
            x2 = int(tad[3])/resolution - s
            if x1 > 0:
                plt.vlines(x1, x1, min(x2, e-s), linestyle='dashed', linewidth=0.6)
            if x2 < e-s:
                plt.vlines(x2, max(0, x1), x2, linestyle='dashed', linewidth=0.6)
            if x1 > 0:
                plt.hlines(x1, x1, min(x2, e-s), linestyle='dashed', linewidth=0.6)
            if x2 < e-s:
                plt.hlines(x2, max(0, x1), x2, linestyle='dashed', linewidth=0.6)

    plt.xlim(0, e-s)
    plt.ylim(e-s, 0)

    plt.title(label)
    if (xticks):
        pltxticks_subplot2grid(0, e-s, figstart, figend, 5)
    else:
        xtickoff_ax()

    ytickoff_ax()

def get_xbin_sqrt2(posi, s, resolution):
    return (posi/resolution - s) * np.sqrt(2)

def calculate_fig_length(e, s):
    sqrt2 = np.sqrt(2)
    return (e - s) * sqrt2

def limit_yaxis_subplot2grid(distance_max, resolution, ynum, _figstart, _figend, *, ax=None):
    if (distance_max > 0):
        ystart = max(0, ynum - distance_max/resolution)
        figstart = 0
        figend = min(distance_max, (_figend - _figstart))
    else:
        ystart = 0
        figstart = 0
        figend = _figend - _figstart

    if ax is None:
        plt.ylim(ynum, ystart)
    else:
        ax.set_ylim(ynum, ystart)
    pltyticks_subplot2grid(ynum, ystart, figstart, figend, 5, ax=ax)

def plot_tads_subplot2grid(tads, s, fig_length, resolution, ynum, *, ax=None):
    if tads is None:
        return

    for tad in tads.itertuples(name=None):
        x1 = get_xbin_sqrt2(tad[2], s, resolution) #(tad[2]/resolution - s) * np.sqrt(2)
        x2 = get_xbin_sqrt2(tad[3], s, resolution) #(tad[3]/resolution - s) * np.sqrt(2)
        if (x1 >0):
            xmed = (x1 + min([x2, fig_length]))/2
            if ax is None:
                plt.plot([x1, min([xmed, x2])],[ynum, ynum-(xmed-x1)],
                        color='k', linestyle='dashed', linewidth=0.6)
            else:
                ax.plot([x1, min([xmed, x2])],[ynum, ynum-(xmed-x1)],
                        color='k', linestyle='dashed', linewidth=0.6)

        if (x2 < fig_length):
            xmed = (max([x1, 0]) + x2) /2
            if ax is None:
                plt.plot([x2, xmed],[ynum, ynum-(x2-xmed)],
                        color='k', linestyle='dashed', linewidth=0.6)
            else:
                ax.plot([x2, xmed],[ynum, ynum-(x2-xmed)],
                        color='k', linestyle='dashed', linewidth=0.6)

def plot_loops_subplot2grid(loops, s, fig_length, resolution, ynum, *, ax=None):
    if loops is None:
        return

    for loop in loops.itertuples(name=None):
        x1 = get_xbin_sqrt2((loop[2] + loop[3])/2, s, resolution) # ((loop[2] + loop[3])/2/resolution - s) * np.sqrt(2)
        x2 = get_xbin_sqrt2((loop[5] + loop[6])/2, s, resolution) # ((loop[5] + loop[6])/2/resolution - s) * np.sqrt(2)

        if (x1 >0 and x2 < fig_length):
            xmed = (x1 + min([x2, fig_length]))/2
            if ax is None:
                plt.scatter(xmed, ynum-(xmed-x1), s=60, marker="o",  facecolor='None', edgecolors='blue')
            else:
                ax.scatter(xmed, ynum-(xmed-x1), s=60, marker="o",  facecolor='None', edgecolors='blue')

def drawHeatmapTriangle_subplot2grid(matrix, resolution, *, figstart=0, figend=0,
                                     distance_max=5000000, vmin=0, vmax=50, label="",
                                     xticks=True, tads=None, loops=None, cmap=None,
                                     logratio=False, control_label="Control",
                                     heatmap_ax=None, colorbar_ax=None):
    if cmap is None:
        cmap = generate_cmap(['#FFFFFF', '#d10a3f'])

    s, e = calculate_start_end(matrix, figstart, figend, resolution)
    fig_length = calculate_fig_length(e, s)

    mat = extract_submatrix(matrix, s, e)
    rotated_mat = ndimage.rotate(mat, 45, order=0, reshape=True, prefilter=False, cval=0)

    if heatmap_ax is None:
        plt.imshow(rotated_mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect='auto')
    else:
        img = heatmap_ax.imshow(rotated_mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect='auto')

    if colorbar_ax is None:
        plt.colorbar()
    else:
        plt.colorbar(img, cax=colorbar_ax)

    ynum = rotated_mat.shape[0] / 2
    plot_tads_subplot2grid(tads, s, fig_length, resolution, ynum, ax=heatmap_ax)
    plot_loops_subplot2grid(loops, s, fig_length, resolution, ynum, ax=heatmap_ax)

    limit_yaxis_subplot2grid(distance_max, resolution, ynum, figstart, figend, ax=heatmap_ax)

    # plot xticks
    if (xticks):
        if (figend - figstart > Mega*10):
            pltxticks_mega_subplot2grid(0, fig_length, figstart, figend, ax=heatmap_ax)
        else:
            nxticks = max(int((figend - figstart)/Mega), 10)
            pltxticks_subplot2grid(0, fig_length, figstart, figend, nxticks, ax=heatmap_ax)
    else:
        xtickoff_ax(ax=heatmap_ax)
    
    ytickoff_ax(ax=heatmap_ax)

    # title
    if logratio:
        str_title = "log(" + label + "/" + control_label + ")"
    else:
        str_title = label
    if heatmap_ax is None:
        plt.title(str_title)
    else:
        heatmap_ax.set_title(str_title)

def drawHeatmapTrianglePair(matrix1, matrix2, resolution, *,
                            tads=None, loops=None, tads2=None, loops2=None,
                            figstart=0, figend=0, distance_max=5000000,
                            vmin=0, vmax=50, label="", xticks=True,
                            cmap=generate_cmap(['#FFFFFF', '#d10a3f'])):
    s, e = calculate_start_end(matrix1, figstart, figend, resolution)
    fig_length = calculate_fig_length(e, s)

    mat = extract_submatrix(matrix1, s, e)

    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)

    rotated_mat = ndimage.rotate(mat, 45, order=0, reshape=True, prefilter=False, cval=0)
    plt.imshow(rotated_mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect='auto')
    plt.colorbar()

    ynum = rotated_mat.shape[0]/2
    plot_tads_subplot2grid(tads, s, fig_length, resolution, ynum)
    plot_loops_subplot2grid(loops, s, fig_length, resolution, ynum)

    limit_yaxis_subplot2grid(distance_max, resolution, ynum, figstart, figend)

    plt.title(label)
    xtickoff_ax()
    ytickoff_ax()

    plt.subplot(2, 1, 2)

    mat2 = extract_submatrix(matrix2, s, e)
    rotated_mat = ndimage.rotate(mat2, 45, order=0, reshape=True, prefilter=False, cval=0)
    plt.imshow(rotated_mat, clim=(vmin, vmax), cmap=cmap, interpolation="nearest", aspect='auto')
    plt.colorbar()

    ynum = rotated_mat.shape[0]/2
    plot_tads_subplot2grid(tads2, s, fig_length, resolution, ynum)
    plot_loops_subplot2grid(loops2, s, fig_length, resolution, ynum)

    if (distance_max > 0):
        ystart = max(0, ynum - distance_max/resolution)
        plt.ylim(ystart, ynum)
        pltyticks_subplot2grid(ynum, ystart, 0, min(distance_max, (figend - figstart)), 5)
    else:
        ystart = 0
        plt.ylim(ystart, ynum)
        pltyticks_subplot2grid(ynum, ystart, 0, figend - figstart, 5)

    if (xticks):
        pltxticks_subplot2grid(0, (e-s)*1.41, figstart, figend, 10)
    else:
        xtickoff_ax()

    ytickoff_ax()

    plt.tight_layout()
