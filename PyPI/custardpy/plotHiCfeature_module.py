import matplotlib.pyplot as plt

from custardpy.HiCmodule import JuicerMatrix
from custardpy.DirectionalityIndex import getDirectionalityIndexOfMultiSample
from custardpy.InsulationScore import getInsulationScoreOfMultiSample
from custardpy.generateCmap import *
from custardpy.loadData import *
from custardpy.PlotModule import *
from custardpy.DirectionalRelativeFreq import *

def plot_legend_in_subplot2grid(ax_origin, ax_legend):
    ax_legend.axis('off')
    _lines, _labels = ax_origin.get_legend_handles_labels()
    ax_legend.legend(_lines, _labels, loc='center')


def set_figsize_x(xsize, figstart, figend):
    if xsize == 0:
        figsize_x = max(int((figend-figstart)/2000000), 10)
    else:
        figsize_x = xsize
    return figsize_x


def plot_PC1(nrow, nrow_now, nrow_eigen, sample, label, s, e, colspan_plot, colspan_full):
    plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_eigen, colspan=colspan_plot)
    plt.plot(sample.getEigen(), color="black")
    plt.xlim([s,e])
    xtickoff_ax()

    plt.title("Compartment PC1 (" + label + ")")


def plot_xy_axis_and_title_of_feature_heatmap(ax, labels, title):
    xtickoff_ax(ax=ax)
    ax.set_yticks(np.arange(len(labels)), labels)
    ax.set_title(title)


def get_drf_array(sample, resolution, drf_right, drf_left):
    drf = DirectionalRelativeFreq(sample, resolution)
    if drf_right:
        return drf.getarrayplus()
    elif drf_left:
        return drf.getarrayminus()
    else:
        return drf.getarraydiff()


def get_drf_matrix(samples, resolution, drf_right, drf_left, *, smooth_median_filter=3):
    if not isinstance(samples, list):
        samples = [samples]
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)
    arrays = [get_drf_array(sample, resolution, drf_right, drf_left) for sample in EnrichMatrices]
    Matrix = np.vstack(arrays)

    return Matrix


def plot_directional_relative_frequency(samples, labels,  nrow, nrow_now, nrow_feature,
                                        s, e, figstart, figend, resolution,
                                        drf_right, drf_left,
                                        colspan_plot, colspan_colorbar, colspan_legend, colspan_full):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    DRFMatrix = get_drf_matrix(samples, resolution, drf_right, drf_left)

    # DRF heatmap
    heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
    colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)

    img = heatmap_ax.imshow(DRFMatrix[:,s:e], cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']), aspect="auto")
    plot_xy_axis_and_title_of_feature_heatmap(heatmap_ax, labels[1:], "Directional Relative Frequency")

    plt.colorbar(img, cax=colorbar_ax)

    nrow_now += nrow_feature

    # DRF barplot
    heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
    legend_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot), rowspan=nrow_feature, colspan=colspan_legend)
    for i, sample in enumerate(DRFMatrix):
        heatmap_ax.plot(sample, label=labels[i+1])

    heatmap_ax.set_xlim([s, e])
    pltxticks_subplot2grid(s, e, figstart, figend, 10, ax=heatmap_ax)
    plot_legend_in_subplot2grid(heatmap_ax, legend_ax)


def plot_triangle_ratio_multi(samples, labels, nrow, nrow_now, nrow_heatmap, nrow_feature,
                              s, e, figstart, figend, distance_max, resolution,
                              vmin_ratio, vmax_ratio,
                              colspan_plot, colspan_colorbar, colspan_legend, colspan_full,
                              *, smooth_median_filter=3):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)

    for i, sample in enumerate(EnrichMatrices):
        # Hi-C logratio Map
        heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_heatmap, colspan=colspan_plot)
        colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_heatmap, colspan=colspan_colorbar)

        drawHeatmapTriangle_subplot2grid(sample, resolution,
                                         figstart=figstart, figend=figend,
                                         vmax=vmax_ratio, vmin=vmin_ratio,
                                         cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                                         distance_max=distance_max,
                                         label=labels[i+1], xticks=True,
                                         logratio=True, control_label=labels[0],
                                         heatmap_ax=heatmap_ax, colorbar_ax=colorbar_ax)
        nrow_now += nrow_heatmap

        # DRF barplot (left & right)
        heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        legend_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_legend)

        drf = DirectionalRelativeFreq(sample, resolution)

        heatmap_ax.plot(drf.getarrayplus(), label="Right")
        heatmap_ax.plot(drf.getarrayminus(), label="Left")
        heatmap_ax.set_xlim([s,e])
        xtickoff_ax(ax=heatmap_ax)
        plot_legend_in_subplot2grid(heatmap_ax, legend_ax)

        nrow_now += nrow_feature

        # DRF histogram
        heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)

        diff = drf.getarraydiff()
        heatmap_ax.bar(range(len(diff)), diff)
        heatmap_ax.set_xlim([s,e])
        xtickoff_ax(ax=heatmap_ax)
        heatmap_ax.set_title("Directional Relative Frequency (Right - Left)")

        nrow_now += nrow_feature

#    plt.tight_layout()


def plot_directionality_index(samples, labels, nrow, nrow_now, nrow_feature,
                              s, e, figstart, figend, distance,
                              colspan_plot, colspan_colorbar, colspan_legend, colspan_full,
                                 *, vmin=-1000, vmax=1000, heatmap=False, lineplot=True):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    vDI = getDirectionalityIndexOfMultiSample(samples, labels, distance=distance)
    if len(samples) == 1:
        vDI = vDI.reshape((1, -1))

    if(heatmap):
        heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)

        img = heatmap_ax.imshow(vDI[:,s:e],
                                cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                                aspect="auto")
        plot_xy_axis_and_title_of_feature_heatmap(heatmap_ax, labels, "Directionality Index")

        plt.colorbar(img, cax=colorbar_ax)

        nrow_now += nrow_feature

    if(lineplot):
        heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        legend_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot), rowspan=nrow_feature, colspan=colspan_legend)

        for i, sample in enumerate(samples):
            heatmap_ax.plot(vDI[i], label=labels[i])
        heatmap_ax.set_xlim([s, e])
        heatmap_ax.set_ylim([vmin, vmax])
        pltxticks_subplot2grid(s, e, figstart, figend, 10, ax=heatmap_ax)

        plot_legend_in_subplot2grid(heatmap_ax, legend_ax)


def plot_virtual_4c(samples, labels, nrow, nrow_now, nrow_feature,
                    s, e, figstart, figend, anchor, resolution, vmin, vmax,
                    colspan_plot, colspan_full):
    a = int(anchor/resolution)

    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    for i, sample in enumerate(samples):
        heatmap_ax = plt.subplot2grid((nrow, colspan_full), (i + nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)

        # Dataframeが値を持っているので、s-eではなくfigstart-figendになってしまう
        # それを避けるためvaluesをつける
        heatmap_ax.plot(sample.getmatrix().iloc[a].values, color="r")
        heatmap_ax.set_xlim([s, e])
        heatmap_ax.set_ylim(vmin, vmax)
        heatmap_ax.set_title(labels[i])

        if s < a and a < e:
            arrow_x = a
            arrow_y = np.sin(arrow_x)
            heatmap_ax.axvline(x=arrow_x, color="black", linestyle="--", linewidth=0.5)
            plt.annotate("", xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y - 0.5),
                         arrowprops=dict(facecolor="black", edgecolor="black"))

        pltxticks_subplot2grid(s, e, figstart, figend, 10, ax=heatmap_ax)


def plot_compartment_heatmap(samples, labels, nrow, nrow_now, nrow_feature,
                             s, e, figstart, figend,
                             colspan_plot, colspan_colorbar, colspan_legend, colspan_full):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    for i, sample in enumerate(samples):
        if i==0: Matrix = sample.getEigen()
        else:    Matrix = np.vstack((Matrix, sample.getEigen()))

    if len(samples) == 1:
        Matrix = Matrix.reshape((1, -1))

    # PC1 heatmap
    heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
    colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)
    img = heatmap_ax.imshow(Matrix[:,s:e], clim=(-0.05, 0.05),
                            cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']),
                            aspect="auto")
    plot_xy_axis_and_title_of_feature_heatmap(heatmap_ax, labels, "Compartment PC1")

    plt.colorbar(img, cax=colorbar_ax)

    nrow_now += nrow_feature

    # PC1 barplot
    heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
    legend_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot), rowspan=nrow_feature, colspan=colspan_legend)

    for i, sample in enumerate(samples):
        heatmap_ax.plot(sample.getEigen(), label=labels[i])
    heatmap_ax.set_xlim([s, e])
    heatmap_ax.set_ylim([-0.05, 0.05])
    pltxticks_subplot2grid(s, e, figstart, figend, 10, ax=heatmap_ax)

    plot_legend_in_subplot2grid(heatmap_ax, legend_ax)


def plot_differential_multi_insulation_score(samples, labels, nrow, nrow_now, nrow_feature,
                                                 figstart, figend, s, e,
                                                 colspan_plot, colspan_colorbar, colspan_full):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    for i, sample in enumerate(samples):
        if i == 0:
            MIref = sample.getMultiInsulationScore()
            continue

        MI = sample.getMultiInsulationScore()

        heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (i-1 + nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        colorbar_ax = plt.subplot2grid((nrow, colspan_full), (i-1 + nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)

        img = heatmap_ax.imshow(MI.T.iloc[:,s:e] - MIref.T.iloc[:,s:e], clim=(-0.4, 0.4),
                                cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                                aspect="auto")
        heatmap_ax.set_title(labels[i])

        if i < len(samples) -1:
            xtickoff_ax(ax=heatmap_ax)
        else:
            pltxticks_subplot2grid(0, e-s, figstart, figend, 10, ax=heatmap_ax)

        plt.colorbar(img, cax=colorbar_ax)


def plot_multi_insulation_score(samples, labels, nrow, nrow_now, nrow_feature,
                                     figstart, figend, s, e,
                                     colspan_plot, colspan_colorbar, colspan_full):
    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    for i, sample in enumerate(samples):
        MI = sample.getMultiInsulationScore()

        heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (i + nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        colorbar_ax = plt.subplot2grid((nrow, colspan_full), (i + nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)

        img = heatmap_ax.imshow(MI.T.iloc[:,s:e], clim=(0.4, 1.0),
                         cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                         aspect="auto")
        heatmap_ax.set_title(labels[i])

        if i < len(samples) -1:
            xtickoff_ax(ax=heatmap_ax)
        else:
            pltxticks_subplot2grid(0, e-s, figstart, figend, 10, ax=heatmap_ax)

        plt.colorbar(img, cax=colorbar_ax)


def plot_single_insulation_score(samples, labels, nrow, nrow_now, nrow_feature,
                                 figstart, figend, s, e,
                                 colspan_plot, colspan_colorbar, colspan_legend, colspan_full,
                                 *, heatmap=False, lineplot=True):
    Matrix = getInsulationScoreOfMultiSample(samples, labels)

    if not isinstance(samples, list):
        samples = [samples]
    if not isinstance(labels, list):
        labels = [labels]

    if(heatmap):
        heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1), rowspan=nrow_feature, colspan=colspan_colorbar)

        img = heatmap_ax.imshow(Matrix.T.iloc[:,s:e], clim=(0.4, 1.0),
                            cmap=generate_cmap(['#d10a3f', '#FFFFFF', '#1310cc']),
                            aspect="auto")
        plot_xy_axis_and_title_of_feature_heatmap(heatmap_ax, labels, "Insulation score")

        plt.colorbar(img, cax=colorbar_ax)

        nrow_now += nrow_feature

    if(lineplot):
        heatmap_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0), rowspan=nrow_feature, colspan=colspan_plot)
        legend_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot), rowspan=nrow_feature, colspan=colspan_legend)

        for i in range(len(samples)):
            heatmap_ax.plot(Matrix.iloc[s:e,i], label=labels[i])
        heatmap_ax.set_xlim([figstart, figend])
        pltxticks_subplot2grid(figstart, figend, figstart, figend, 10, ax=heatmap_ax)

        plot_legend_in_subplot2grid(heatmap_ax, legend_ax)


def plot_HiC_Map(nrow, nrow_now, nrow_heatmap, sample, label, chr,
                 type, resolution, vmax, vmin, figstart, figend, distance_max,
                 colspan_plot, colspan_colorbar, colspan_full,
                 *, tadfile="", loopfile=""):
    # load TADs
    if(tadfile != ""):
        print(tadfile)
        tads = loadTADs(tadfile, chr, start=figstart, end=figend)
    else:
        tads=None

    # load loops
    if(loopfile != ""):
        print(loopfile)
        loops = loadloops(loopfile, chr, start=figstart, end=figend)
    else:
        loops=None

    # plot Hi-C
    heatmap_ax  = plt.subplot2grid((nrow, colspan_full), (nrow_now, 0),
                                    rowspan=nrow_heatmap, colspan=colspan_plot)
    colorbar_ax = plt.subplot2grid((nrow, colspan_full), (nrow_now, colspan_plot +1),
                                    rowspan=nrow_heatmap, colspan=colspan_colorbar)
    drawHeatmapTriangle_subplot2grid(sample.getmatrix(), resolution, figstart=figstart, figend=figend,
                                     tads=tads, loops=loops, vmax=vmax, vmin=vmin, distance_max=distance_max,
                                     label="Contact map (" + label + ")",
                                     xticks=True, heatmap_ax=heatmap_ax, colorbar_ax=colorbar_ax)

def plot_SamplePair_triu(sample1, sample2, label1, label2,
                        resolution, figstart, figend, vmax, vmin,
                        *, log=False):
    if (log):
        m1 = sample1.getlog()
        m2 = sample2.getlog()
    else:
        m1 = sample1.getmatrix()
        m2 = sample2.getmatrix()

    matrix = np.triu(m1) + np.tril(m2, k=-1)

    if (label1 != "" and label2 != ""):
        label = label1 + "(upper) - " + label2 + "(lower)"
    else:
        label = ""

    drawHeatmapSquare(matrix, resolution,
                      figstart=figstart, figend=figend,
                      vmax=vmax, vmin=vmin,
                      label=label, xticks=True)