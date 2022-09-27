#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from custardpy.HiCmodule import JuicerMatrix
from custardpy.DirectionalityIndex import getDirectionalityIndexOfMultiSample
from custardpy.InsulationScore import getInsulationScoreOfMultiSample
from custardpy.loadData import *
from custardpy.PlotModule import *
from custardpy.DirectionalFreqRatio import *

def mergeState(df):
    length = df.shape[0]
    print (length)
    i_prev = 0
    for i in range(1, length):
        if (df.state[i] == df.state[i_prev]):
            df.iloc[i_prev,2] = df.iloc[i,2]
            df.iloc[i,2] = np.nan
        else:
            i_prev = i

    return df.dropna(axis = 0, how = 'any')

#import pdb
def plotDirectionalFreqRatio(plt, samples, resolution, figstart, figend, labels,
                             nrow, nrow_now, nrow_feature, args):
    s = int(figstart / resolution)
    e = int(figend   / resolution)

    smooth_median_filter = 3
    EnrichMatrices = make3dmatrixRatio(samples, smooth_median_filter)
    for i, sample in enumerate(EnrichMatrices):
        dfr = DirectionalFreqRatio(sample, resolution)
        if (args.dfr_right == True):
            array = dfr.getarrayplus()
        elif (args.dfr_left == True):
            array = dfr.getarrayminus()
        else:
            array = dfr.getarraydiff()

        if i==0: Matrix = array
        else:    Matrix = np.vstack((Matrix, array))

    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_feature, colspan=5)
    plt.imshow(Matrix[:,s:e], cmap=generate_cmap(['#1310cc', '#FFFFFF', '#d10a3f']), aspect="auto")
    plt.colorbar()
    plt.yticks(np.arange(len(labels)-1), labels[1:len(labels)])
    xtickoff(plt)

    nrow_now += nrow_feature
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=3, colspan=4)
    for i, sample in enumerate(Matrix):
        plt.plot(sample, label=labels[i+1])
#        if (args.dfr_right == True):  plt.title("Right")
 #       elif (args.dfr_left == True): plt.title("Left")
  #      else: plt.title("Right - Left")

    plt.xlim([s, e])
    pltxticks(s, e, figstart, figend, 10)

    nrow_now += 3

    from hmmlearn import hmm
    np.random.seed(1234)
    Matrix[np.isnan(Matrix)] = 0
    ncluster = 8
    model = hmm.GaussianHMM(n_components=ncluster, n_iter=10000,
                            covariance_type="diag").fit(Matrix.T)
    state = model.predict(Matrix.T)
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=3, colspan=4)
    plt.plot(state)
    plt.xlim([s, e])
    pltxticks(s, e, figstart, figend, 10)

    plt.tight_layout()

    df = pd.DataFrame(state)
    df.columns = ["state"]
    df["chr"] = args.chr
    df["start"] = np.arange(len(state)) * resolution
    df["end"] = df["start"] + args.resolution
    df = df.loc[:,["chr","start","end","state"]]
    df["score"] = "."
    df["strand"] = "."
    df["thickstart"] = df["start"]
    df["thickend"] = df["end"]

    RGBarray = ["0,0,255", "0,153,204", "0,255,255", "51,255,153",
                "0,204,51", "0,102,0", "0,0,51", "204,0,51",
                "255,0,204", "204,153,255", "204,204,102", "255,255,0"]

    df["Rgb"] = "0,0,255"
    for i in range(12):
        df.loc[df['state'] == i, 'Rgb'] = RGBarray[i]

    df = mergeState(df)
    df["end"] = df['end'].astype(int)

    df.to_csv(args.output + ".bed", sep="\t", header=False, index=False)


#@pysnooper.snoop()
def main():
    parser = argparse.ArgumentParser()
    tp = lambda x:list(map(str, x.split(':')))
    parser.add_argument("input",  help="<Input directory>:<label>", type=tp, nargs='*')
    parser.add_argument("output", help="Output prefix", type=str)
    parser.add_argument("chr",    help="chromosome", type=str)
    parser.add_argument("--type", help="normalize type", type=str, default="KR")
    parser.add_argument("--distance", help="distance for DI", type=int, default=500000)
    parser.add_argument("-r", "--resolution", help="resolution", type=int, default=25000)
    parser.add_argument("-s", "--start", help="start bp", type=int, default=0)
    parser.add_argument("-e", "--end",   help="end bp", type=int, default=1000000)
    parser.add_argument("--multi",       help="plot MultiInsulation Score", action='store_true')
    parser.add_argument("--compartment", help="plot Compartment (eigen)", action='store_true')
    parser.add_argument("--di",   help="plot Directionaly Index", action='store_true')
    parser.add_argument("--dfr",   help="plot DirectionalFreqRatio", action='store_true')
    parser.add_argument("--dfr_right",   help="(with --dfr) plot DirectionalFreqRatio (Right)", action='store_true')
    parser.add_argument("--dfr_left",   help="(with --dfr) plot DirectionalFreqRatio (Left)", action='store_true')
    parser.add_argument("--vmax", help="max value of color bar", type=int, default=50)
    parser.add_argument("--vmin", help="min value of color bar", type=int, default=0)
    parser.add_argument("-d", "--vizdistancemax", help="max distance in heatmap", type=int, default=0)
    parser.add_argument("--xsize", help="xsize for figure", type=int, default=10)
#    parser.add_argument("--ysize", help="ysize (* times of samples)", type=int, default=3)

    args = parser.parse_args()
#    print(args)

    dirs = []
    labels = []
    for input in args.input:
        dirs.append(input[0])
        if (len(input) >1):
            labels.append(input[1])
        else:
            labels.append("")

    chr = args.chr
    resolution = args.resolution
    type = args.type
    figstart = args.start
    figend = args.end
    s = int(figstart / resolution)
    e = int(figend   / resolution)
    length = figend - figstart
    binnum = e-s
    vmax = args.vmax
    vmin = args.vmin

    print ("width: " + str(length) + ", " + str(binnum) + " bins.")
    if (length <= 0):
        print ("Error: end < start.")
        exit(1)
    print (chr)
    print (resolution)

    samples = []
    for dir in dirs:
        observed = dir + "/Matrix/intrachromosomal/" + str(resolution) + "/observed."  + type + "." + chr + ".matrix.gz"
        eigen = dir + "/Eigen/" + str(resolution) + "/gd_eigen."  + type + "." + chr + ".txt"
        samples.append(JuicerMatrix("RPM", observed, resolution, eigenfile=eigen))

    nrow_heatmap = 3
    nrow_eigen = 1
    nrow_now = 0

    ### Plot
    plt.figure(figsize=(args.xsize, 10))
    nrow_feature = int(len(samples)/3)
    nrow = nrow_heatmap + nrow_eigen + nrow_feature + 6

    # Hi-C Map
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_heatmap, colspan=5)

    tadfile = dirs[0] + "/contact_domain/" + str(resolution) + "_blocks.bedpe"
    print(tadfile)
    tads = loadTADs(tadfile, chr[3:], start=figstart, end=figend)
    loopfile = dirs[0] + "/loops/merged_loops.bedpe"
    print(loopfile)
    loops = loadloops(loopfile, chr[3:], start=figstart, end=figend)

    drawHeatmapTriangle(plt, samples[0].getmatrix(), resolution,
                        figstart=figstart, figend=figend,
                        tads=tads, loops=loops,
                        vmax=vmax, vmin=vmin, distancemax=args.vizdistancemax,
                        label=labels[0], xticks=False)

    nrow_now += nrow_heatmap
    # Compartment
    plt.subplot2grid((nrow, 5), (nrow_now, 0), rowspan=nrow_eigen, colspan=4)
    plt.plot(samples[0].getEigen())
    plt.xlim([s,e])
    xtickoff(plt)

    nrow_now += nrow_eigen

    plotDirectionalFreqRatio(plt, samples, resolution, figstart, figend, labels,
                             nrow, nrow_now, nrow_feature, args)

    plt.savefig(args.output + ".pdf")

if(__name__ == '__main__'):
    main()
