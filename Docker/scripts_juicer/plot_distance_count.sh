#!/bin/bash

cell=$1
odir=$2
mkdir -p $odir/distance

if test -e $odir/aligned/merged_nodups.txt.gz; then
    input=$odir/aligned/merged_nodups.txt.gz
else
    input=$odir/aligned/merged_nodups.txt
fi

winsize=50000
qthre=30
prefix=distance_vs_count.MAPQ$qthre
distance_vs_count.Juicer $input $winsize $qthre > $odir/distance/$prefix.txt
Rscript $(cd $(dirname $0) && pwd)/plot_distance_count.R $odir/distance/$prefix.txt $odir/distance/$prefix.pdf $cell

distance_vs_count.Juicer.log $input $qthre > $odir/distance/$prefix.log.txt
Rscript $(cd $(dirname $0) && pwd)/plot_distance_count.log.R $odir/distance/$prefix.log.txt $odir/distance/$prefix.log.pdf $cell
