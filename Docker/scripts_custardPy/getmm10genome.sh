#!/bin/bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
mkdir -p tempdir
splitmultifasta mm10.fa --dir tempdir

for chr in $(seq 1 19) X Y M
do
    fa=tempdir/chr$chr.fa
    s="$s $fa"
done
cat $s > genome.mm10.fa
rm -rf tempdir mm10.fa
