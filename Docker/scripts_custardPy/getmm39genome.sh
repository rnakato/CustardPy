#!/bin/bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz
mkdir -p tempdir
splitmultifasta mm39.fa --dir tempdir

for chr in $(seq 1 19) X Y M
do
    fa=tempdir/chr$chr.fa
    s="$s $fa"
done
cat $s > genome.mm39.fa
rm -rf tempdir mm39.fa
makegenometable.pl genome.mm39.fa > genometable.mm39.txt

wget https://nakatolab.iqb.u-tokyo.ac.jp/Datafolder_for_sharing/DockerDatabase/refFlat/refFlat.mm39.txt.gz
unpigz refFlat.mm39.txt.gz
