#!/bin/bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mkdir -p tempdir
splitmultifasta hg38.fa --dir tempdir

for chr in $(seq 1 22) X Y M
do
    fa=tempdir/chr$chr.fa
    s="$s $fa"
done
cat $s > genome.hg38.fa
rm -rf tempdir hg38.fa
makegenometable.pl genome.hg38.fa > genometable.hg38.txt

wget https://nakatolab.iqb.u-tokyo.ac.jp/Datafolder_for_sharing/DockerDatabase/refFlat/refFlat.hg38.txt.gz
unpigz refFlat.hg38.txt.gz
