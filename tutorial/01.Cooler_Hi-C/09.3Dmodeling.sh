#!/bin/bash

build=hg38

sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.3.0.sif"

cell=Control
odir=CustardPyResults/Cooler_$build/$cell
chr=chr21
s=24000000
e=32000000
resolution=25000
matrix=$odir/Matrix/$resolution/balanced.$chr.matrix.gz
eigen=$odir/Eigen/$resolution/PC1.cis.genome.All.bed12
dir=$odir/pastis/${resolution}/$chr

mkdir -p $dir

awk -F'\t' '$1=="chr1"' $eigen | cut -f1-4 > $dir/PC1.cis.genome.All.$chr.bed

$sing custardpy_pastis $odir $matrix $dir/PC1.cis.genome.All.$chr.bed $s $e $resolution
