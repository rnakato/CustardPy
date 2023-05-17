#!/bin/bash

norm=$1
hic=$2
odir=$3
build=$4

for res in 25000 50000 100000
do
    makeMatrix_intra.sh $norm $odir $hic $res $build
    #    makeMatrix_inter.sh $matrixdir "$hic" $res $build 0.0
done
