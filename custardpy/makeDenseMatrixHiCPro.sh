#!/bin/bash

scriptname=${0##*/}
if test $# -ne 6; then
    echo "$scriptname <human|mouse> <outputdir> <expname> <rawmatrix> <icematrix> <absbed>"
    exit 0
fi

odir=$2
expname=$3
rawmatrix=$4
icematrix=$5
absbed=$6

dmdir=$odir/denseMatrix
abdir=$odir/absBedChr
mkdir -p $dmdir $abdir

HiCProdir=$(cd $(dirname $0) && pwd)/../HiC-Pro_2.10.0/bin/utils

if test ! -e $dmdir/chr1_$expname.raw_dense.matrix; then
    python $HiCProdir/sparseToDense.py $rawmatrix -b $absbed -c -o $expname.raw_dense.matrix
    mv *$expname.raw_dense.matrix $dmdir
fi
if test ! -e $dmdir/chr1_$expname.iced_dense.matrix; then
    python $HiCProdir/sparseToDense.py $icematrix -b $absbed -c -o $expname.iced_dense.matrix
    mv *$expname.iced_dense.matrix $dmdir
fi

if test $1 = "human"; then
    chrarray="$(seq 1 22) X"
elif test $1 = "mouse"; then
    chrarray="$(seq 1 19) X"
else
    echo "Error: specify <human|mouse>"
fi

for i in $chrarray; do
    if test ! -e $abdir/chr${i}_abs.bed; then
	cat $absbed | awk -v chr=chr$i '$1 == chr {print}' > $abdir/chr${i}_abs.bed
    fi
done
