#!/bin/bash

cell=$1
fqdir=`pwd`/fastq
odir=`pwd`/JuicerResults/$cell  # 絶対パス
#echo "Results will be output in '$odir'"

mkdir -p $odir
if test ! -e $odir/fastq; then
    ln -s $fqdir/$cell/ $odir/fastq;
fi

echo $odir
