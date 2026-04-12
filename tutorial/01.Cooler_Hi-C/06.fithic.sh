#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.3.0.sif"

cell=Control
odir=CustardPyResults/Cooler_$build/$cell/

mkdir -p $odir/log
resolutions="10000 25000"

for resolution in $resolutions
do
    echo "Run Fit-Hi-C.."
    pairfile=$odir/pairs/dedup.bwa.q30.pairs.gz
    odir_fithic=$odir/loops/fithic/$resolution
    $sing run_fithic.sh -g $gt $pairfile $odir_fithic $cell $resolution | tee $odir/log/fithic.$resolution.log
done
