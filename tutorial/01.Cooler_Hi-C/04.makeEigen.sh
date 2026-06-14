#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

cell=Control
cool=CustardPyResults/Cooler_$build/$cell/cool/cooler.dedup.q30.multires.cool
odir=CustardPyResults/Cooler_$build/$cell/
genome=genome.$build.fa

resolutions="25000 50000 100000"

for resolution in $resolutions
do
    echo "calculate PC1 with resolution $resolution.."
    $sing makeEigen_cool.sh -p $ncore $odir/Eigen/$resolution $cool $resolution $gt $genome
done
