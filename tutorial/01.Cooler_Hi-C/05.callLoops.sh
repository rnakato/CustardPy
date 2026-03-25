#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.1.sif"

cell=Control
cool=CustardPyResults/Cooler_$build/$cell/cool/cooler.dedup.q30.multires.cool
odir=CustardPyResults/Cooler_$build/$cell/

resolutions="5000 10000 25000"

for resolution in $resolutions
do
    echo "call loops by cooltools.."
    $sing cooltools_dots.py $cool $odir $gt --resolution $resolution -p $ncore
done
