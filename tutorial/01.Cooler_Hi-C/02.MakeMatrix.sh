#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.1.sif"

cell=Control
cool=CustardPyResults/Cooler_$build/$cell/cool/cooler.dedup.q30.multires.cool
odir=CustardPyResults/Cooler_$build/$cell/

resolutions="25000 50000 100000"

for resolution in $resolutions
do
    echo "generate Matrix..."
    matdir=$odir/Matrix/$resolution
    mkdir -p $matdir
    for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
    do
        echo -en "$chr.."
        $sing cooler_dump.py --resolution $resolution $cool $chr $matdir 
    done
    echo "done."
done
