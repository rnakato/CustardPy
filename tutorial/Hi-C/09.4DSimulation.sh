#!/bin/bash

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec custardpy.sif"

build=hg38
odir=CustardPyResults_Hi-C/Juicer_$build/Hap1-A
hic=$odir/aligned/inter_30.hic
chr=chr21
start=24000000
end=32000000
resolution=25000
norm=SCALE

$sing custardpy_phic $odir $hic $chr $start $end $resolution $norm
