#!/bin/bash

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

odir=JuicerResults_$build/Hap1-A
hic=$odir/aligned/inter_30.hic
chr=chr21
start=24000000
end=32000000
resolution=100000
norm=SCALE

$sing custardpy_phic $odir $hic $chr $start $end $resolution $norm
