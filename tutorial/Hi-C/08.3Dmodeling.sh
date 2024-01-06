#!/bin/bash

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.6.0.sif"
#sing="singularity exec custardpy.sif"

build=hg38
odir=CustardPyResults_Hi-C/Juicer_$build/Control
chr=chr21
s=24000000
e=32000000
resolution=25000
norm=SCALE

$sing custardpy_pastis $odir $chr $s $e $resolution $norm
