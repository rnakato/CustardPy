#!/bin/bash

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

odir=JuicerResults_hg38/Hap1-A
chr=chr21
s=24000000
e=32000000
resolution=25000
norm=SCALE

$sing custardpy_pastis $odir $chr $s $e $resolution $norm
