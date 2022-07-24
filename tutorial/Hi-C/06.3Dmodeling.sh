#!/bin/bash

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

odir=$(pwd)/JuicerResults/Hap1-A
chr=chr21
s=24000000
e=32000000
resolution=100000
norm=SCALE

$sing run-pastis.sh $odir $chr $s $e $resolution $norm
