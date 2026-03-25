#!/bin/bash

build=hg38

sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.1.1.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults/Juicer_$build/$cell
chr=chr21
s=24000000
e=32000000
resolution=25000
norm=SCALE

$sing custardpy_pastis $odir $chr $s $e $resolution $norm
