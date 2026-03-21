#!/bin/bash

build=hg38

#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.2.0.0.sif"
sing="singularity exec custardpy.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults_Hi-C/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

chr=chr21
start=24000000
end=32000000
resolution=25000
norm=SCALE

$sing custardpy_phic $odir $hic $chr $start $end $resolution $norm
