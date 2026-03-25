#!/bin/bash

build=hg38

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.1.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

chr=chr21
start=24000000
end=32000000
resolution=25000
norm=SCALE

$sing custardpy_phic $odir $hic $chr $start $end $resolution $norm
