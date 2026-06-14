#!/bin/bash

build=hg38

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

cell=Control # siCTCF siRad21
cdir=CustardPyResults/Juicer_$build/$cell
odir=$cdir/phic
hic=$cdir/aligned/inter_30.hic

### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

chr=chr21  #chr21
start=24000000
end=32000000
resolution=100000
norm=KR
#tolerance=0.4

$sing custardpy_phic -n $norm $odir $hic $chr $start $end $resolution
