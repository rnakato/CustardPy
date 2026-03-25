#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.1.1.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

norm=SCALE
$sing juicer_callTAD.sh $norm $odir $hic $gt
