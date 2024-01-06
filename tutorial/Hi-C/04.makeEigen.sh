#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt

#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.5.0.sif"
sing="singularity exec custardpy.sif"

cell=Control # siCTCF siRad21 siNIPBL
odir=CustardPyResults_Hi-C/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

norm=SCALE
resolution=25000
ncore=24  # number of cores

$sing makeEigen.sh -p $ncore $norm $odir $hic $resolution $gt $gene
