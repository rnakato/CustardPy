#!/bin/bash

build=hg38
gt=genometable.$build.txt

#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.5.0.sif"
sing="singularity exec custardpy.sif"

cell=Control # siCTCF siRad21 
odir=CustardPyResults_Hi-C/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

norm=SCALE  # KR VC SQRT VC_SQRT NONE
resolution=25000

# Contact matrix
echo "generate Matrix..."
$sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt

# Insulation Score
echo "calculate Insulation score.."
$sing makeInslationScore.sh $norm $odir $resolution $gt
