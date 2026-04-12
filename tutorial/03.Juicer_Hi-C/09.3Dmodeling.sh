#!/bin/bash

build=hg38

sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.3.0.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults/Juicer_$build/$cell
chr=chr21
s=24000000
e=32000000
resolution=25000
norm=SCALE
matrix=$odir/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz
eigen=$odir/Eigen/$resolution/Compartment.$norm.$chr.All.bed
dir=$odir/pastis/${resolution}/$chr

$sing custardpy_pastis $odir $matrix $eigen $s $e $resolution
