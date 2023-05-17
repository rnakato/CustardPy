#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI
gt=genometable.$build.txt
gene=refFlat.$build.txt

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

odir=JuicerResults_$build/Hap1-A
#odir=JuicerResults_$build/WaplKO_3.3-A
hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=25000

$sing makeEigen.sh -p 24 $norm $odir $hic $resolution $gt $gene
