#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=HindIII
gt=genometable.$build.txt
gene=refFlat.$build.txt

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

odir=JuicerResults/Hap1-A
hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000

$sing makeEigen.sh $norm $odir $hic $resolution $gt $gene
