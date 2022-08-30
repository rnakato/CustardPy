#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI

gt=/work/Database/UCSC/$build/genome_table
sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

odir=JuicerResults_$build/Hap1-A
hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000
$sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt
$sing makeInslationScore.sh $norm $odir $resolution $gt
