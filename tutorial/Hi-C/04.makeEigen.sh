#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI

gt=/work/Database/UCSC/$build/genome_table
gene=/work/Database/UCSC/$build/refFlat.txt
sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

odir=$(pwd)/JuicerResults/Hap1-A

hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000

# Pearson
#$sing makeEigen.sh Pearson $norm $odir $hic $resolution $gt $gene
# Eigen
#$sing makeEigen.sh Eigen $norm $odir $hic $resolution $gt $gene
$sing makeEigen.sh $norm $odir $hic $resolution $gt $gene
