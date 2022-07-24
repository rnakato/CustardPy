#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=HindIII

gt=/work/Database/UCSC/$build/genome_table
gene=/work/Database/UCSC/$build/refFlat.txt
sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

odir=JuicerResults/Hap1-A

hic=$odir/aligned/inter_30.hic
norm=SCALE
$sing juicer_callTAD.sh $norm $odir $hic $gt
