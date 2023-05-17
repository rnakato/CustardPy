#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI
#gt=/work/Database/UCSC/$build/genome_table
gt=genometable.$build.txt

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

odir=JuicerResults_$build/Hap1-A

hic=$odir/aligned/inter_30.hic
norm=SCALE
$sing juicer_callTAD.sh $norm $odir $hic $gt
