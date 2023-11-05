#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
#gt=/work/Database/UCSC/$build/genome_table
sing="singularity exec custardpy.sif"
gt=genometable.$build.txt

odir=CustardPyResults_Hi-C/Juicer_$build/Hap1-A

hic=$odir/aligned/inter_30.hic
norm=SCALE
$sing juicer_callTAD.sh $norm $odir $hic $gt
