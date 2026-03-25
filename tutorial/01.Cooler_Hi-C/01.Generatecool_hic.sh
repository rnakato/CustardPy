#!/bin/bash

build=hg38
gt=genometable.$build.txt
index_bwa=bwa-indexes/$build
genome=genome.$build.fa
ncore=64
enzyme=MboI
fastq_post="_"  # "_" or "_R"

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.1.sif"

cell=Control #siCTCF siRad21 siNIPBL

$sing custardpy_cooler -g $gt -f $genome -b $build -e $enzyme -z $fastq_post -i $index_bwa -p $ncore fastq/$cell $cell