#!/bin/bash

build=hg38
gt=genometable.hg38.txt
index_bwa=bwa-indexes/hg38
ncore=64

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

#fq1=fastq/Hap1-A/SRR5266584_1.fastq.gz
#fq2=fastq/Hap1-A/SRR5266584_2.fastq.gz
#prefix=Hap1-A

fq1=fastq/Control/SRR5952305_1.fastq.gz
fq2=fastq/Control/SRR5952305_2.fastq.gz
prefix=Control
enzyme=MboI

# generate .hic and .cool files from fastq
$sing custardpy_mappingHiC -g $gt -i $index_bwa \
      -b $build -e $enzyme -p $ncore \
      $fq1 $fq2 $prefix
