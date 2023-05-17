#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
bwa_index=bwa-indexes/$build

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.0.0.sif"
sing="singularity exec custardpy.sif"

fq1=fastq/SRR8954797_1.fastq.gz
fq2=fastq/SRR8954797_2.fastq.gz
prefix=ESC_WT01
#$sing custardpy_mappingMicroC -t bwa -i $bwa_index -g $gt -f $genome -p $ncore $fq1 $fq2 $prefix

odir=Cooler_MicroC_bwa/$prefix
hic=$odir/hic/contact_map.q30.hic
norm=SCALE

$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
