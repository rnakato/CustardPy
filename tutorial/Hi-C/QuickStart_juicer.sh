#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt
bwaindex=bwa-indexes/$build
ncore=64

sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

cell=Hap1-A
#cell=Control
fastq_post="_"  # "_" or "_R"
enzyme=MboI

fqdir=fastq/$cell
$sing custardpy_juicer -p $ncore -a $gene -b $build -g $gt -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell
