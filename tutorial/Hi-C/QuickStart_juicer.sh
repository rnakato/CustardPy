#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt
bwaindex=bwa-indexes/$build
ncore=64

sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.1.0.sif"
#sing="singularity exec --nv custardpy.sif"

fastq_post="_"  # "_" or "_R"
enzyme=MboI

for cell in Hap1-A SCC4KO-A WaplKO_3.3-A
do
    fqdir=fastq/$cell
    $sing custardpy_juicer -p $ncore -a $gene -b $build -g $gt -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell
done
