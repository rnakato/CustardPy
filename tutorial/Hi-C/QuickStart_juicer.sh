#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt
bwaindex=bwa-indexes/$build
ncore=64

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.4.0.sif"
#sing="singularity exec --nv custardpy_juicer.sif"
sing="singularity exec --nv --bind /work,/work2 ../../Docker/custardpy.1.0.0.sif"

fastq_post="_"  # "_" or "_R"
enzyme=MboI

for cell in Hap1-A #SCC4KO-A #Hap1-A WaplKO_3.3-A
do
    fqdir=fastq/$cell
    $sing custardpy_juicer -p $ncore -a $gene -b $build -g $gt -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell
done
