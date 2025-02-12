#!/bin/bash

build=hg38
gt=genometable.hg38.txt
index_bwa=bwa-indexes/hg38
gene=refFlat.$build.txt
genome=genome.$build.fa
ncore=64
enzyme=MboI

sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.2.2.1.sif"
#sing="singularity exec --nv custardpy.sif"

for cell in Control siCTCF siRad21 siNIPBL
do
    $sing custardpy_cooler_HiC -g $gt -b $build -f $genome -i $index_bwa -p $ncore fastq/$cell $cell
    odir=CustardPyResults_Hi-C/Cooler_$build/$cell
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
