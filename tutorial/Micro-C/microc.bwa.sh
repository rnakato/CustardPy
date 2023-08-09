#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
bwa_index=bwa-indexes/$build

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.4.0.sif"
sing="singularity exec custardpy.sif"

for cell in ESC_WT01 ESC_WT09
do
    $sing custardpy_cooler_MicroC -g $gt -f $genome -i $bwa_index -p $ncore fastq/$cell $cell

    odir=CustardPyResults_MicroC/$cell/bwa
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
