#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
index_bwa=bwa-indexes/$build

sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.1.0.sif"
#sing="singularity exec custardpy.sif"

cell=C36_rep1
#$sing custardpy_cooler -g $gt -f $genome -b $build -e None -i $index_bwa -p $ncore fastq/$cell $cell
#exit
odir=CustardPyResults/Cooler_mm39/$cell
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
