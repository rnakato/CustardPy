#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
bwa_index=bwa-indexes/$build

sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.6.0.sif"
#sing="singularity exec custardpy.sif"

cell=C36_rep1
$sing custardpy_cooler_MicroC -g $gt -f $genome -i $bwa_index -p $ncore fastq/$cell $cell

odir=CustardPyResults_MicroC/Cooler_bwa/$cell
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
