#!/bin/bash

build=mm39
ncore=64
gt=genometable.mm39.txt
gene=refFlat.mm39.txt
genome=genome.$build.fa
chromap_index=chromap-indexes/$build

sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.6.0.sif"
#sing="singularity exec custardpy.sif"

cell=C36_rep1
$sing custardpy_cooler_MicroC -g $gt -f $genome -i $chromap_index -p $ncore -t chromap fastq/$cell $cell

odir=CustardPyResults_MicroC/$cell/chromap
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
