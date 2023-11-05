#!/bin/bash

build=mm39
ncore=64
gt=genometable.mm39.txt
gene=refFlat.mm39.txt
genome=genome.$build.fa
chromap_index=chromap-indexes/$build

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec custardpy.sif"

for cell in ESC_WT01 ESC_WT09
do
    $sing custardpy_cooler_MicroC -g $gt -f $genome -i $chromap_index -p $ncore -t chromap fastq/$cell $cell

    odir=CustardPyResults_MicroC/$cell/chromap
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
