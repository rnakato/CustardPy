#!/bin/bash

build=hg38
gt=genometable.$build.txt
index_bwa=bwa-indexes/$build
genome=genome.$build.fa
ncore=64
enzyme=MboI

sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.0.sif"
#sing="singularity exec custardpy.sif"

for cell in Control siCTCF siRad21 siNIPBL
do
    $sing custardpy_cooler -g $gt -f $genome -b $build -e $enzyme -i $index_bwa -p $ncore fastq/$cell $cell

    odir=CustardPyResults_Hi-C/Cooler_$build/$cell

    cool=$odir/cool/cooler.dedup.q30.multires.cool
    pair=$odir/pairs/dedup.bwa.q30.pairs.gz
    $sing custardpy_process_cool -t $ncore -g $gt -p $pair -f $genome $cool $odir $cell

    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
#    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done
