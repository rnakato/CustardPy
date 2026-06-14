#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
index_bwa=bwa-indexes/$build

sing="apptainer exec --nv --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.1.sif"

cell=C36_rep1
$sing custardpy_cooler -g $gt -f $genome -b $build -e None -i $index_bwa -p $ncore fastq/$cell $cell

# Downstream analysis of .cool files (eigenvector, insulation score, loop calling, etc.)
# Because of the small number of mapped reads, we used loop resolutions of 25 kb here.
odir=CustardPyResults/Cooler_$build/$cell
cool=$odir/cool/contact.bwa.q30.mcool
$sing custardpy_process_cool -t $ncore -s "25000" -g $gt -f $genome $cool $odir $cell

## loop calling with fithic (take long time, so run separately)
resolution=25000
odir_fithic=$odir/loops/fithic/$resolutions_fithic
$sing run_fithic.sh -g $gt $cool $odir_fithic $cell $resolution \
    | tee $odir/log/fithic.$resolution.log

mkdir -p $odir/loops/mustache
$sing mustache -f $cool -r 5kb -norm weight -pt 0.05 -p $ncore \
      -o $odir/loops/mustache/loop.cool.5kb.tsv


odir=CustardPyResults/Cooler_mm39/$cell
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
#$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
