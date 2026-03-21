#!/bin/bash

build=mm39
ncore=64
gt=genometable.$build.txt
gene=refFlat.$build.txt
genome=genome.$build.fa
index_bwa=bwa-indexes/$build

sing="singularity exec --nv --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.0.sif"
#sing="singularity exec custardpy.sif"

cell=C36_rep1
#$sing custardpy_cooler -g $gt -f $genome -b $build -e None -i $index_bwa -p $ncore fastq/$cell $cell
#exit

# downstream analysis of .cool files (eigenvector, insulation score, loop calling, etc.)
odir=CustardPyResults/Cooler_$build/$cell
cool=$odir/cool/cooler.dedup.q30.multires.cool
$sing custardpy_process_cool -t $ncore -g $gt -f $genome $cool $odir $cell
exit

## loop calling with fithic (take long time, so run separately)
resolutions_fithic=25000
pair=$odir/pairs/dedup.bwa.q30.pairs.gz
odir_fithic=$odir/loops/fithic/$resolutions_fithic
#    $sing run_fithic.sh -g $gt $pairfile $odir_fithic $label $resolution | tee $odir/log/fithic.$resolution.log
exit

odir=CustardPyResults/Cooler_mm39/$cell
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
