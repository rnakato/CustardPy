#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt
bwaindex=bwa-indexes/$build
ncore=64

sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.2.2.1.sif"
#sing="singularity exec --nv custardpy.sif"

fastq_post="_"  # "_" or "_R"
enzyme=MboI

for cell in Control siCTCF siRad21 siNIPBL
do
    fqdir=fastq/$cell
    $sing custardpy_juicer -p $ncore -a $gene -b $build -g $gt -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell
    odir=CustardPyResults_Hi-C/Juicer_$build/$cell
    hic=$odir/aligned/inter_30.hic
    norm=SCALE
    $sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
done

odir=CustardPyResults_Hi-C/Juicer_$build/
norm=SCALE
$sing Juicerstats.sh $odir $norm

$sing execute_R plot_distance_count_all.R $odir $odir/plot_distance_count_all.pdf
$sing execute_R plot_distance_count_all.log.R $odir $odir/plot_distance_count_all.log.pdf
