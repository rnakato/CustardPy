#!/bin/bash

build=hg38
gt=genometable.$build.txt
gene=refFlat.$build.txt
bwaindex=bwa-indexes/$build
ncore=64

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

cell=Hap1-A
#cell=Control
fastq_post="_"  # "_" or "_R"
enzyme=HindIII

fqdir=fastq/$cell

$sing custardpy_juicer -p $ncore -a $gene -b $build -g $gt -i $bwaindex -e $enzyme -z $fastq_post  $fqdir $cell
exit



odir=JuicerResults/$cell
rm -rf $odir
$sing juicer_map.sh -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

# Optional
$sing juicer_pigz.sh $odir
$sing plot_distance_count.sh $cell $odir

hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000
# Contact matrix
$sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt
# InsulationScore
$sing makeInslationScore.sh $norm $odir $resolution $gt
# TAD
$sing juicer_callTAD.sh $norm $odir $hic $gt
# Eigen
$sing makeEigen.sh $norm $odir $hic $resolution $gt $gene
# Loop
$sing call_HiCCUPS.sh $norm $odir $hic
