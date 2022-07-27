#!/bin/bash

build=mm39
ncore=64
gt=genometable.mm39.txt
gene=refFlat.mm39.txt

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"
sing_juicer="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

fq1=fastq/SRR8954797_1.fastq.gz
fq2=fastq/SRR8954797_2.fastq.gz
prefix=ESC_WT01

bwa_index=bwa-indexes/mm39
$sing custardpy_mappingMicroC -t bwa -i $bwa_index -g $gt -p $ncore $fq1 $fq2 $prefix

odir=Cooler_MicroC_bwa/ESC_WT01/
hic=$odir/hic/contact_map.q30.hic
norm=SCALE

resolution=100000

# Contact matrix
$sing_juicer makeMatrix_intra.sh $norm $odir $hic $resolution $gt
# InsulationScore
$sing_juicer makeInslationScore.sh $norm $odir $resolution $gt
# TAD
$sing_juicer juicer_callTAD.sh $norm $odir $hic $gt
# Eigen
$sing_juicer makeEigen.sh $norm $odir $hic $resolution $gt $gene
# Loop
singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif call_HiCCUPS.sh $norm $odir $hic
