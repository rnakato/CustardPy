#!/bin/bash

build=hg38
gt=genometable.hg38.txt
index_bwa=bwa-indexes/hg38
gene=refFlat.$build.txt
ncore=64

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.2.0.sif"
sing="singularity exec custardpy.sif"

fq1=fastq/Hap1-A/SRR5266584_1.fastq.gz
fq2=fastq/Hap1-A/SRR5266584_2.fastq.gz
prefix=Hap1-A
enzyme=MboI

# generate .hic and .cool files from fastq
$sing custardpy_mappingHiC -g $gt -i $index_bwa \
      -b $build -e $enzyme -p $ncore \
      $fq1 $fq2 $prefix

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.2.0.sif"
sing="singularity exec --nv custardpy_juicer.sif"

odir=Cooler_$build/$prefix
hic=$odir/hic/contact_map.q30.hic
norm=SCALE
#resolution=100000

$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

# Contact matrix
#$sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt
# InsulationScore
#$sing makeInslationScore.sh $norm $odir $resolution $gt
# TAD
#$sing juicer_callTAD.sh $norm $odir $hic $gt
# Eigen
#$sing makeEigen.sh $norm $odir $hic $resolution $gt $gene
# Loop
#$sing call_HiCCUPS.sh $norm $odir $hic
