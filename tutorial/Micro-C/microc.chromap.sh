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

chromap_index=chromap-indexes/$build
genome=genome.$build.fa
#$sing custardpy_mappingMicroC -t chromap -i $chromap_index -g $gt -f $genome -p $ncore $fq1 $fq2 $prefix

odir=Cooler_MicroC_chromap/$prefix
hic=$odir/hic/contact_map.q30.hic
norm=SCALE

$sing_juicer custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

#$sing_juicer call_HiCCUPS.sh $norm $idir $hic