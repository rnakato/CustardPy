#!/bin/bash

build=mm39
ncore=64
gt=genometable.mm39.txt

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"
sing_juicer="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

fq1=fastq/SRR8954797_1.fastq.gz
fq2=fastq/SRR8954797_2.fastq.gz
prefix=ESC_WT01

bwa_index=bwa-indexes/mm39
#$sing custardpy_mappingMicroC -t bwa -i $bwa_index -g $gt -p $ncore $fq1 $fq2 $prefix

idir=Cooler_MicroC_bwa/ESC_WT01/
hic=$idir/hic/contact_map.q30.hic
norm=SCALE
#$sing_juicer call_HiCCUPS.sh $norm $idir $hic

#idir=Results_chromap/$prefix
chromap_index=chromap-indexes/$build
genome=genome.$build.fa
#$sing custardpy_mappingMicroC -t chromap -i $chromap_index -g $gt -f $genome -p $ncore $fq1 $fq2 $prefix
hic=$idir/hic/contact_map.q30.hic
norm=SCALE
$sing_juicer call_HiCCUPS.sh $norm $idir $hic
