#!/bin/bash

build=hg38
gt=genometable.$build.txt
bwaindex=/work/Database/bwa-indexes/UCSC-$build
tmpdir=/tmp/custardPy.`date +%Y%m%d%H%M`
fastq_post="_"  # "_" or "_R"
enzyme=HindIII
ncore=64

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"
cell=Hap1-A
fqdir=fastq/$cell
odir=JuicerResults/$cell

rm -rf $odir
$sing juicer_map.sh -m $tmpdir -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

# Optional
$sing juicer_pigz.sh $odir
$sing plot_distance_count.sh $cell $odir
