#!/bin/bash

build=hg38
gt=genometable.$build.txt
bwaindex=bwa-indexes/$build
tmpdir=/tmp/custardPy.`date +%Y%m%d%H%M`
fastq_post="_"  # "_" or "_R"
enzyme=MboI
ncore=64

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec custardpy.sif"

cell=Hap1-A
fqdir=fastq/$cell
odir=CustardPyResults_Hi-C/Juicer_$build/$cell

rm -rf $odir
$sing juicer_map.sh -m $tmpdir -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

# Optional
$sing juicer_pigz.sh $odir
$sing plot_distance_count.sh $cell $odir
