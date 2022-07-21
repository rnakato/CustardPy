#!/bin/bash

build=hg38
gt=genometable.$build.txt
bwaindex=/work/Database/bwa-indexes/UCSC-$build
tmpdir=/tmp/custardPy.`date +%Y%m%d%H%M`
fastq_post="_"  # "_" or "_R"
enzyme=HindIII
ncore=64

sing="singularity exec --bind /mnt/d custardpy_juicer.sif"
cell=Hap1-A
#fqdir=$(pwd)/fastq/$cell/
#odir=$(pwd)/JuicerResults/$cell
fqdir=fastq/$cell/
odir=JuicerResults/$cell

rm -rf $odir
$sing juicer_map.sh -m $tmpdir -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post
exit
# Optional
$sing juicer_pigz.sh $odir
$sing plot_distance_count.sh $cell $odir



exit
##### chromap
sing="singularity exec rnakato_4dn.img"
enzyme=MboI
enzymelen=4
restrictionsite=/work/Database/HiC-restriction_sites/${enzyme}_resfrag_$build.bed
gunzip -c $input_pairs | $sing /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - out.ff.pairs $restrictionsite
