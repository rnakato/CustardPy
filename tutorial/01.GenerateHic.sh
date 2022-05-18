#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI

gt=/work/Database/UCSC/$build/genome_table
gene=/work/Database/UCSC/$build/refFlat.txt
sing="singularity exec --bind /work,/work2 /work/SingularityImages/rnakato_juicer.1.6.2.sif"
tmpdir=/tmp/juicer.1.6.1.`date +%Y%m%d%H%M`

cell=Hap1-A
fqdir=$(pwd)/fastq/$cell/
odir=$(pwd)/JuicerResults/$cell

rm -rf $odir
$sing juicer_map.sh -m $tmpdir $odir $build $fqdir $enzyme $fastq_post

# Optional
$sing juicer_pigz.sh $odir
$sing plot_distance_count.sh $cell $odir

exit
##### chromap
sing="singularity exec --bind /work,/work2 /work/SingularityImages/rnakato_4dn.img"
enzyme=MboI
enzymelen=4
restrictionsite=/work/Database/HiC-restriction_sites/${enzyme}_resfrag_$build.bed
gunzip -c $input_pairs | $sing /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - out.ff.pairs $restrictionsite
