#!/bin/bash

build=hg38
gt=/work/Database/UCSC/$build/genome_table
bwaindex=/work/Database/bwa-indexes/UCSC-$build
gene=/work/Database/UCSC/$build/refFlat.txt
tmpdir=/tmp/juicer.1.6.1.`date +%Y%m%d%H%M`
fastq_post="_"  # "_" or "_R"
enzyme=MboI
ncore=64

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"
cell=Hap1-A
fqdir=$(pwd)/fastq/$cell/
odir=$(pwd)/JuicerResults/$cell

#rm -rf $odir
#$sing juicer_map.sh -m $tmpdir -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

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
