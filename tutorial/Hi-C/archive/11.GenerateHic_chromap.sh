#!/bin/bash
build=hg38
ncore=64
chromap_index=chromap-indexes/$build
genome=genome.$build.fa
gt=genometable.$build.txt
gene=refFlat.$build.txt

#sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.2.2.2.sif"
sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.0.0.sif"

#sing="singularity exec custardpy.sif"

cell=Control
fqdir=fastq/$cell
odir=chromap_results_$build/$cell

fq1=$fqdir/SRR17870718_1.fastq.gz
fq2=$fqdir/SRR17870718_2.fastq.gz

enzyme=MboI
restrictionsite=/Cooler-restriction_sites/${enzyme}_resfrag_$build.bed

mkdir -p $odir/pairs $odir/log $odir/cool $odir/hic

echo "start mapping..."
#$sing chromap --preset hic -t $ncore --remove-pcr-duplicates -x $chromap_index -r $genome -1 $fq1 -2 $fq2 -o $odir/pairs/mapped.rmdup.pairs 2> $odir/log/chromap.rmdup.log
unpigz $odir/pairs/mapped.rmdup.pairs.gz
sort -k2,2d -k4,4d $odir/pairs/mapped.rmdup.pairs > $odir/pairs/mapped.sorted.pairs

#$sing pigz $odir/pairs/mapped.rmdup.pairs

exit

pairs=$odir/pairs/mapped.rmdup.pairs.gz
echo "start generating .cool..."
#$sing hictk load --format 4dn --bin-size 5kbp $pairs $odir/cool/mapped.5000.cool -t 4 --force --chrom-sizes $gt
#$sing hictk zoomify $odir/cool/mapped.5000.cool $odir/cool/mapped.mcool -t 4 --force
echo "start generating .hic..."
#$sing hictk load --format 4dn --bin-size 5kbp $pairs $odir/hic/mapped.5000.hic -t 4 --force --chrom-sizes $gt 
#$sing hictk zoomify $odir/hic/mapped.5000.hic $odir/hic/mapped.hic -t 4 --force

#exit

echo "Applying custardpy_process_hic..."
norm=SCALE
hic=$odir/hic/mapped.hic
dir=$odir/TAD/$norm
mkdir -p $dir
$sing custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $odir/hic/mapped.hic $odir
exit

#java -Xms512m -Xmx64384m -jar ../../Docker/juicer_tools.2.20.00.jar dump observed SCALE $hic chr22 chr22 BP 50000 testmatrix.txt"
exit

resolutions="10000 25000 50000"
for res in 50000 #$resolutions;
do
$sing java -Xms512m -Xmx64384m -jar ../../Docker/juicer_tools.2.20.00.jar arrowhead -m 2000 -r $res --threads $ncore -k $norm $hic $dir
done
exit
hicdir=$odir/loops/$norm
mkdir -p $hicdir
resolutions="5000,10000,25000"
$sing java -Xms512m -Xmx64384m -jar ../../Docker/juicer_tools.2.20.00.jar hiccups -r $resolutions -k $norm $hic $hicdir


exit 

singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.2.2.2.sif hictk load --format 4dn --bin-size 5kbp CustardPyResults_Hi-C/Cooler_hg38/Control/pairs/dedup.bwa.q30.pairs.gz cooler.mapped.5000.hic -t 4 --force --chrom-sizes genometable.hg38.txt --assume-unsorted

singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.2.2.2.sif hictk zoomify cooler.mapped.5000.hic cooler.mapped.hic --force

java -Xms512m -Xmx64384m -jar ../../Docker/juicer_tools.2.20.00.jar arrowhead -m 2000 -r 50000 --threads 64 -k SCALE cooler.mapped.hic testfff