#!/bin/bash
build=hg38
gt=/work/Database/UCSC/$build/genome_table
index_bwa=/work/Database/bwa-indexes/UCSC-$build
ncore=64

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

fq1=fastq/Hap1-A/SRR5266584_1.fastq.gz
fq2=fastq/Hap1-A/SRR5266584_2.fastq.gz
idir=Cooler_$build/Hap1-A
mkdir -p $idir/bam $idir/pairs $idir/log $idir/temp
#$sing bwa mem -5SP -T0 -t $ncore $index_bwa $fq1 $fq2 | samtools sort -@2 > $idir/bam/mapped.bwa.sort.bam 2> $idir/log/bwa.log

#$sing pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path $gt $idir/bam/mapped.bwa.sort.bam \
#       | pairtools sort --nproc 8 --tmpdir=$idir/temp \
#       | pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats $idir/pairs/mapped.bwa.stats.txt \
#       | pairtools split --nproc-in 8 --nproc-out 8 --output-pairs $idir/pairs/mapped.bwa.pairs --output-sam - \
#       | samtools sort -@16 -T $idir/temp/temp.bam > $idir/bam/mapped.final.sort.bam
#rm -rf $idir/temp

#echo "samtools indexing..."
#samtools index $idir/bam/mapped.final.sort.bam
#python /work/sakata/python/get_qc.py -p $idir/pairs/mapped.bwa.stats.txt > $idir/pairs.stats.txt

$sing bgzip $idir/pairs/mapped.bwa.pairs
mkdir -p $idir/hic $idir/cool
pair=$idir/pairs/mapped.bwa.pairs.gz
$sing pairix -p pairs $pair
$sing gunzip -c $input_pairs | $sing /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - out.ff.pairs $restrictionsite

exit

gt=/work/Database/UCSC/$build/genome_table
enzyme=MboI
enzymelen=4
restrictionsite=/work/Database/HiC-restriction_sites/${enzyme}_resfrag_$build.bed

####$sing pairtools sort --nproc 32 --output sorted.pair.gz CTCFKD_1.pairs

#$sing pairtools dedup --mark-dups --output-dups --output-unmapped --output-stats stats.txt --output marked.gz CTCFKD_1.pairs
#$sing pairix -p pairs marked.gz
#$sing python3 /usr/local/bin/pairsqc/pairsqc.py -p marked.gz -c $gt -tP -s test -O pairqc -M 8.4

#$sing Rscript /usr/local/bin/pairsqc/plot.r $enzymelen $(pwd)/pairqc_report
input_pairs=marked.gz
#gunzip -c $input_pairs | $sing /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - out.ff.pairs $restrictionsite
#$sing bgzip  -f out.ff.pairs
#$sing pairix -f out.ff.pairs.gz
