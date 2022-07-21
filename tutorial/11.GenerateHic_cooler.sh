#!/bin/bash
build=hg38
gt=/work/Database/UCSC/$build/genome_table
index_bwa=/work/Database/bwa-indexes/UCSC-$build
ncore=64

#fq1=fastq/Hap1-A/SRR5266584_1.fastq.gz
#fq2=fastq/Hap1-A/SRR5266584_2.fastq.gz
#idir=Cooler_$build/Hap1-A
fq1=fastq/Control/SRR5952305_1.fastq.gz
fq2=fastq/Control/SRR5952305_2.fastq.gz
odir=Cooler_$build/Control

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"
ncore=64
Ddir=/work/Database/Database_fromDocker/Ensembl-GRCh38
chromap_index=$Ddir/chromap-indexes/genome
genome=$Ddir/genome.fa

mkdir -p $odir/pairs $odir/log

echo "start mapping..."
#$sing chromap --preset hic -t $ncore --remove-pcr-duplicates -x $chromap_index -r $genome -1 $fq1 -2 $fq2 -o $odir/pairs/mapped.rmdup.pairs 2> $odir/log/chromap.rmdup.log
echo "parse by pairtools..."
enzyme=HindIII
restrictionsite=/Cooler-restriction_sites/${enzyme}_resfrag_$build.bed

$sing gunzip -c $odir/pairs/mapped.rmdup.pairs.gz | $sing fragment_4dnpairs.pl -a - $odir/pairs/mapped.ff.pairs $restrictionsite
$sing bgzip  -f $odir/pairs/mapped.ff.pairs
$sing pairix -f $odir/pairs/mapped.ff.pairs.gz
pair=$odir/pairs/mapped.ff.pairs.gz
$sing genCoolHiC_from_pairs Control $odir $pair $gt

exit

enzyme=HindIII
mkdir -p $idir/hic $idir/cool

$sing bgzip -d -c $input_pairs | $sing fragment_4dnpairs.pl -a - $idir/pairs/mapped.bwa.ff.pairs $restrictionsite
$sing fragment_4dnpairs.pl -a $idir/pairs/mapped.bwa.pairs $idir/pairs/mapped.bwa.ff.pairs $restrictionsite
$sing bgzip $idir/pairs/mapped.bwa.ff.pairs
pair=$idir/pairs/mapped.bwa.ff.pairs.gz
$sing pairix -p pairs $pair

max_split=2
binsize_min=5000
binsize_multi="5000,10000,25000,50000,100000,500000,1000000,2500000,5000000,10000000"
tool=chromap
$sing python /opt/scripts/pairsqc.py -p $pair -c $gt -tP -s $tool -O $idir/qc.$tool -M 8.4
$sing Rscript plot.r 4 $idir/qc.${tool}_report
$sing cooler cload pairix -p $ncore -s $max_split $gt:$binsize_min $pair $idir/cool/$prefix.$tool.cool
$sing cooler balance -p $ncore $idir/cool/$prefix.$tool.cool
$sing run-cool2multirescool.sh -i $idir/cool/$prefix.$tool.cool -p $ncore -o $idir/cool/$prefix.$tool -u $binsize_multi


exit
### bwa
mkdir -p $idir/bam $idir/pairs $idir/temp
$sing bwa mem -5SP -T0 -t $ncore $index_bwa $fq1 $fq2 | samtools sort -@2 > $idir/bam/mapped.bwa.sort.bam

$sing pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path $gt $idir/bam/mapped.bwa.sort.bam \
       | pairtools sort --nproc 8 --tmpdir=$idir/temp \
       | pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats $idir/pairs/mapped.bwa.stats.txt \
       | pairtools split --nproc-in 8 --nproc-out 8 --output-pairs $idir/pairs/mapped.bwa.pairs --output-sam - \
       | samtools sort -@16 -T $idir/temp/temp.bam > $idir/bam/mapped.final.sort.bam
rm -rf $idir/temp

echo "samtools indexing..."
samtools index $idir/bam/mapped.final.sort.bam
python /work/sakata/python/get_qc.py -p $idir/pairs/mapped.bwa.stats.txt > $idir/pairs.stats.txt

mkdir -p $idir/hic $idir/cool
pair=$idir/pairs/mapped.bwa.pairs.gz
$sing pairix -p pairs $pair

#$sing genCoolHiC_from_pairs $idir/hic/q30 $idir $pair $gt


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
