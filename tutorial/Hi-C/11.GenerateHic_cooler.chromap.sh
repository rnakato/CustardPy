#!/bin/bash
build=hg38
ncore=64
chromap_index=chromap-indexes/$build
genome=genome.$build.fa
gt=genometable.$build.txt

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

fq1=fastq/Hap1-A/SRR5266584_1.fastq.gz
fq2=fastq/Hap1-A/SRR5266584_2.fastq.gz
odir=Cooler_$build/Hap1-A
#fq1=fastq/Control/SRR5952305_1.fastq.gz
#fq2=fastq/Control/SRR5952305_2.fastq.gz
#odir=Cooler_$build/Control
enzyme=HindIII
restrictionsite=/Cooler-restriction_sites/${enzyme}_resfrag_$build.bed

mkdir -p $odir/pairs $odir/log

echo "start mapping..."
$sing chromap --preset hic -t $ncore --remove-pcr-duplicates -x $chromap_index -r $genome -1 $fq1 -2 $fq2 -o $odir/pairs/mapped.rmdup.pairs 2> $odir/log/chromap.rmdup.log

exit
echo "start splitting pairsam by pairtools..."
TEMPFILE=$tempdir/temp.gz
TEMPFILE1=$tempdir/temp1.gz
# Select UU, UR, RU reads
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
          --output-rest $odir/pairs/bwa.unmapped.sam.pairs.gz \
          --output ${TEMPFILE} \
          $odir/pairs/bwa.sam.pairs.gz
pairtools split --output-pairs ${TEMPFILE1} ${TEMPFILE}
pairtools select 'True' --chrom-subset $gt -o $odir/pairs/bwa.dedup.pairs.gz ${TEMPFILE1}
pairix $odir/pairs/bwa.dedup.pairs.gz  # sanity check & indexing
rm ${TEMPFILE} ${TEMPFILE1}

    echo "add juicer-style fragment information..."
    # use fragment_4dnpairs.pl in pairix/util instead of juicer/CPU/common
    ffpairs=$odir/pairs/bwa.ff.pairs
    gunzip -c $odir/pairs/bwa.marked.sam.pairs.gz | fragment_4dnpairs.pl -a - $ffpairs $restrictionsite
    bgzip  -f $ffpairs
    pairix -f $ffpairs.gz

    rm -rf $tempdir



#$sing gunzip -c $odir/pairs/mapped.rmdup.pairs.gz | $sing fragment_4dnpairs.pl -a - $odir/pairs/mapped.ff.pairs $restrictionsite
#$sing bgzip  -f $odir/pairs/mapped.ff.pairs
#$sing pairix -f $odir/pairs/mapped.ff.pairs.gz
#pair=$odir/pairs/mapped.ff.pairs.gz
#$sing genCoolHiC_from_pairs Control $odir $pair $gt

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


gt=/work/Database/UCSC/$build/genome_table
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



exit
##### chromap
sing="singularity exec rnakato_4dn.img"
enzyme=MboI
enzymelen=4
restrictionsite=/work/Database/HiC-restriction_sites/${enzyme}_resfrag_$build.bed
gunzip -c $input_pairs | $sing /usr/local/bin/pairix/util/fragment_4dnpairs.pl -a - out.ff.pairs $restrictionsite
