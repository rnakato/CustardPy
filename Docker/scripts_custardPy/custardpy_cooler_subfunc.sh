#!/bin/bash

ex(){ echo $1; eval $1; }

parse_pairtools(){
    odir=$1
    gt=$2
    qthre=$3
    restrictionsite=$4
    enzymelen=$5

    bamdir=$odir/mapfile
    pairdir=$odir/pairs
    logdir=$odir/log
    tempdir=$odir/temp
    qcdir=$odir/qc_report/
    mkdir -p $bamdir $pairdir $logdir $tempdir $qcdir

    echo "start parsing by pairtools..."
    # '--walks-policy mask' rescues single ligations: https://github.com/open2c/pairtools/blob/master/doc/parsing.rst
    # '--no-flip' option of pairtools parse causes an error "RuntimeError: Pairs are not triangular: found blocks 'chr2|chr11'' and 'chr11|chr2'"
    samtools view -h $bamdir/mapped.bwa.cram \
    | pairtools parse --chroms-path $gt --drop-sam --drop-seq \
      --output-stats $qcdir/pairtools_parse.stats.txt  \
      --walks-policy mask --add-columns mapq --min-mapq $qthre --max-inter-align-gap 30 \
    | pairtools sort --nproc 5 --tmpdir=$tempdir --output $pairdir/mapped.bwa.q$qthre.pairs.gz

    echo "pairtools dedup..."
    pairtools dedup --mark-dups --output-stats $qcdir/pairtools_dedup.stats.txt \
            --output-dups $pairdir/dup.bwa.q$qthre.pairs.gz --output-unmapped $pairdir/unmapped.bwa.q$qthre.pairs.gz \
            --output $pairdir/dedup.bwa.q$qthre.pairs.gz $pairdir/mapped.bwa.q$qthre.pairs.gz
    pairix -f $pairdir/dedup.bwa.q$qthre.pairs.gz # sanity check

    get_qc.py -p $qcdir/pairtools_dedup.stats.txt > $qcdir/mapping_stats.txt

    # Note that (for now) the pairtools module for MultiQC is only available in the open2C fork and not in the main MultiQC repository.
#    pairtools stats $pairdir/dedup.bwa.q$qthre.pairs.gz -o $qcdir/pairtools.stats.txt
#    multiqc -o qcdir/ $qcdir/pairtools.stats.txt 
#    rm $qcdir/pairtools.stats.txt

#    echo "start splitting pairsam by pairtools..."
#    TEMPFILE=$tempdir/temp.gz
#    TEMPFILE1=$tempdir/temp1.gz
    ## Select UU, UR, RU reads
#    pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
#            --output-rest $pairdir/bwa.unmapped.sam.pairs.gz \
#            --output ${TEMPFILE} \
#            $pairdir/bwa.dedup.sam.pairs.gz
#    pairtools split --output-pairs ${TEMPFILE1} ${TEMPFILE}
#    pairtools select 'True' --chrom-subset $gt -o $pairdir/bwa.marked.sam.pairs.gz ${TEMPFILE1}
#    pairix $pairdir/bwa.marked.sam.pairs.gz # sanity check & indexing
#    rm ${TEMPFILE} ${TEMPFILE1} $pairdir/bwa.sam.pairs.gz $pairdir/bwa.dedup.sam.pairs.gz

    if test "$restrictionsite" != ""; then
        echo "add restrictionsite information..."
    #    ffpairs=$pairdir/bwa.ff.pairs
    #    gunzip -c $pairdir/dedup.bwa.q$qthre.pairs.gz | fragment_4dnpairs.pl -a - $ffpairs $restrictionsite
    #    bgzip  -f $ffpairs
    #    pairix -f $ffpairs.gz
        pairs=$pairdir/dedup.restricted.q$qthre.pairs.gz
        pairtools restrict -f $restrictionsite $pairdir/dedup.bwa.q$qthre.pairs.gz -o $pairs
        pairix -f $pairs

    #   echo "pairsqc.py -p $pair -c $gt -tP -s $prefix -O $odir/qc -M $max_distance"
        echo "Implement pairsqc..."
        python /opt/pairsqc/pairsqc.py -p $pairs -c $gt -tP -s $prefix -O $odir/qc -M $max_distance
        Rscript /opt/pairsqc/plot.r $enzymelen $odir/qc_report
    else
        pairs=$pairdir/dedup.bwa.q$qthre.pairs.gz
    fi
    rm -rf $tempdir

    echo "pairtools finished!"
    echo "Output pairs file: $pairs"
}

parse_pairtools_MicroC_old(){
    odir=$1
    gt=$2
    qthre=$3

    bamdir=$odir/mapfile
    pairdir=$odir/pairs
    logdir=$odir/log
    tempdir=$odir/temp
    qcdir=$odir/qc_report/
    mkdir -p $bamdir $pairdir $logdir $tempdir $qcdir

    echo "start parsing by pairtools..."
    samtools view -h $bamdir/mapped.bwa.cram \
    | pairtools parse --chroms-path $gt --output-stats $qcdir/pairtools_parse.stats.txt \
                      --walks-policy mask --add-columns mapq --min-mapq $qthre --max-inter-align-gap 30 \
    | pairtools sort --nproc 4 --tmpdir=$tempdir --output $pairdir/bwa.sam.pairs.gz

    echo "pairtools dedup..."
    pairtools dedup --nproc-in 4 --nproc-out 4 --mark-dups --output-dups - --output-unmapped - \
                --output-stats $qcdir/mapped.bwa.stats.txt \
                --output $pairdir/mapped.bwa.dedup.pairs.gz $pairdir/bwa.sam.pairs.gz
    get_qc.py -p $qcdir/mapped.bwa.stats.txt > $qcdir/mappingstats.txt
    pairix -f $pairdir/mapped.bwa.dedup.pairs.gz # sanity check

    echo "start splitting pairsam by pairtools..."
    ## Select UU, UR, RU reads
    TEMPFILE=$tempdir/temp.gz
    TEMPFILE1=$tempdir/temp1.gz
    pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
            --output-rest $pairdir/bwa.unmapped.sam.pairs.gz \
            --output ${TEMPFILE} \
            $pairdir/mapped.bwa.dedup.pairs.gz
    pairtools split --nproc-in 4 --nproc-out 4 --output-pairs ${TEMPFILE1} --output-sam - ${TEMPFILE} \
           | samtools view -C - -T $genome \
           | samtools sort -m 8G -O cram \
           >  $bamdir/mapped.final.sort.cram
    samtools index $bamdir/mapped.final.sort.cram

    pairtools select 'True' --chrom-subset $gt -o $pairdir/bwa.marked.sam.pairs.gz ${TEMPFILE1}
    pairix $pairdir/bwa.marked.sam.pairs.gz # sanity check & indexing
    rm ${TEMPFILE} ${TEMPFILE1} $pairdir/bwa.sam.pairs.gz $pairdir/mapped.bwa.dedup.pairs.gz
    rm -rf $tempdir

    echo "pairtools finished!"
    echo "Output pairs file: $pairdir/bwa.marked.sam.pairs.gz"
}

gen_cool_hic(){
    odir=$1
    gt=$2
    prefix=$3
    binsize_multi=$4
    max_distance=$5
    max_split=$6
    qthre=$7
    pair=$8

    echo "generate .cool file..."
    mkdir -p $odir/cool $odir/log
    cooler cload pairix -p $ncore -s $max_split $gt:$binsize_min $pair $odir/cool/$prefix.cool >$odir/log/cooler_cload_pairix.log
    cooler balance -p $ncore $odir/cool/$prefix.cool >$odir/log/cooler_balance.log

    for binsize in 25000 50000 100000; do
        cfile=$odir/cool/$prefix.$binsize.cool
        cooler cload pairix -p $ncore -s $max_split $gt:$binsize $pair $cfile >$odir/log/cooler_cload_pairix.$binsize.log
        cooler balance -p $ncore $cfile >$odir/log/cooler_balance.$binsize.log
    done
    run-cool2multirescool.sh -i $odir/cool/$prefix.cool -p $ncore -o $odir/cool/$prefix -u $binsize_multi

    echo "generate .hic..."
    mkdir -p $odir/hic
    juicertools.sh pre -q $qthre $pair $odir/hic/contact_map.q$qthre.hic $gt

    echo "postprocess finished!"
    echo "Output pairs file: $odir/cool/$prefix.cool and $odir/hic/contact_map.q$qthre.hic"
}
