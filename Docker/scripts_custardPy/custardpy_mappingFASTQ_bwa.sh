#!/bin/bash

ex(){ echo $1; eval $1; }

mapping_reads_bwa(){
    dir=$1
    odir=$2
    prefix=$3
    index_bwa=$4
    ncore=$5
    postfix1=$6
    postfix2=$7

    if [ ! -d "$dir" ]; then
        echo "Error: $dir does not exist."
        return
    fi
    if test "$index_bwa" = ""; then
        echo "Error: specify BWA index (-i)."
        exit 1
    fi

    echo "start mapping by BWA..."

    fq1_list=()
    fq2_list=()

    for fq1 in "$dir"/*"$postfix1"; do
        if [ ! -f "$fq1" ]; then
            echo "$dir: *$postfix1 does not exist. Skipping"
            continue
        fi
        fq2="${fq1%$postfix1}$postfix2"
        if [ ! -f "$fq2" ]; then
            echo "$dir: FASTQ_2 "$fq2" does not exist. Skipping"
            continue
        fi

        fq1_list+=("$fq1")
        fq2_list+=("$fq2")
    done

    bamdir=$odir/mapfile
    logdir=$odir/log
    tempdir=$odir/temp
    mkdir -p $bamdir $logdir $tempdir

    bams=""
    for ((i=0; i<${#fq1_list[@]}; i++)); do
        fq1=${fq1_list[$i]}
        fq2=${fq2_list[$i]}
        name=`basename $fq1 $postfix1`

        ex "bwa mem -t $ncore -T0 -SP $index_bwa $fq1 $fq2 2> $logdir/bwa_mapping_sam_$name > $tempdir/$name.bwa.sam"

        bams+=" $tempdir/$name.bwa.sam"
    done
    if test "$bams" != ""; then
        echo "Merge SAM files to BAM..."
        ex "samtools merge -f -@ $ncore - $bams | samtools view -C - -T $genome > $bamdir/mapped.bwa.cram"
#        ex "samtools merge -f -@ $ncore $bamdir/mapped.bwa.bam $bams"
        ex "rm $bams"
    else
        echo "Error: no FASTQ files are available."
        echo "Check the FASTQ directory ($fqdir/*) and the postfix of FASTQ files (*$postfix1 and *$postfix2)."
        return
    fi
    rm -rf $tempdir

    echo "mapping finished!"
    echo "Output file: $bamdir/mapped.bwa.cram"
}


mapping_reads_chromap(){
    dir=$1
    odir=$2
    prefix=$3
    index_chromap=$4
    genome=$5
    ncore=$6
    postfix1=$7
    postfix2=$8

    if [ ! -d "$dir" ]; then
        echo "Error: $dir does not exist."
        return
    fi
    if test "$index_chromap" = ""; then
        echo "Error: specify chromap index (-i)."
        exit 1
    fi
    if test "$genome" = ""; then
        echo "Error: specify genome file (-f)."
        exit 0
    fi

    echo "start mapping by chromap..."
    fq1_list=()
    fq2_list=()

    for fq1 in "$dir"/*"$postfix1"; do
        if [ ! -f "$fq1" ]; then
            echo "$dir: *$postfix1 does not exist. Skipping"
            continue
        fi
        fq2="${fq1%$postfix1}$postfix2"
        if [ ! -f "$fq2" ]; then
            echo "$dir: FASTQ_2 "$fq2" does not exist. Skipping"
            continue
        fi

        fq1_list+=("$fq1")
        fq2_list+=("$fq2")
    done
    fq1_list=$(IFS=","; echo "${fq1_list[*]}")
    fq2_list=$(IFS=","; echo "${fq2_list[*]}")

    pairdir=$odir/pairs
    logdir=$odir/log
    mkdir -p $pairdir $logdir

    chromap --preset hic -t $ncore --remove-pcr-duplicates -x $index_chromap -r $genome \
            -1 $fq1_list -2 $fq2_list -o $pairdir/mapped.chromap.rmdup.pairs \
            2> $logdir/chromap.rmdup.log
            
    bgzip -f $pairdir/mapped.chromap.rmdup.pairs
    pair=$pairdir/mapped.chromap.rmdup.pairs.gz
    pairix $pair # sanity check

    echo "mapping finished!"
    echo "Output file: $pair"
}

