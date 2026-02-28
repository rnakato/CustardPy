#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [options] <odir> <cool> <resolution> <genometable> <genome>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <cool>: .cool file' 1>&2
    echo '   <resolution>: resolution of matrix' 1>&2
    echo '   <genometable>: genometable file' 1>&2
    echo '   <genome>: genome file (fasta format)' 1>&2
    echo '   Options:' 1>&2
    echo '     -p <int>: the number of CPUs (default: 6)' 1>&2
}

ncore=6
while getopts p: option; do
    case ${option} in
        p) ncore=${OPTARG} ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -ne 5 ]; then
  usage
  exit 1
fi

odir=$1
cool=$2
binsize=$3
gt=$4
genome=$5

pwd=$(cd $(dirname $0) && pwd)
chrlist=$(getchr_from_genometable.sh $gt)

ex(){ echo $1; eval $1; }

mkdir -p $odir
cooltools genome binnify $gt $binsize > $odir/bins.$binsize.bed
cooltools genome gc $odir/bins.$binsize.bed $genome > $odir/gc.$binsize.bedGraph

cooltools eigs-cis $cool::/resolutions/$binsize \
    --phasing-track $odir/gc.$binsize.bedGraph::GC \
    --n-eigs 3 \
    -o $odir/eigs

pigz -f $odir/eigs.cis.vecs.tsv

for chr in ${chrlist[@]}
do
    zcat $odir/eigs.cis.vecs.tsv | awk -v chr="$chr" '$1==chr { print ($5 == "" ? "nan" : $5) }' > $odir/eigs.cis.$chr.txt
    pigz -f $odir/eigs.cis.$chr.txt
done

cooltools_eigs.py $cool::/resolutions/$binsize $odir/eigs.cis.vecs.tsv.gz

toBed12(){
    prefix=$1
    cat $prefix.StrongA.bed | awk -v 'OFS=\t' '{print $1, $2, $3, "StrongA\t0\t+", $2, $3, "255,0,50"}' > $prefix.StrongA.bed12
    cat $prefix.WeakA.bed   | awk -v 'OFS=\t' '{print $1, $2, $3, "WeakA\t0\t+", $2, $3, "255,255,50"}' > $prefix.WeakA.bed12
    cat $prefix.WeakB.bed   | awk -v 'OFS=\t' '{print $1, $2, $3, "WeakB\t0\t+", $2, $3, "50,150,50"}' > $prefix.WeakB.bed12
    cat $prefix.StrongB.bed | awk -v 'OFS=\t' '{print $1, $2, $3, "StrongB\t0\t+", $2, $3, "50,50,255"}' > $prefix.StrongB.bed12

    cat $prefix.StrongA.bed12 \
        $prefix.WeakA.bed12 \
        | sort -k1,1 -k2,2n \
        > $prefix.A.bed12
    cat $prefix.WeakB.bed12 \
        $prefix.StrongB.bed12 \
        | sort -k1,1 -k2,2n \
        > $prefix.B.bed12
    cat $prefix.StrongA.bed12 \
        $prefix.WeakA.bed12 \
        $prefix.WeakB.bed12 \
        $prefix.StrongB.bed12 \
        | sort -k1,1 -k2,2n \
        > $prefix.All.bed12
}
export -f toBed12

getCompartment(){
    chr=$1
    dir=$2
    binsize=$3
    if test $chr != "chrY" -a $chr != "chrM" -a $chr != "chrMT" ; then
        if test -e $dir/eigs.cis.$chr.txt.gz; then
            classifyCompartment.py $dir/eigs.cis.$chr.txt.gz $dir/PC1.cis.$chr $chr $binsize
            toBed12 $dir/PC1.cis.$chr
       fi
    fi
}
export -f getCompartment

echo ${chrlist[@]} | tr ' ' '\n' | xargs -n1 -I {} -P $ncore bash -c "getCompartment {} $odir $binsize"

for str in A B All StrongA WeakA WeakB StrongB
do
    if find $odir -maxdepth 1 -name "Compartment.SCALE.chr*.$str.bed" | read; then
        cat $odir/PC1.cis.chr*.$str.bed > $odir/PC1.cis.genome.$str.bed
    fi
    rm $odir/PC1.cis.chr*.$str.bed

    if find $odir -maxdepth 1 -name "PC1.cis.chr*.$str.bed12" | read; then
        cat $odir/PC1.cis.chr*.$str.bed12 > $odir/PC1.cis.genome.$str.bed12
    fi
    rm $odir/PC1.cis.chr*.$str.bed12 
done

#cat $odir/PC1.cis.chr*.Weak.bed > $odir/PC1.cis.genome.Weak.bed
#rm $odir/PC1.cis.chr*.Weak.bed