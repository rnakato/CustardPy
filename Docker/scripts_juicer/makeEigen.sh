#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [options] <norm> <odir> <hic> <resolution> <genometable> <refFlat>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <hic>: .hic file' 1>&2
    echo '   <resolution>: resolution of matrix' 1>&2
    echo '   <genometable>: genometable file' 1>&2
    echo '   <refFlat>: gene annotation file (refFlat format)' 1>&2
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

if [ $# -ne 6 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
hic=$3
binsize=$4
gt=$5
gene=$6

pwd=$(cd $(dirname $0) && pwd)
chrlist=$(getchr_from_genometable.sh $gt)

XDG_RUNTIME_DIR=$HOME/.xdg

dir=$matrixdir/Eigen/$binsize
mkdir -p $dir

ex(){ echo $1; eval $1; }

getEigen(){
    h1d one PC1 $matrixdir/Matrix/intrachromosomal/$binsize/observed.$norm.chr$chr.matrix.gz $binsize chr$chr -p $dir/geneDensity.txt -o $dir/eigen.$norm.chr$chr
    
    if test -e $dir/eigen.$norm.chr$chr.bedGraph; then
        cut -f4 $dir/eigen.$norm.chr$chr.bedGraph | sed -e 's/^$/nan/g' > $dir/eigen.$norm.chr$chr.txt
        gzip -f $dir/eigen.$norm.chr$chr.txt
    fi
}
export -f getEigen

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

tmpfile=$(mktemp)
grep -v geneName $gene > $tmpfile
ex "h1d basic gd $tmpfile $binsize $gt -o $dir/geneDensity"
rm $tmpfile

func(){
    chr=$1
    dir=$2
    norm=$3
    binsize=$4
    matrixdir=$5
    if test $chr != "chrY" -a $chr != "chrM" -a $chr != "chrMT" ; then
       chr=$(echo $chr | sed -e 's/chr//g')

       getEigen

       if test -e $dir/eigen.$norm.chr$chr.txt.gz; then
           classifyCompartment.py $dir/eigen.$norm.chr$chr.txt.gz $dir/Compartment.$norm.chr$chr chr$chr $binsize
           toBed12 $dir/Compartment.$norm.chr$chr
       fi
    fi
}
export -f func

echo ${chrlist[@]} | tr ' ' '\n' | xargs -n1 -I {} -P $ncore bash -c "func {} $dir $norm $binsize $matrixdir"

for str in A B All StrongA WeakA WeakB StrongB
do
    if find $dir -maxdepth 1 -name "Compartment.SCALE.chr*.$str.bed" | read; then
        cat $dir/Compartment.$norm.chr*.$str.bed > $dir/Compartment.$norm.genome.$str.bed
    fi

    if find $dir -maxdepth 1 -name "Compartment.$norm.chr*.$str.bed12" | read; then
        cat $dir/Compartment.$norm.chr*.$str.bed12 > $dir/Compartment.$norm.genome.$str.bed12
    fi
done
