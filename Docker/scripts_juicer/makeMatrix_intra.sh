#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [-l] <norm> <odir> <hic> <resolution> <gt>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <hic>: .hic file' 1>&2
    echo '   <resolution>: resolution of the matrix' 1>&2
    echo '   <gt>: genome table' 1>&2
    echo '   Options:' 1>&2
    echo '     -l: output contact matrix as a list (default: dense matrix)' 1>&2
}

list="no"
while getopts l option; do
    case ${option} in
        l) list="yes" ;;
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

norm=$1
matrixdir=$2
hic=$3
binsize=$4
gt=$5

pwd=$(cd $(dirname $0) && pwd)
chrlist=$($pwd/getchr_from_genometable.sh $gt)

dir=$matrixdir/Matrix/intrachromosomal/$binsize
mkdir -p $dir
for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT"; then continue; fi

    echo "$chr"
    for type in observed oe
    do
	tempfile=$dir/$type.$norm.$chr.txt
        juicertools.sh dump $type $norm $hic $chr $chr BP $binsize $tempfile
	if test $list = "no" -o -s $tempfile; then
            convert_JuicerDump_to_dense.py $tempfile $dir/$type.$norm.$chr.matrix.gz $gt $chr $chr -r $binsize
	    rm $tempfile
	fi
    done
    #    for type in expected norm
    #    do
    #        juicertools.sh dump $type $norm $hic.hic $chr BP $binsize $dir/$type.$norm.$chr.matrix -d
    #    done
done
