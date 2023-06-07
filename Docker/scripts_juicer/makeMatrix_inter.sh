#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [-l] <norm> <odir> <hic> <resolution> <chr1> <chr2>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <hic>: .hic file' 1>&2
    echo '   <resolution>: resolution of the matrix' 1>&2
    echo '   <chr1, chr2>: two input chromosomes' 1>&2
    echo '   Options:' 1>&2
    echo '     -l: output contact matrix as a list (default: dense matrix)' 1>&2
}

list="no"
while getopts l option; do
    case ${option} in
        l) list="yes" ;;
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
chr1=$5
chr6=$6

pwd=$(cd $(dirname $0) && pwd)

echo "$chr1-$chr2"
dir=$matrixdir/Matrix/interchromosomal/$binsize/$chr1-$chr2
mkdir -p $dir

for type in observed oe
do
    tempfile=$dir/$type.$norm.txt
    juicertools.sh dump $type $norm $hic $chr1 $chr2 BP $binsize $tempfile
    if test $list = "no" -o -s $tempfile; then
        convert_JuicerDump_to_dense.py $tempfile $dir/$type.$norm.matrix.gz $gt $chr1 $chr2 -r $binsize
        rm $tempfile
    fi
done

#for str in observed #oe
#do
#    merge_JuicerMatrix_to_Genome.py $dir/interchromosomal \
#				    $dir/interchromosomal/$binsize/genome.$str.full.$lim_pzero.pickle \
#				    $binsize $str $lim_pzero $chrnum
 #   merge_JuicerMatrix_to_Genome.py $dir/interchromosomal \
#				    $dir/interchromosomal/$binsize/genome.$str.evenodd.$lim_pzero.pickle \
#				    $binsize $str $lim_pzero $chrnum --evenodd
#done
