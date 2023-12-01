#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <odir> <hic>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <hic>: .hic file' 1>&2
    echo '   Options:' 1>&2
    echo '     -r resolutions: the resolutions (default: "5000,10000,25000", should be quoted and separated by comma)' 1>&2
    echo '     -o: Use older version of juicer_tools.jar (juicer_tools.1.9.9_jcuda.0.8.jar, default: juicer_tools.1.22.01.jar)' 1>&2
}

resolutions="5000,10000,25000"
useoldversion="no"
while getopts r:o option; do
    case ${option} in
        r) resolutions=${OPTARG} ;;
        o) useoldversion="yes" ;;
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

if [ $# -ne 3 ]; then
  usage
  exit 1
fi

norm=$1
odir=$2
hic=$3

ex(){ echo $1; eval $1; }

pwd=$(cd $(dirname $0) && pwd)

hicdir=$odir/loops/$norm
mkdir -p $hicdir

if test $useoldversion = "no"; then
    # juicer_tools.1.22.01.jar
    ex "juicertools.sh hiccups -r $resolutions -k $norm $hic $hicdir --ignore-sparsity"
else 
    ex "juicertools.sh -o hiccups -r $resolutions -k $norm $hic $hicdir --ignore_sparsity"
    # add "chr" to the 1st and 4th columns
    mv $hicdir/merged_loops.bedpe $hicdir/merged_loops.bedpe.original
    cat $hicdir/merged_loops.bedpe.original | awk -F'\t' 'BEGIN {OFS="\t"} {if ($0 ~ /^#/) print; else {$1="chr"$1; $4="chr"$4; print}}'  > $hicdir/merged_loops.bedpe
fi

grep -v \# $hicdir/merged_loops.bedpe | awk '{OFS="\t"} {printf "%s\t%d\t%d\t%s\t%d\t%d\n", $1, $2, $3, $4, $5, $6 }' > $hicdir/merged_loops.simple.bedpe