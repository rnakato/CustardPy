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
}

resolutions="5000,10000,25000"
while getopts r: option; do
    case ${option} in
        r) resolutions=${OPTARG} ;;
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
juicertool="juicertools.sh"

hicdir=$odir/loops/$norm
mkdir -p $hicdir
ex "$juicertool hiccups -r $resolutions -k $norm $hic $hicdir --ignore-sparsity"
grep -v \# $hicdir/merged_loops.bedpe | awk '{OFS="\t"} {printf "%s\t%d\t%d\t%s\t%d\t%d\n", $1, $2, $3, $4, $5, $6 }' > $hicdir/merged_loops.simple.bedpe
