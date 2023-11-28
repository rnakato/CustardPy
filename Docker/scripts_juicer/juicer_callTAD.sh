#!/bin/bash -e
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [options] <norm> <odir> <hic> <gt>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <hic>: .hic file' 1>&2
    echo '   <gt>: genome table' 1>&2
    echo '   Options:' 1>&2
    echo '     -r resolutions: the resolutions for ArrowHead (default: "10000 25000 50000", should be quoted and separated by spaces)' 1>&2
    echo '     -p ncore: number of CPUs (default: 24)' 1>&2
    echo '     -o: Use older version of juicer_tools.jar (juicer_tools.1.9.9_jcuda.0.8.jar, default: juicer_tools.1.22.01.jar)' 1>&2
}

resolutions="10000 25000 50000"
ncore=24
useoldversion="no"
while getopts r:p:o option; do
    case ${option} in
        r) resolutions=${OPTARG} ;;
        p) ncore=${OPTARG} ;;
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

if [ $# -ne 4 ]; then
  usage
  exit 1
fi

norm=$1
odir=$2
hic=$3
gt=$4

dir=$odir/TAD/$norm
mkdir -p $dir

for res in $resolutions; do
    if test ! -e $dir/${res}_blocks.bedpe; then
        if test $useoldversion = "no"; then
            # juicer_tools.1.22.01.jar
            juicertools.sh arrowhead -m 2000 -r $res --threads $ncore -k $norm --ignore-sparsity $hic $dir 
        else 
            juicertools.sh -o arrowhead -m 2000 -r $res --threads $ncore -k $norm --ignore_sparsity $hic $dir 
        fi
    fi

    # make TAD bed
    grep -v \# $dir/${res}_blocks.bedpe | awk '{OFS="\t"} {printf "%s\t%d\t%d\tTAD%d\n", $1, $2, $3, NR }' > $dir/${res}_blocks.bed
    sort -k1,1 -k2,2n $dir/${res}_blocks.bed | bedtools merge > $dir/${res}_blocks.merged.bed

    # bedpe to boundary bed
    grep -v \# $dir/${res}_blocks.bedpe \
        | awk -v res=$res '{OFS="\t"} NR>1 {printf "%s\t%d\t%d\n%s\t%d\t%d\n", $1, $2, $2+res, $1, $3-res, $3}' \
        | sort -k1,1 -k2,2n \
        | uniq \
        | sort -k1,1 -k2,2n \
        | bedtools merge > $dir/${res}_blocks.boundaries.bed

    # TAD coverage
    bedtools genomecov -bga -i $dir/${res}_blocks.bed -g $gt > $dir/${res}_blocks.TADcoverage.bed

    # intra-TAD regions
    cat $dir/${res}_blocks.TADcoverage.bed \
        | awk '{if($4>0 && $3-$2>=100000) printf "%s\t%d\t%d\t%d\n", $1, $2,$3, $3-$2}' \
              >  $dir/${res}_blocks.TADregions.bed

    # non-TAD region
    cat $dir/${res}_blocks.TADcoverage.bed | awk '{if($4==0) print}' | cut -f1-3 | sort -k1,1 -k2,2n > $dir/${res}_blocks.nonTADregions.bed
done
