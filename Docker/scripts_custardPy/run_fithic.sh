#!/bin/bash -e
cmdname=`basename $0`
function usage()
{
    echo "Usage: ${0##*/} <pairfile> <odir> <label> <resolution>"
    echo '   <pairfile>: .pair file' 1>&2
    echo "   <odir> : Output directory"
    echo "   <label> : Label of the sample"
    echo "   <resolution> : Resolution for fithic (e.g., 25000, 50000, 100000)"
    echo '   Options:' 1>&2
    echo "      -g genometable: genome table file (describing the chromosome length)"
    echo '   Example:' 1>&2
    echo "      $cmdname dedup.bwa.q30.pairs.gz fithicresult/ Control 25000" 1>&2
}

while getopts g: option; do
    case ${option} in
    g) gt=${OPTARG} ;;
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

pairfile=$1
odir=$2
prefix=$3
resolution=$4

if test "$gt" = ""; then
    echo "Error: specify genome table (-g)."
    exit 0
fi

ex(){ echo $1; eval $1; }

echo -e "\nRun fithic..."
mkdir -p $odir

ex "zcat $pairfile | grep -v \# | awk -F'\t' 'BEGIN{OFS=\"\t\"} {print \$1, \$2, \$3, \$6, \$4, \$5, \$7, \$5 - \$3}' \
    > $odir/dedup.bwa.q30.validpair"

ex "/opt/fithic/fithic/utils/validPairs2FitHiC-fixedSize.sh \
     $resolution \
     ${prefix}_$resolution \
     $odir/dedup.bwa.q30.validpair \
     $odir"

ex "rm $odir/dedup.bwa.q30.validpair"

ex "/opt/fithic/fithic/utils/createFitHiCFragments-fixedsize.py \
       --chrLens $gt \
       --resolution $resolution \
       --outFile $odir/fragmentsfile.$resolution.gz"

ex "/opt/fithic/fithic/utils/HiCKRy.py \
       -i $odir/${prefix}_${resolution}_fithic.contactCounts.gz \
       -f $odir/fragmentsfile.$resolution.gz \
       -o $odir/biasValues.$resolution.gz"

fitoutput=$odir

ex "/opt/fithic/fithic/fithic.py \
       -i $odir/${prefix}_${resolution}_fithic.contactCounts.gz \
       -f $odir/fragmentsfile.$resolution.gz \
       -t $odir/biasValues.$resolution.gz \
       -r $resolution \
       -o $fitoutput \
       -l ${prefix}_${resolution} \
       -U 1000000 \
       -L 15000 \
       -v"

ex "/opt/fithic/fithic/utils/createFitHiCHTMLout.sh \
     ${prefix}_${resolution} \
     1 \
     $fitoutput"

ex "/opt/fithic/fithic/utils/merge-filter.sh \
     $fitoutput/${prefix}_${resolution}.spline_pass1.res$resolution.significances.txt.gz \
     5000 \
     $fitoutput/${prefix}_${resolution}.merged.gz \
     0.01"

ex "rm $odir/${prefix}_${resolution}_fithic.contactCounts.gz $odir/biasValues.${resolution}.gz $odir/fragmentsfile.${resolution}.gz"