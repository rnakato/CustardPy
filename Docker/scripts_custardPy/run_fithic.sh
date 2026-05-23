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
    echo "      -q <float>: q-value threshold (default: 0.01)"
    echo '   Example:' 1>&2
    echo "      $cmdname dedup.bwa.q30.pairs.gz fithicresult/ Control 25000" 1>&2
}

qval=0.01
while getopts g:q: option; do
    case ${option} in
    g) gt=${OPTARG} ;;
    q) qval=${OPTARG} ;;
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

cool=$1
odir=$2
prefix=$3
resolution=$4

if test "$gt" = ""; then
    echo "Error: specify genome table (-g)."
    exit 0
fi

ex(){ echo $1; eval $1; }

fithicdir=/opt/fithic/fithic
#fithicdir=./fithic/fithic

echo -e "\nRun fithic..."
mkdir -p $odir $odir/intermediate

ex "cool2fithicdata.py \
  --cool $cool::/resolutions/$resolution \
  --resolution $resolution \
  --lower $resolution \
  --upper 1000000 \
  --out-prefix $odir/intermediate/${prefix}_$resolution.fromcool"

ex "$fithicdir/fithic.py \
  -i $odir/intermediate/${prefix}_${resolution}.fromcool.contactCounts.gz \
  -f $odir/intermediate/${prefix}_${resolution}.fromcool.fragments.gz \
  -t $odir/intermediate/${prefix}_${resolution}.fromcool.biasValues.gz \
  -r $resolution \
  -o $odir \
  -l ${prefix}_${resolution} \
  -L $resolution \
  -U 1000000 \
  -v"

ex "$fithicdir/utils/createFitHiCHTMLout.sh \
     ${prefix}_${resolution} \
     1 \
     $odir"

ex "$fithicdir/utils/merge-filter.sh \
     $odir/${prefix}_${resolution}.spline_pass1.res$resolution.significances.txt.gz \
     $resolution \
     $odir/${prefix}_${resolution}.merged.gz \
     $qval > $odir/${prefix}_${resolution}.merged.log"

echo -e "\nNumber of candidate interactions:"
zcat $odir/${prefix}_${resolution}.spline_pass1.res$resolution.significances.txt.gz |wc -l 
echo "Number of significant interactions after merging (FDR < $qval):"
zcat $odir/${prefix}_${resolution}.merged.gz |wc -l

