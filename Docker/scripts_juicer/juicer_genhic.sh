#!/bin/bash -e
cmdname=`basename $0`
function usage()
{
    echo "$cmdname [-e enzyme] [-m tmpdir] [-p ncore] <fastqdir> <odir> <build> <gt> <bwaindex> <enzyme> <fastq_post>" 1>&2
    echo '   <fastqdir>: directory that contains input fastq files (e.g., "fastq/sample1")' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <build>: genome build (e.g., hg38)' 1>&2
    echo '   <gt>: genome table' 1>&2
    echo '   <bwaindex>: index file of BWA' 1>&2
    echo '   <enzyme>: enzyme (e.g., HindIII, MboI)' 1>&2
    echo '   <fastq_post [_|_R]>: if the filename of fastqs is *_[1|2].fastq, supply "_". if *_[R1|R2].fastq, choose "_R".' 1>&2
    echo '   Options:' 1>&2
    echo '      -p ncore: number of CPUs (default: 32)' 1>&2
    echo '      -m tmpdir: tempdir' 1>&2
    echo '      -L: Allocate larger memory ("-Xms1024m -Xmx655360m", default: "-Xms512m -Xmx65536m", for deep-sequenced samples; e.g., Rao 2014)' 1>&2
    echo '   Example:' 1>&2
    echo "      $cmdname $(pwd)/fastq/Hap1-A/ $(pwd)/JuicerResults/Hap1-A hg38 genometable.hg38.txt bwaindex/hg38 HindIII _R" 1>&2
}

tmpdir=""
ncore=32
memoryparam=""
while getopts p:m:L option; do
    case ${option} in
        p) ncore=${OPTARG} ;;
        m) tmpdir=${OPTARG} ;;
	    L) memoryparam="-L" ;;
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

if [ $# -ne 7 ]; then
  usage
  exit 1
fi

fqdir=$1
odir=$2
label=$(basename $odir)
build=$3
gt=$4
bwaindex=$5
enzyme=$6
fastq_post=$7

jdir=/opt/juicer

if [ -n "$tmpdir" ]; then
  param="-p $tmpdir"
fi

if [[ ${fqdir} =~ ^/.+$ ]]; then
    fqdir=$fqdir
else
    fqdir=$(pwd)/$fqdir
fi

if [[ ${odir} =~ ^/.+$ ]]; then
    odir=$odir
else
    odir=$(pwd)/$odir
fi

ex(){ echo $1; eval $1; }

pwd=`pwd`
ex "bash $jdir/CPU/juicer.sh -t $ncore -g $build -d $odir $param $memoryparam \
     -s $enzyme -a $label -p $gt \
     -z $bwaindex -D $jdir -e $fastq_post -S final \
     2>&1 | tee $odir/juicer_genhic.log"
