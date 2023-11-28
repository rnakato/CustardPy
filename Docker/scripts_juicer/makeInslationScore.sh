#!/bin/bash
cmdname=`basename $0`
function usage()
{
    echo "$cmdname <norm> <odir> <resolution> <gt>" 1>&2
    echo '   <norm>: normalization type (NONE|VC|VC_SQRT|KR|SCALE)' 1>&2
    echo '   <odir>: output directory (e.g., "JuicerResults/sample1")' 1>&2
    echo '   <resolution>: resolution of the matrix' 1>&2
    echo '   <gt>: genome table' 1>&2
}

if [ $# -ne 4 ]; then
  usage
  exit 1
fi

norm=$1
matrixdir=$2
binsize=$3
gt=$4

pwd=$(cd $(dirname $0) && pwd)
chrlist=$(getchr_from_genometable.sh $gt)

dir=$matrixdir/InsulationScore/$norm/$binsize
mkdir -p $dir

ex(){ echo $1; eval $1; }

for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    i=$(echo $chr | sed -e 's/chr//g')
    echo "chr$i"
    matrix=$matrixdir/Matrix/intrachromosomal/$binsize/observed.$norm.chr$i.matrix.gz

    if test -s $matrix; then
	    ex "InsulationScore.py            $matrix $dir/Insulationscore.chr$i.100k chr$i $binsize --distance 100000"
	    ex "InsulationScore.py            $matrix $dir/Insulationscore.chr$i.500k chr$i $binsize --distance 500000"
	    ex "InsulationScore.py            $matrix $dir/Insulationscore.chr$i.1M   chr$i $binsize --distance 1000000"
	    ex "plotInsulationScore           $matrix $dir/Insulationscore.chr$i.IS.png $binsize"
	    ex "plotMultiScaleInsulationScore $matrix $dir/Insulationscore.chr$i.MI.png $binsize"

        # obtain insulated boundaries and TADs
        thre=0.7
        for d in 100k 500k 1M; do
            ex "getBoundaryfromInsulationScore --thre $thre \
                $dir/Insulationscore.chr$i.$d.$binsize.bedGraph \
                $dir/Boundary.chr$i.$d.$binsize.thre${thre}.bed"
            ex "bedtools complement -L -i $dir/Boundary.chr$i.$d.$binsize.thre${thre}.bed -g $gt \
                > $dir/TAD.chr$i.$d.$binsize.thre${thre}.bed"
        done
    else
	    echo "Warning: $matrix does not exist."
    fi
done

for d in 100k 500k 1M; do
    cat $dir/Insulationscore.chr*.$d.$binsize.bedGraph > $dir/Insulationscore.genome.$d.$binsize.bedGraph
    cat $dir/Boundary.chr*.$d.$binsize.thre${thre}.bed > $dir/Boundary.genome.$d.$binsize.thre${thre}.bed
    cat $dir/TAD.chr*.$d.$binsize.thre${thre}.bed > $dir/TAD.genome.$d.$binsize.thre${thre}.bed
done
