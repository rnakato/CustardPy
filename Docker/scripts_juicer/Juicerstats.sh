#!/bin/bash

cmdname=`basename $0`
function usage()
{
    echo "$cmdname <CustardPy dir> <normalization type>" 1>&2
}

if [ $# -ne 2 ]; then
  usage
  exit 1
fi


odir=$1
norm=$2
statsfile=inter_30.txt

outputstats(){
    odir=$1
    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
	echo -e "Sample\t`parseJuicerstats.pl -h $odir/$cell/aligned/$statsfile`\tNumber of TADs\tNumber of loops"
        break
    done

    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
        echo -en "`basename $cell`\t"
	echo -en "`parseJuicerstats.pl $odir/$cell/aligned/$statsfile`\t"
	echo -en "`cat $odir/$cell/TAD/$norm/25000_blocks.bed | wc -l`\t"
	grep -v \# $odir/$cell/loops/$norm/merged_loops.bedpe | wc -l
    done
}

outputstats $odir > $odir/Juicerstats.tsv
#csv2xlsx.pl -i $odir/Juicerstats.tsv -o $odir/Juicerstats.xlsx
echo "Juicer stats are written to $odir/Juicerstats.tsv."
