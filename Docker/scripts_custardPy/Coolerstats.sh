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

outputstats(){
    odir=$1
    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
	echo -en "Sample\t"
	echo -en "`parseCoolerStats.py --header $odir/$cell/qc_report/mapping_stats.txt`\t%\t"
	echo -e "Number of TADs\tNumber of loops"
	break
    done

    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
	echo -en "$cell\t"
	echo -en "`parseCoolerStats.py $odir/$cell/qc_report/mapping_stats.txt `\t"
	echo -en "`cat $odir/$cell/TAD/$norm/25000_blocks.bed | wc -l`\t"
	grep -v \# $odir/$cell/loops/$norm/merged_loops.bedpe | wc -l
    done
}

outputstats $odir > $odir/Coolerstats.tsv
echo "Cooler stats are written to $odir/Coolerstats.tsv."
