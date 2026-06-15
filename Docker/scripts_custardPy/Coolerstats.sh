#!/bin/bash

cmdname=`basename $0`
function usage()
{
    echo "$cmdname <CustardPy dir>" 1>&2
    echo "  Example: $cmdname  CustardPyResults/Cooler_hg38/" 1>&2
}

if [ $# -ne 1 ]; then
  usage
  exit 1
fi

odir=$1

outputstats(){
    odir=$1
    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
	echo -en "Sample\t"
	echo -en "`parseCoolerStats.py --header $odir/$cell/qc_report/mapping_stats.txt`\t%\t"
	echo -e "Number of TADs (100kbp)\tNumber of loops (dots 10kbp)\tNumber of loops (fithic 25kbp)\tNumber of loops (mustache 5kb)"
	break
    done

    for cell in `ls -l $odir/ | grep ^d | awk '{print $9}'`
    do
	echo -en "$cell\t"
	echo -en "`parseCoolerStats.py $odir/$cell/qc_report/mapping_stats.txt `\t"
	echo -en "`cat $odir/$cell/InsulationScore/25000/TADs.w100000.bed | wc -l`\t"
	echo -en "`cat $odir/$cell/loops/cooltools/dots.10000.tsv | wc -l`\t"
	echo -en "`zcat $odir/$cell/loops/fithic/25000/${cell}_25000.merged.gz | wc -l`\t"
	echo -en "`cat $odir/$cell/loops/mustache/loop.cool.5kb.tsv | wc -l`\n"
    done
}

outputstats $odir > $odir/Coolerstats.tsv
echo "Cooler stats are written to $odir/Coolerstats.tsv."
