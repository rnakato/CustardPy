#!/bin/bash

ex(){ echo $1; eval $1; }

odir=$1

sam="$odir/*/*.sam" 
gzsort="$odir/splits/*.fastq.gz.sort.txt"
nodups="$odir/aligned/merged_nodups.txt"
merged_sort="$odir/aligned/merged_sort.txt"
dups="$odir/aligned/dups.txt"
opt="$odir/aligned/opt_dups.txt"
collisions="$odir/aligned/collisions*.txt"

echo "The following files will be gzipped:"
echo "   $sam"
echo "   $gzsort"
echo "   $nodups"
echo "   $merged_sort"
echo "   $dups"
echo "   $opt"
echo "   $collisions"

ex "pigz $sam $gzsort $nodups $merged_sort $dups $opt $collisions"

#ex "pigz $odir/*/*.sam $odir/splits/*.fastq.gz.sort.txt $odir/aligned/merged_nodups.txt $odir/aligned/merged_sort.txt $odir/aligned/dups.txt $odir/aligned/opt_dups.txt $odir/aligned/collisions*.txt"
