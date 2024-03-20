#!/bin/bash

ex(){ echo $1; eval $1; }

odir=$1
ex "pigz $odir/*/*.sam $odir/splits/*.fastq.gz.sort.txt $odir/aligned/merged_nodups.txt $odir/aligned/merged_sort.txt $odir/aligned/dups.txt $odir/aligned/opt_dups.txt $odir/aligned/collisions*.txt"
