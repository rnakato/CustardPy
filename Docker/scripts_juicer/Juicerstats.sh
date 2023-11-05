#!/bin/bash

odir=$1
statsfile=inter_30.txt

for cell in $odir/*
do
    echo "Sample\t`parseJuicerstats.pl -h $cell/aligned/$statsfile`"
    break
done

for cell in $odir/*
do
    echo -n "`basename $cell`\t"
    parseJuicerstats.pl $cell/aligned/$statsfile
done
