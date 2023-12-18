#!/bin/bash

for build in hg19 hg38 mm10 mm39 rn7 galGal6 ce11 danRer11 dm6 xenLae2 sacCer3
do
    genome=/work/Database/Database_fromDocker/Referencedata_${build}/genome.fa
    gt=/work/Database/Database_fromDocker/Referencedata_${build}/genometable.txt
    for enzyme in AluI #HindIII DpnII MboI Sau3AI #Arima not supported
    do
        cooler digest $gt $genome $enzyme > ${enzyme}_resfrag_$build.bed
        pigz ${enzyme}_resfrag_$build.bed
    done
done

for build in S.pombe HVAEP T2T
do
    genome=/work/Database/Database_fromDocker/${build}/genome.fa
    gt=/work/Database/Database_fromDocker/${build}/genometable.txt
    for enzyme in AluI #HindIII DpnII MboI Sau3AI #Arima not supported
    do
        cooler digest $gt $genome $enzyme > ${enzyme}_resfrag_$build.bed
        pigz ${enzyme}_resfrag_$build.bed
    done
done
