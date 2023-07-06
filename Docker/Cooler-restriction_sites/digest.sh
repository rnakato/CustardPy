#!/bin/bash

for build in #hg19 hg38 mm10 mm39 rn7 galGal6 ce11 danRer11 dm6 xenLae2 sacCer3
do
    genome=/work/Database/Database_fromDocker/Referencedata_${build}/genome.fa
    gt=/work/Database/Database_fromDocker/Referencedata_${build}/genometable.txt
    for enzyme in HindIII DpnII MboI Sau3AI #Arima not supported
    do
        cooler digest $gt $genome $enzyme > ${enzyme}_resfrag_$build.bed
        pigz ${enzyme}_resfrag_$build.bed
    done
done


for build in S.pombe HVAEP T2T
do
    genome=/work/Database/Database_fromDocker/${build}/genome.fa
    gt=/work/Database/Database_fromDocker/${build}/genometable.txt
    for enzyme in HindIII DpnII MboI Sau3AI #Arima not supported
    do
        cooler digest $gt $genome $enzyme > ${enzyme}_resfrag_$build.bed
        pigz ${enzyme}_resfrag_$build.bed
    done
done

exit

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

build=mm39
fa=/work/Database/UCSC/${build}/genome.fa

for re in MboI HindIII DpnII
do
    $sing python /opt/juicer/misc/generate_site_positions.py $re $build $fa
done

#re="DpnII"
#resite="^GATC"
#build="sacCer3"
#type="UCSC"
#oname="${re}_resfrag_${build}.bed"
#fa="/home/Database/new/${type}/${build}/genome.fa"
#/home/git/HiC-Pro_2.9.0/bin/utils/digest_genome.py -r $resite -o $oname $fa
