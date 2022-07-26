#!/bin/bash

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
