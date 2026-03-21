#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="singularity exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.0.0.sif"
#sing="singularity exec custardpy.sif"

cell=Control
cool=CustardPyResults_Hi-C/Cooler_$build/$cell/cool/cooler.dedup.q30.multires.cool
odir=CustardPyResults_Hi-C/Cooler_$build/$cell/InsulationScore

resolutions="25000 50000 100000"

for resolution in $resolutions
do
    echo "calculate Insulation score.."
#    $sing cooltools_insulation.py $cool $odir $gt --resolution $resolution
    $sing /work3/DockerFiles/CustardPy/Docker/scripts_custardPy/cooltools_insulation.py $cool $odir $gt --resolution $resolution
done
