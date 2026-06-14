#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

odir=CustardPyResults/Cooler_$build/

$sing Coolerstats.sh $odir