#!/bin/bash

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

build=hg38
odir=CustardPyResults/Cooler_$build/

$sing Coolerstats.sh $odir
