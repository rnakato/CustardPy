#!/bin/bash

build=hg38
fastq_post="_"  # "_" or "_R"  before .fastq.gz
enzyme=MboI

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
#gt=/work/Database/UCSC/$build/genome_table
sing="singularity exec custardpy.sif"
gt=genometable.$build.txt

odir=CustardPyResults_Hi-C/Juicer_$build/Hap1-A
#odir=CustardPyResults_Hi-C/Juicer_$build/WaplKO_3.3-A/

hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=25000

# Contact matrix
echo "generate Matrix..."
$sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt

# InsulationScore
echo "calculate Insulation score.."
$sing makeInslationScore.sh $norm $odir $resolution $gt
