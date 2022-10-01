#!/bin/bash

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.2.0.sif"

chr=chr7
start=15000000
end=40000000
norm=SCALE
cell=Hap1-A
binsize=100000
matrix=JuicerResults_hg38/$cell/Matrix/intrachromosomal/$binsize/observed.$norm.$chr.matrix.gz

$sing plotHiCMatrix $matrix ContactMap.$cell.$chr.$start-$end.png $start $end $cell

dir1=JuicerResults_hg38/Hap1-A
$sing plotHiCfeature $dir1:Hap1-A Compartment.chr9.1M-38M --compartment \
      chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000


$sing drawSquarePair \
         Control/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         Rad21KD_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         drawSquarePair.chr21 --start 24000000 --end 32000000

$sing drawSquareMulti Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL SquareMulti.chr9 \
          chr9 --start 1000000 --end 38000000  --type VC_SQRT --vmax 20
