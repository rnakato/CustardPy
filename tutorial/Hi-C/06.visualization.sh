#!/bin/bash

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.2.0.sif"

chr=chr20
start=8000000
end=16000000
norm=SCALE
cell=Hap1-A
matrix=JuicerResults_hg38/$cell/Matrix/intrachromosomal/100000/observed.$norm.$chr.matrix.gz

$sing plotHiCMatrix $matrix ContactMap.$cell.$chr.$start-$end.png $start $end $cell

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --compartment --type $norm -d 5000000 \
      -o Compartment.$chr.$start-$end

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --multi --type $norm -d 5000000 \
      -o MultiIS.$chr.$start-$end

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --multidiff --type $norm -d 5000000 \
      -o MultiISdiff.$chr.$start-$end

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --dfr --type $norm -d 5000000 \
      -o DFR.$chr.$start-$end

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --di --type $norm -d 5000000 \
      -o DI.$chr.$start-$end

$sing drawSquarePair \
         Control/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         Rad21KD_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         drawSquarePair.chr21 --start 24000000 --end 32000000

$sing drawSquareMulti Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL SquareMulti.chr9 \
          chr9 --start 1000000 --end 38000000  --type VC_SQRT --vmax 20
