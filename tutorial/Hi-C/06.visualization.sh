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
      --type $norm -d 5000000 \
      -o IS.$chr.$start-$end

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
      --drf --type $norm -d 5000000 \
      -o DRF.$chr.$start-$end

$sing plotHiCfeature \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --di --type $norm -d 5000000 \
      -o DI.$chr.$start-$end

$sing drawSquarePair \
      JuicerResults_hg38/Hap1-A/Matrix/intrachromosomal/25000/observed.$norm.$chr.matrix.gz \
      JuicerResults_hg38/WaplKO_3.3-A/Matrix/intrachromosomal/25000/observed.$norm.$chr.matrix.gz \
      SquarePair.$chr --start $start --end $end -r 25000

$sing drawSquareMulti \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -o SquareMulti.$chr \
      -c $chr --start $start --end $end --type $norm

$sing drawSquareRatioMulti \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -o SquareRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm

$sing drawTriangleMulti \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -o TriangleMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000

$sing drawTrianglePair \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -o TrianglePair.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000

$sing drawTriangleRatioMulti \
      JuicerResults_hg38/Hap1-A:Control \
      JuicerResults_hg38/WaplKO_3.3-A:WaplKO \
      -o TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000