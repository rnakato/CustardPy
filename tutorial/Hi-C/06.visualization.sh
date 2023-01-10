#!/bin/bash

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.2.0.sif"
sing="singularity exec custardpy.sif"

chr=chr20
start=8000000
end=16000000
norm=SCALE
cell=Hap1-A
resolution=25000
Resdir=JuicerResults_hg38
matrix=$Resdir/$cell/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz

$sing plotHiCMatrix $matrix ContactMap.$cell.$chr.$start-$end.png $start $end $cell


$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --type $norm -d 5000000 \
      -o IS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --compartment --type $norm -d 5000000 \
      -o Compartment.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --multi --type $norm -d 5000000 \
      -o MultiIS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --multidiff --type $norm -d 5000000 \
      -o MultiISdiff.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --drf --type $norm -d 5000000 \
      -o DRF.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -c $chr --start $start --end $end \
      --di --type $norm -d 5000000 \
      -o DI.$chr.$start-$end

$sing drawSquarePair \
      $Resdir/Hap1-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz \
      $Resdir/WaplKO_3.3-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz \
      SquarePair.$chr --start $start --end $end -r $resolution

$sing drawSquareMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o SquareMulti.$chr \
      -c $chr --start $start --end $end --type $norm

$sing drawSquareRatioMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o SquareRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm

$sing drawTriangleMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o TriangleMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000

$sing drawTrianglePair \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o TrianglePair.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000

$sing drawTriangleRatioMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000
