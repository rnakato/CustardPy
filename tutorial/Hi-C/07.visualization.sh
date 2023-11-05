#!/bin/bash

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec custardpy.sif"

chr=chr20
start=8000000
end=16000000
norm=SCALE
cell=Hap1-A
resolution=25000
Resdir=ustardPyResults_Hi-C/Juicer_hg38
matrix=$Resdir/$cell/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz


$sing plotHiCMatrix \
      $matrix \
      $Resdir/ContactMap.$cell.$chr.$start-$end.pdf \
      $start $end $cell

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --type $norm -d 5000000 \
      -o $Resdir/IS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --multi --type $norm -d 5000000 \
      -o $Resdir/MultiIS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --multidiff --type $norm -d 5000000 \
      -o $Resdir/MultiISdiff.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --compartment --type $norm -d 5000000 \
      -o $Resdir/Compartment.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --di --type $norm -d 5000000 \
      -o $Resdir/DI.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -c $chr --start $start --end $end -r $resolution \
      --drf --type $norm -d 5000000 \
      -o $Resdir/DRF.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -o $Resdir/TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --triangle_ratio_multi --type $norm -d 5000000

$sing plotHiCfeature \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -o $Resdir/virtual4C.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --v4c --anchor 10400000 --vmax 100 --type $norm

$sing drawSquarePair \
      $Resdir/Hap1-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/WaplKO_3.3-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:WaplKO \
      -o $Resdir/SquarePair.$chr --start $start --end $end -r $resolution

$sing drawSquareRatioPair \
      $Resdir/Hap1-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/WaplKO_3.3-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:WaplKO \
      $Resdir/Hap1-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/SCC4KO-A/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:SCC4KO \
      -o $Resdir/drawSquareRatioPair.$chr --start $start --end $end -r $resolution

$sing drawSquareMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -o $Resdir/SquareMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution

$sing drawSquareRatioMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -o $Resdir/SquareRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution

$sing drawTriangleMulti \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      $Resdir/SCC4KO-A:SCC4KO \
      -o $Resdir/TriangleMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution

$sing drawTrianglePair \
      $Resdir/Hap1-A:Control \
      $Resdir/WaplKO_3.3-A:WaplKO \
      -o $Resdir/TrianglePair.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution
