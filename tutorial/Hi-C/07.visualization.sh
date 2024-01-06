#!/bin/bash

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.6.0.sif"
#sing="singularity exec custardpy.sif"

chr=chr20
start=8000000
end=16000000
norm=SCALE
cell=Control
resolution=25000
Resdir=CustardPyResults_Hi-C/Juicer_hg38
matrix=$Resdir/$cell/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz

$sing plotHiCMatrix \
      $matrix \
      $Resdir/ContactMap.$cell.$chr.$start-$end.pdf \
      $start $end $cell

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --type $norm -d 5000000 \
      -o $Resdir/IS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --multi --type $norm -d 5000000 \
      -o $Resdir/MultiIS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --multidiff --type $norm -d 5000000 \
      -o $Resdir/MultiISdiff.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --compartment --type $norm -d 5000000 \
      -o $Resdir/Compartment.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --di --type $norm -d 5000000 \
      -o $Resdir/DI.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -c $chr --start $start --end $end -r $resolution \
      --drf --type $norm -d 5000000 \
      -o $Resdir/DRF.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -o $Resdir/TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --triangle_ratio_multi --type $norm -d 5000000

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -o $Resdir/virtual4C.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --v4c --anchor 10400000 --vmax 100 --type $norm

$sing drawSquarePair \
      $Resdir/Control/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/siRad21/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:siRad21 \
      -o $Resdir/SquarePair.$chr --start $start --end $end -r $resolution

$sing drawSquareRatioPair \
      $Resdir/Control/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/siRad21/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:siRad21 \
      $Resdir/Control/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:Control \
      $Resdir/siCTCF/Matrix/intrachromosomal/$resolution/observed.$norm.$chr.matrix.gz:siCTCF \
      -o $Resdir/drawSquareRatioPair.$chr --start $start --end $end -r $resolution

$sing drawSquareMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -o $Resdir/SquareMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution

$sing drawSquareRatioMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -o $Resdir/SquareRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution

$sing drawTriangleMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      $Resdir/siNIPBL:siNIPBL \
      -o $Resdir/TriangleMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution

$sing drawTrianglePair \
      $Resdir/Control:Control \
      $Resdir/siRad21:siRad21 \
      -o $Resdir/TrianglePair.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution
