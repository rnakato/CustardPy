#!/bin/bash
sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.2.0.sif"
#sing="apptainer exec custardpy.sif"

outputdir=figure
mkdir -p $outputdir

Resdir=CustardPyResults/Cooler_hg38
resolution=25000
norm=Cooler

chr=chr20
start=8000000
end=16000000

cell=Control
$sing plotHiCMatrix \
      $Resdir/$cell/Matrix/$resolution/balanced.$chr.matrix.gz \
      $outputdir/ContactMap.$cell.$chr.$start-$end.pdf \
      $start $end $cell

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --type $norm -d 5000000 \
      -o $outputdir/IS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control/cool/cooler.dedup.q30.multires.cool:Control \
      $Resdir/siCTCF/cool/cooler.dedup.q30.multires.cool:siCTCF \
      $Resdir/siRad21/cool/cooler.dedup.q30.multires.cool:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --type $norm -d 5000000 \
      -o $outputdir/IS2.$chr.$start-$end


$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --multi --type $norm -d 5000000 \
      -o $outputdir/MultiIS.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --multidiff --type $norm -d 5000000 \
      -o $outputdir/MultiISdiff.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --compartment --type $norm -d 5000000 \
      -o $outputdir/Compartment.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --di --type $norm -d 5000000 \
      -o $outputdir/DI.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -c $chr --start $start --end $end -r $resolution \
      --drf --type $norm -d 5000000 \
      -o $outputdir/DRF.$chr.$start-$end

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --triangle_ratio_multi --type $norm -d 5000000

# If you want to get the tsv file of the logfoldchange matrices
$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/TriangleRatioMulti.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --triangle_ratio_multi --type $norm -d 5000000 \
      --output_logfc_matrix

$sing plotHiCfeature \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/virtual4C.$chr \
      -c $chr --start $start --end $end -r $resolution \
      --v4c --anchor 10400000 --vmax 100 --type $norm

$sing drawSquarePair \
      $Resdir/Control/Matrix/$resolution/balanced.$chr.matrix.gz:Control \
      $Resdir/siRad21/Matrix/$resolution/balanced.$chr.matrix.gz:siRad21 \
      -o $outputdir/SquarePair.$chr --start $start --end $end -r $resolution

$sing drawSquareRatioPair \
      $Resdir/Control/Matrix/$resolution/balanced.$chr.matrix.gz:Control \
      $Resdir/siRad21/Matrix/$resolution/balanced.$chr.matrix.gz:siRad21 \
      $Resdir/Control/Matrix/$resolution/balanced.$chr.matrix.gz:Control \
      $Resdir/siCTCF/Matrix/$resolution/balanced.$chr.matrix.gz:siCTCF \
      -o $outputdir/SquareRatioPair.$chr --start $start --end $end -r $resolution

$sing drawSquareMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/SquareMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution \
      --cooler

$sing drawSquareMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/SquareMulti_logscale.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution --log \
      --cooler

$sing drawSquareRatioMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/SquareRatioMulti.$chr \
      -c $chr --start $start --end $end --type $norm -r $resolution \
      --cooler

$sing drawTriangleMulti \
      $Resdir/Control:Control \
      $Resdir/siCTCF:siCTCF \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/TriangleMulti.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution \
      --cooler

$sing drawTrianglePair \
      $Resdir/Control:Control \
      $Resdir/siRad21:siRad21 \
      -o $outputdir/TrianglePair.$chr \
      -c $chr --start $start --end $end --type $norm -d 5000000 -r $resolution \
      --cooler
