#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.1.1.sif"

cell=Control
cool=CustardPyResults/Cooler_$build/$cell/cool/cooler.dedup.q30.cool
boundary=CustardPyResults/Cooler_$build/$cell/InsulationScore/25000/Boundaries.w500000.strong.bed
odir=CustardPyResults/Cooler_$build/$cell/coolpup
mkdir -p $odir

# TAD boundary average plot
$sing coolpup.py $cool $boundary --local --n_proc 64 --outname $odir/$cell.boundary.clpy
$sing plotpup.py --input_pups $odir/$cell.boundary.clpy --rownames $cell --output $odir/$cell.boundary.png

# TAD boundary normalized by random shifts
$sing coolpup.py $cool $boundary --local --n_proc 64 --nshifts 1  --ignore_diags 0 --flank 300000 --outname $odir/${cell}_pileup_1shift.clpy
$sing plotpup.py --input_pups $odir/${cell}_pileup_1shift.clpy --rownames $cell --output $odir/${cell}_pileup_1shift.png

# Loop average plot (APA plot)
loop=CustardPyResults/Cooler_$build/$cell/loops/cooltools/dots.10000.tsv
$sing coolpup.py $cool $loop --n_proc 64 --outname $odir/APA.$cell.clpy  --mindist 10000
$sing plotpup.py --input_pups $odir/APA.$cell.clpy --rownames $cell --output $odir/APA.${cell}.png --vmin 0.5 --vmax 12.0

# See more examples in https://coolpuppy.readthedocs.io/en/latest/Examples/Walkthrough_CLI.html
