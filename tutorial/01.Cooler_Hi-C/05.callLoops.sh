#!/bin/bash

build=hg38
gt=genometable.$build.txt

sing="apptainer exec --bind /work,/work2,/work3 /work/SingularityImages/custardpy.3.5.2.sif"

cell=Control
odir=CustardPyResults/Cooler_$build/$cell/
cool=$odir/cool/contact.bwa.q30.mcool

mkdir -p $odir/log

# cooltools dots (HICCUPS-like)
for resolution in 5000 10000 25000
do
    echo "call loops by cooltools.."
    $sing cooltools_dots.py $cool $odir $gt --resolution $resolution -p $ncore
done

# Fit-Hi-C
for resolution in 10000 25000
do
    echo "Run Fit-Hi-C.."
    odir_fithic=$odir/loops/fithic/$resolution
    $sing run_fithic.sh -g $gt $cool $odir_fithic $cell $resolution | tee $odir/log/fithic.$resolution.log
done

# Mustache
$sing mustache -f $cool \
	    -r 5kb -norm weight \
	    -pt 0.05 -p 12 \
	    -o mustache.cool.5kb.tsv
