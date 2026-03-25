sing="apptainer exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.1.1.sif"

odir=CustardPyResults/Juicer_$build/
norm=SCALE
$sing Juicerstats.sh $odir $norm
