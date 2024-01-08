#sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.5.0.sif"
sing="singularity exec custardpy.sif"

odir=CustardPyResults_Hi-C/Juicer_$build/
norm=SCALE
$sing Juicerstats.sh $odir $norm
