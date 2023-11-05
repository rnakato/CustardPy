build=hg38

#sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec custardpy.sif"

odir=CustardPyResults_Hi-C/Juicer_$build/
norm=SCALE
$sing Juicerstats.sh $odir $norm
