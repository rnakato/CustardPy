build=hg38

#sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy.1.4.3.sif"
sing="singularity exec --nv custardpy.sif"

odir=CustardPyResults_Hi-C/Juicer_$build/Hap1-A
#odir=CustardPyResults_Hi-C/Juicer_$build/WaplKO_3.3-A
hic=$odir/aligned/inter_30.hic
norm=SCALE
$sing call_HiCCUPS.sh $norm $odir $hic
#$sing call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe
