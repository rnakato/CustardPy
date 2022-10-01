build=hg38
odir=JuicerResults_$build/Hap1-A

sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.2.0.sif"

hic=$odir/aligned/inter_30.hic
norm=SCALE
$sing call_HiCCUPS.sh $norm $odir $hic
#$sing call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe
