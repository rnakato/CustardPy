build=hg38

#sing="singularity exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.1.5.0.sif"
sing="singularity exec --nv custardpy.sif"

cell=Control # siCTCF siRad21 siNIPBL
odir=CustardPyResults_Hi-C/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

norm=SCALE
$sing call_HiCCUPS.sh $norm $odir $hic

motifdir=
#$sing call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe
