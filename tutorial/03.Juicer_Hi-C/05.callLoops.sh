build=hg38

sing="apptainer exec --nv --bind /work,/work2,/work3 /work3/SingularityImages/custardpy.3.5.2.sif"

cell=Control # siCTCF siRad21
odir=CustardPyResults/Juicer_$build/$cell
hic=$odir/aligned/inter_30.hic
### In case of starting from .hic files:
#hic=hic/$cell/GSE196034_${cell}_merged.hic

for norm in VC VC_SQRT KR SCALE
do
    $sing call_HiCCUPS.sh $norm $odir $hic
done
#motifdir=
#$sing call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe

# Mustache
hic=CustardPyResults/Cooler_hg38/Control/hic/contact.bwa.q30.hic
$sing mustache -f $hic -norm SCALE -p 12 -r 5kb -pt 0.05 -o mustache.hic.5kb.tsv
