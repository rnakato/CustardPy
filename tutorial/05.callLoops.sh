build=hg38
gt=/work/Database/UCSC/$build/genome_table
gene=/work/Database/UCSC/$build/refFlat.txt
odir=$(pwd)/JuicerResults/Hap1-A

hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000

singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif call_HiCCUPS.sh $norm $odir $hic $build
