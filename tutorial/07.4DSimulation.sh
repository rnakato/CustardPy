#!/bin/bash

singjuicer="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"
sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

chr=chr21
start=24000000
end=32000000
resolution=100000
norm=SCALE
odir=JuicerResults/Hap1-A
hic=$odir/aligned/inter_30.hic

cell=Hap1-A
dir=$odir/4dmodel/$resolution/$chr/$start-$end/
mkdir -p $dir
#$singjuicer juicertools.sh dump observed $norm $hic $chr:$start:$end $chr:$start:$end BP $resolution $dir/$cell.txt

hic=$odir/aligned/inter_30.hic
$sing strawMatrixFromHiC.py $hic $dir/$cell.txt $chr --start $start --end $end
$sing phic preprocessing --input $dir/$cell.txt --res $resolution --plt-max-c 0.1

exit

func(){
    cell=$1
    chr=$2
    start=$3
    end=$4

    pdir="/work/git/PHi-C/Tutorial/"

    # Parameters for plot
    PLT_MAX_K_BACKBONE=0.5
    PLT_MAX_K=0.010
    PLT_K_DIS_BINS=100
    PLT_MAX_K_DIS=1500
    OFFSET=2
    PLT_MAX_LOG_C=0
    PLT_MIN_LOG_C=-3

    # Hyperparameters for optimization
    ALPHA1=0.002
    ALPHA2=0.0001
    STEP1=2000
    STEP2=5000
    ITERATION=100
    INIT_K_BACKBONE=0.3

    # Number of optimization sample
    SAMPLE=1
    cores=64

    $sing juicertools.sh dump observed $norm $hic $chr:$start:$end $chr:$start:$end BP $res $cell.txt
    $pdir/1_conversion.py $cell.txt $res
    $pdir/2_normalization.py $cell $res $OFFSET $PLT_MAX_LOG_C $PLT_MIN_LOG_C
    $pdir/3_optimization.py --gpu $cell -p $cores
    python $pdir/4_validation.py $cell $res $SAMPLE $PLT_MAX_LOG_C $PLT_MIN_LOG_C $PLT_MAX_K_BACKBONE $PLT_MAX_K $PLT_K_DIS_BINS $PLT_MAX_K_DIS
    # Name
    KFILE=$cell/optimized_data/00_K.txt

    # Parameters
    FRAME=1000
    SAMPLE=100

#---------------------------------------------------------------------------------------------------
    # Run python codes
    python $pdir/5_4d_simulation.py $KFILE $FRAME $cell
    python $pdir/6_conformation.py $KFILE $SAMPLE $cell
}

cell=Hap1-A
odir=$(pwd)/JuicerResults/$cell

hic=$odir/aligned/inter_30.hic
norm=SCALE
resolution=100000

s=16000000
e=24000000
func $cell chr21 $s $e


exit
N=`echo $s $e $res | awk '{printf ("%d",($2-$1)/$3 +1)}'`

comp=/work2/Hi-C/Sakata_RPE/JuicerResults_merged_20201031/$cell/mega/Eigen/$res/Compartment.$norm.chr$chr.All.bed
../../convert_compartment.pl $cell/polymer_N$N.psf $comp $res 16000000 24000000 > $cell/polymer_N$N.compartment.psf

FRAME=1000
python ../plot_distance_map.py $cell/dynamics_00_K.xyz $N $FRAME $cell/distance_map/

#    ffmpeg -framerate 24 -i $cell/distance_map/frame_%04d.svg -pix_fmt yuv420p \
    #           -c:v libx264 -preset slow -crf 24 -r 24 \
    #           -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
    #           $cell/distance_map/movie.mp4
