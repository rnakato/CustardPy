#!/bin/bash

sing="singularity exec --bind /work,/work2 /work/SingularityImages/custardpy.0.1.0.sif"

odir=$(pwd)/JuicerResults/Hap1-A

resolution=100000
str=Hap1-A
chr=chr21
s=24000000
e=32000000
s1=$(( $s / $resolution ))
e1=$(( $e / $resolution ))

Matrix=$odir/Matrix/intrachromosomal/$resolution/observed.SCALE.$chr.matrix.gz
Eigen=$ddir/Eigen/$resolution/Compartment.SCALE.$chr.All.bed

dir=$odir/pastis/${resolution}/$chr/$s-$e
mkdir -p $dir
$sing matrix2npy.py $Matrix $dir/${str}.npy \
      --start1 $s1 --end1 $e1 \
      --start2 $s1 --end2 $e1

config=$dir/config.ini
#$sing cp /opt/scripts/pastis/config.ini $dir/config.ini
echo "[all]"               > $config
echo "verbose: 1"          >> $config
echo "max_iter: 1000"      >> $config
echo "normalize: True"     >> $config
echo "output_name: ${str}" >> $config
echo "counts: ${str}.npy"  >> $config

$sing pastis-pm2 $dir
$sing addCompartment4pastis.pl $dir/PM2.${str}.pdb $Eigen $s1
$sing addChain4pastis.pl $dir/PM2.${str}.pdb.addCompartment.pdb
