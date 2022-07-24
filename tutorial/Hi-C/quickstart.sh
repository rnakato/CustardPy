#!/bin/bash

build=hg38
gt=/work/Database/UCSC/$build/genome_table
bwaindex=/work/Database/bwa-indexes/UCSC-$build
gene=/work/Database/UCSC/$build/refFlat.txt
tmpdir=/tmp/juicer.1.6.1.`date +%Y%m%d%H%M`
fastq_post="_"  # "_" or "_R"
enzyme=MboI
norm=SCALE
ncore=64 # number of CPUs

sing="singularity exec --nv --bind /work,/work2 /work/SingularityImages/custardpy_juicer.0.1.0.sif"

for cell in `ls fastq/* -d`  # for all directories in fastq/
do
    cell=$(basename $cell)
    fqdir=$(pwd)/fastq/$cell/
    odir=$(pwd)/JuicerResults/$cell
    echo $cell

    # generate .hic file from fastq by Juicer
    rm -rf $odir
    $sing juicer_map.sh -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

    # Compress intermediate files
    $sing juicer_pigz.sh $odir

    # plot contact frequency
    if test ! -e $odir/distance; then $sing plot_distance_count.sh $cell $odir; fi

    hic=$odir/aligned/inter_30.hic
    # call TADs (arrowHead)
    $sing juicer_callTAD.sh $norm $hic $odir $gt

    # call loops (HICCUPS, add '--nv' option to use GPU)
    $sing call_HiCCUPS.sh $norm $odir $hic $build
    # motif analysis
    $sing juicertools.sh motifs $build $motifdir $odir/loops/$norm/merged_loops.bedpe hg38.motifs.txt

    for resolution in 25000 50000 100000
    do
        # make contact matrix for all chromosomes
        $sing makeMatrix_intra.sh $norm $odir $hic $resolution $gt
        # calculate Pearson coefficient and Eigenvector
        $sing makeEigen.sh Pearson $norm $odir $hic $resolution $gt $gene
        $sing makeEigen.sh Eigen $norm $odir $hic $resolution $gt $gene
        # calculate insulation score
        $sing makeInslationScore.sh $norm $odir $resolution $gt
    done
done
