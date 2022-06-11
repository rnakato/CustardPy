Quickstart
=====================

Full command example
----------------------------------------------------------------

These scripts assume that the fastq files are stored in ``fastq/$cell`` (e.g., ``fastq/Control_1``).
The outputs are stored in ``JuicerResults/$cell``.

The whole commands using the Singularity image (``rnakato_juicer.sif``) are as follows:

.. code-block:: bash

    build=hg38
    fastq_post="_R"  # "_" or "_R"  before .fastq.gz
    enzyme=MboI      # enzyme type
    norm=SCALE       # normalization type

    gt=genome_table.$build.txt  # genome_table file
    bwaindex=bwa-indexes/UCSC-$build  # BWA index file
    gene=refFlat.$build.txt # gene annotation (refFlat format)
    ncore=64 # number of CPUs

    sing="singularity exec --nv --bind /work custardpy_juicer.sif" # singularity command

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
        $sing call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe

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
