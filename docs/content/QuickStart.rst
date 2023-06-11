Quickstart
=====================

.. A common problem in Hi-C analysis is the strict requirement of specific input formats. Many tools require input data to be in a specific format, and consequently, their use is hindered if the data under investigation does not conform to these specifications.

.. Since CustardPy covers the processing of Hi-C data from FASTQ and uses the generated data for the subsequent analysis, users can avoid the potential format incompatibility.

.. note::

    As the CustardPy commands below are included in the CustardPy docker image, you need to add docker or singularity commands as shown below.

.. code-block:: bash

    # This example command will mount the /work directory of the host machine
    # For docker
    singularity exec [--nv] --bind /work custardpy.sif <command>
    # For singularity
    docker run --rm -it [--gpus all] -v /work:/work rnakato/custardpy <command>

See also the sample scripts in the `tutorial/ <https://github.com/rnakato/CustardPy/tree/main/tutorial>`_ directory.


Hi-C analysis using Juicer
---------------------------------------------

Hi-C analysis from FASTQ files
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can implement whole commands for Juicer analysis from FASTQ files using ``custardpy_juicer`` command.

.. code-block:: bash

    build=hg38   # genome build
    gt=genometable.$build.txt # genome_table file
    gene=refFlat.$build.txt   # gene annotation (refFlat format)
    bwaindex=bwa-indexes/$build  # BWA index file
    ncore=64  # number of CPUs

    cell=Hap1-A
    fastq_post="_"  # "_" or "_R"
    enzyme=MboI

    fqdir=fastq/$cell
    custardpy_juicer -p $ncore -a $gene -b $build -g $gt \
        -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell

- ``custardpy_juicer`` assumes that the fastq files are stored in ``fastq/$cell`` (here ``fastq/Hap1-A``). The outputs are stored in ``JuicerResults_$build/$cell``.
- ``$fastq_post`` indicates the filename of input fastqs is ``*_[1|2].fastq.gz`` or ``*_[R1|R2].fastq.gz``.
- Avaible genome build: hg19, hg38, mm10, mm39, rn7, galGal5, galGal6, ce10, ce11, danRer11, dm6, xenLae2, sacCer3
- Available Enzymes: HindIII, DpnII, MboI, Sau3AI, Arima


Hi-C analysis from a .hic file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Use ``custardpy_process_hic`` command to start Hi-C analysis from a ``.hic`` file.

.. code-block:: bash

    build=hg38   # genome build
    gt=genometable.$build.txt # genome_table file
    gene=refFlat.$build.txt   # gene annotation (refFlat format)
    ncore=64  # number of CPUs
    cell=Hap1-A
    hic=sample.hic

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $cell

- The outputs are stored in ``$cell``.


Running commands separately
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can also execute commands separately. 

.. code-block:: bash

    build=hg38
    fastq_post="_R"  # "_" or "_R"  before .fastq.gz
    enzyme=MboI      # enzyme type
    norm=SCALE       # normalization type

    gt=genome_table.$build.txt  # genome_table file
    bwaindex=bwa-indexes/UCSC-$build  # BWA index file
    gene=refFlat.$build.txt # gene annotation (refFlat format)
    ncore=64 # number of CPUs

    cell=Hap1-A
    fqdir=fastq/$cell/
    odir=JuicerResults/$cell

    # generate .hic file from fastq by Juicer
    rm -rf $odir
    juicer_map.sh -p $ncore $fqdir $odir $build $gt $bwaindex $enzyme $fastq_post

    # Compress intermediate files
    juicer_pigz.sh $odir

    # plot contact frequency
    if test ! -e $odir/distance; then plot_distance_count.sh $cell $odir; fi

    hic=$odir/aligned/inter_30.hic
    # call TADs (arrowHead)
    juicer_callTAD.sh $norm $odir $hic $gt

    # call loops (HICCUPS, add '--nv' option to use GPU)
    call_HiCCUPS.sh $norm $odir $hic
    # motif analysis
    call_MotifFinder.sh $build $motifdir $odir/loops/$norm/merged_loops.bedpe

    for resolution in 25000 50000 100000
    do
        # make contact matrix for all chromosomes
        makeMatrix_intra.sh $norm $odir $hic $resolution $gt
        # calculate Eigenvector
        makeEigen.sh -p 32 $norm $odir $hic $resolution $gt $gene
        # calculate insulation score
        makeInslationScore.sh $norm $odir $resolution $gt
    done
    


Micro-C analysis by Cooler
--------------------------------------------------

Micro-C analysis by `Cooler <https://cooler.readthedocs.io/en/latest/index.html>`_ and `cooltools <https://github.com/open2c/cooltools>`_.

Micro-C using BWA
+++++++++++++++++++++++++++++++++

This command maps reads by BWA, make ``.cool`` and ``.hic`` files and call loops using Juicer.

.. code-block:: bash

    build=mm10
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    bwa_index=bwa-indexes/UCSC-$build

    prefix=ESC_WT01   # modify this for your FASTQ data
    fq1=fastq/${prefix}_1.fastq.gz
    fq2=fastq/${prefix}_2.fastq.gz

    # Generate .hic file from FASTQ
    custardpy_cooler_MicroC -t bwa -i $bwa_index -g $gt -p $ncore $fq1 $fq2 $prefix

    # Juicer analysis with the .hic file
    odir=Cooler_MicroC_bwa/$prefix
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

    
Micro-C using chromap
+++++++++++++++++++++++++++++++

**CustardPy** also supports chromap for read mapping.

.. code-block:: bash

    build=mm10
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    genome=genome.$build.fa     # genome fasta file
    chromap_index=chromap-indexes/UCSC-$build

    prefix=ESC_WT01   # modify this for your FASTQ data
    fq1=fastq/${prefix}_1.fastq.gz
    fq2=fastq/${prefix}_2.fastq.gz

    # Generate .hic file from FASTQ
    custardpy_cooler_MicroC -t chromap -i $chromap_index -g $gt -f $genome -p $ncore $fq1 $fq2 $prefix

    # Juicer analysis with the .hic file
    odir=Cooler_MicroC_chromap/$prefix
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
