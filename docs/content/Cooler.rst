Cooler_based analysis
===========================

**CustardPy** Docker image also supports the Hi-C/Micro-C analysis by `Cooler <https://cooler.readthedocs.io/en/latest/index.html>`_ and `cooltools <https://github.com/open2c/cooltools>`_.

Micro-C using BWA
-------------------------

This command maps reads by BWA, make .cool and .hic files and call loops using Juicer HiCCUPS.

.. code-block:: bash

    build=mm10
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    bwa_index=bwa-indexes/UCSC-$build

    prefix=SRR10480692
    fq1=fastq/${prefix}_1.fastq.gz
    fq2=fastq/${prefix}_2.fastq.gz

    idir=Results_bwa/$prefix
    singularity exec --bind /work custardpy.sif \
        microc_bwa -p $ncore $prefix $idir $fq1 $fq2 $bwa_index $gt

    hic=$idir/hic/contact_map.q30.hic
    norm=SCALE
    singularity exec --nv custardpy_juicer.sif \
        call_HiCCUPS.sh $norm $idir $hic $build

    
Micro-C using chromap
-------------------------

**CustardPy** also allows chromap for read mapping.

.. code-block:: bash

    build=mm10
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    genome=genome.$build.fa     # genome fasta file
    chromap_index=chromap-indexes/UCSC-$build

    sing=""
    sing_juicer=""

    prefix=SRR10480692
    fq1=fastq/${prefix}_1.fastq.gz
    fq2=fastq/${prefix}_2.fastq.gz

    idir=Results_chromap/$prefix
    singularity exec custardpy.sif \
        microc_chromap -p $ncore $prefix $idir $fq1 $fq2 $chromap_index $gt $genome

    hic=$idir/hic/contact_map.q30.hic
    norm=SCALE
    singularity exec --nv custardpy_juicer.sif call_HiCCUPS.sh $norm $idir $hic $build
