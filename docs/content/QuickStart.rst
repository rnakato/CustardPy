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

    # Example of custardpy_juicer
    # For docker
    docker run --rm -it --gpus all -v /work:/work rnakato/custardpy \
        custardpy_juicer -p $ncore -a $gene -b $build -g $gt \
        -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell

    # For singularity
    singularity exec --nv --bind /work custardpy.sif \
        custardpy_juicer -p $ncore -a $gene -b $build -g $gt \
        -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell

See also the sample scripts in the `tutorial <https://github.com/rnakato/CustardPy/tree/main/tutorial>`_ on GitHub.


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

    cell=Control
    fastq_post="_"  # "_" or "_R"
    enzyme=MboI

    fqdir=fastq/$cell
    custardpy_juicer -p $ncore -a $gene -b $build -g $gt \
        -i $bwaindex -e $enzyme -z $fastq_post $fqdir $cell

- ``custardpy_juicer`` assumes that the fastq files are stored in ``fastq/$cell`` (here ``fastq/Control``). The outputs are stored in ``CustardPyResults_Hi-C/Juicer_$build/$cell``.
- ``$fastq_post`` indicates the filename of input fastqs is ``*_[1|2].fastq.gz`` or ``*_[R1|R2].fastq.gz``.
- Avaible genome build: hg19, hg38, mm10, mm39, rn7, galGal5, galGal6, ce10, ce11, danRer11, dm6, xenLae2, sacCer3
- Available Enzymes: HindIII, DpnII, MboI, Sau3AI, Arima, AluI


Hi-C analysis from a .hic file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you start the Hi-C analysis from a ``.hic`` file, use ``custardpy_process_hic`` command.

.. code-block:: bash

    build=hg38   # genome build
    gt=genometable.$build.txt # genome_table file
    gene=refFlat.$build.txt   # gene annotation (refFlat format)
    ncore=64  # number of CPUs
    cell=Control
    hic=sample.hic

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $cell

- The outputs are stored in ``$cell``.

.. note::

    Due to the backward incompatibility of Juicertools, ``custardpy_process_hic`` fails with an error when processing .hic files created by older Juicertools. In this case, use the ``-o`` option which uses older versions of Juicertools in CustardPy.


Hi-C analysis using Cooler
---------------------------------------------

CustardPy allows the Hi-C analysis by `Cooler <https://cooler.readthedocs.io/en/latest/index.html>`_ and `cooltools <https://github.com/open2c/cooltools>`_. 
``custardpy_cooler_HiC`` generates a ``.cool`` file and converts it to a ``.hic`` file. You can apply ``custardpy_process_hic`` command to it.
The outputs are stored in ``CustardPyResults_MicroC/Cooler_$build//$cell``.

.. code-block:: bash

    build=hg38
    gt=genometable.hg38.txt
    index_bwa=bwa-indexes/hg38
    gene=refFlat.$build.txt
    genome=genome.$build.fa
    ncore=64

    cell=Control
    enzyme=MboI

    # Generate .cool and .hic files from FASTQ
    custardpy_cooler_HiC -g $gt -b $build -f $genome -i $index_bwa -p $ncore fastq/$cell $cell

    # Downstream analysis using .hic
    odir=CustardPyResults_cooler/$build/$cell
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE
    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir


    
Micro-C analysis using Cooler
--------------------------------------------------

Micro-C analysis by `Cooler <https://cooler.readthedocs.io/en/latest/index.html>`_ and `cooltools <https://github.com/open2c/cooltools>`_.

Micro-C using BWA
+++++++++++++++++++++++++++++++++

The command ``custardpy_cooler_MicroC`` maps Micro-C reads by BWA and makes ``.cool`` and ``.hic`` files. The ``.hic`` file is processed using ``custardpy_process_hic``.

.. code-block:: bash

    build=mm39
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    bwa_index=bwa-indexes/UCSC-$build
    genome=genome.$build.fa
    cell=C36_rep1   # modify this for your FASTQ data

    # Generate .hic file from FASTQ
    custardpy_cooler_MicroC -t bwa -g $gt -f $genome -i $bwa_index -p $ncore fastq/$cell $cell

    # Juicer analysis with the .hic file
    odir=CustardPyResults_MicroC/Cooler_bwa/$cell
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

- ``custardpy_cooler_MicroC`` assumes that the fastq files are stored in ``fastq/$cell`` (here ``fastq/C36_rep1``). The outputs are stored in ``CustardPyResults_MicroC/Cooler_bwa/$cell``.
    
.. Micro-C using chromap
.. +++++++++++++++++++++++++++++++

.. **CustardPy** also supports chromap for read mapping.

.. .. code-block:: bash
.. 
..     build=mm10
..     ncore=64
..     gt=genome_table.$build.txt  # genome_table file
..     genome=genome.$build.fa     # genome fasta file
..     chromap_index=chromap-indexes/UCSC-$build

..     cell=ESC_WT01   # modify this for your FASTQ data

..     # Generate .hic file from FASTQ
..     custardpy_cooler_MicroC -t chromap -i $chromap_index -g $gt -f $genome -p $ncore fastq/$cell $cell

..     # Juicer analysis with the .hic file
..     odir=CustardPyResults_MicroC/$cell/chromap
..     hic=$odir/hic/contact_map.q30.hic
..     norm=SCALE
..     custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir
