Quickstart
=====================

.. contents:: 
   :depth: 3

CustardPy using Docker or Apptainer
---------------------------------------------

Since CustardPy is provided as a Docker image, you need to use ``docker`` or ``apptainer`` commands. We recommended to use Apptainer. 

The command below is a general example of how to run CustardPy using Docker or Apptainer.
It will mount the ``/work`` directory of the host machine to the ``/work`` directory of the container.
``--gpus all`` option in docker or ``--nv`` option in Apptainer is required to use GPU in the container. 
If you do not have a GPU or do not want to use it, simply omit this option.

.. code-block:: bash

    # For docker
    $ docker run --rm -it [--gpus all] -v /work:/work rnakato/custardpy <command>
    ## Example
    $ docker run --rm -it [--gpus all] -v /work:/work rnakato/custardpy \
        custardpy_cooler -g $gt -f $genome -b $build -e $enzyme \
        -i $index_bwa -p $ncore fastq/$cell $cell

    # For Apptainer
    $ apptainer exec [--nv] --bind /work custardpy.sif <command>
    ## Example
    $ apptainer exec [--nv] --bind /work custardpy.sif \
        custardpy_cooler -g $gt -f $genome -b $build -e $enzyme \
        -i $index_bwa -p $ncore fastq/$cell $cell


See also the sample scripts in the `CustardPy Tutorial <https://github.com/rnakato/CustardPy/tree/main/tutorial>`_ on GitHub.

Available Genome Builds
-----------------------------------

CustardPy provides the reference genome files and restriction enzyme files for **Cooler/Juicer/HiC-Pro** the following genome builds.

.. csv-table::
   :class: align-center

   "**Species**", "**Builds**"
   "Human", "hg38, hg39, T2T"
   "Mouse", "mm10, mm39"
   "Rat",   "rn7"
   "Zebrafish", "danRer11"
   "Chicken", "galGal5, galGal6"
   "Xenopus_tropicalis", "xenLae2"
   "Fly", "dm6"
   "C.elegans", "ce10, ce11"
   "S.serevisiae", "sacCer3"

Hi-C analysis using Juicer
---------------------------------------------

Analysis from FASTQ files
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The ``custardpy_juicer`` command performs the Hi-C analysis from FASTQ files using `Juicer <https://github.com/aidenlab/juicer/wiki>`_.

This command includes the ``custardpy_process_hic`` command described below.

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

- ``custardpy_juicer`` assumes that the fastq files are stored in ``fastq/$cell`` (here ``fastq/Control``). The outputs are stored in ``CustardPyResults/Juicer_$build/$cell``.
- ``$fastq_post`` indicates the filename of input fastqs is ``*_[1|2].fastq.gz`` or ``*_[R1|R2].fastq.gz``.
- Available Enzymes: **HindIII, DpnII, MboI, Sau3AI, Arima, AluI**

.. _process_hic:

Analysis from a .hic file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you start the Hi-C analysis from a ``.hic`` file, use the ``custardpy_process_hic`` command.

.. code-block:: bash

    build=hg38   # genome build
    gt=genometable.$build.txt # genome_table file
    gene=refFlat.$build.txt   # gene annotation (refFlat format)
    ncore=64  # number of CPUs
    odir=outputdir
    hic=sample.hic

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

- The outputs are stored in ``$odir``.

.. note::

    Due to the backward incompatibility of `Juicertools <https://github.com/aidenlab/juicer/wiki>`_ , ``custardpy_process_hic`` fails with an error when processing ``.hic`` files created by older versions of Juicertools. In this case, try the ``-o`` option which uses older versions of Juicertools in CustardPy.


Output data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

These commands create the directories ``Eigen/``, ``InsulationScore/``, ``Matrix/``, ``TAD/``, and ``loops/`` in ``$odir`` and save the results there. 

- The ``Eigen/`` directory contains the eigenvectors and compartment structures for each chromosome.
- The ``InsulationScore/`` directory contains the insulation scores and TAD boundaries for each chromosome.
- The ``Matrix/`` directory contains the normalized contact matrices (observed and O/E).
- The ``TAD/`` directory contains the TADs called by `ArrowHead`.
- The ``loops/`` directory contains the chromatin loops called by `HiCCUPS`. If the GPU is unavailable, the chromatin loops will not be calculated and ``loops/`` directory will be blank.

The ``custardpy_juicer`` command also creates the ``distance/`` directory and saves the distance decay curve plot there.

The ``fastq/``, ``aligned/``, and ``splits/`` directories are direct outputs from **Juicer**, while several large intermediate files are automatically gzipped.


Hi-C/Micro-C analysis using Cooler
---------------------------------------------

.. note::
   Starting from version 3.0.0, ``custardpy_cooler_HiC`` and ``custardpy_cooler_MicroC`` were merged into a single command, ``custardpy_cooler``.

Analysis from FASTQ files
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Although the Juicer pipeline is widely used for Hi-C analysis, it is not suitable for Micro-C analysis.
CustardPy also allows the Hi-C and Micro-C analysis by `Cooler <https://cooler.readthedocs.io/en/latest/index.html>`_ and `cooltools <https://github.com/open2c/cooltools>`_. 

The ``custardpy_cooler`` command performs the analysis from FASTQ files, generates a ``.cool`` file, and converts it to a ``.hic`` file. 
The outputs are stored in ``CustardPyResults/Cooler_$build/$cell``.

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
    custardpy_cooler -g $gt -f $genome -b $build -e $enzyme -i $index_bwa -p $ncore fastq/$cell $cell


For the Micro-C analysis, use ``-e None`` not to specify the enzyme in the ``custardpy_cooler`` command.

.. code-block:: bash

    build=mm39
    ncore=64
    gt=genome_table.$build.txt  # genome_table file
    bwa_index=bwa-indexes/UCSC-$build
    genome=genome.$build.fa
    cell=Control   # modify this for your FASTQ data
    enzyme=None

    # Generate .hic file from FASTQ
    custardpy_cooler -g $gt -f $genome -b $build -e $enzyme -i $index_bwa -p $ncore fastq/$cell $cell


Analysis from a .cool file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The ``custardpy_process_cool`` command to analyze a ``.cool`` file created by the ``custardpy_cooler`` command or any other tools.

.. code-block:: bash

    odir=CustardPyResults/Cooler_$build/$cell
    cool=$odir/cool/cooler.dedup.q30.multires.cool
    pair=$odir/pairs/dedup.bwa.q30.pairs.gz
    $sing custardpy_process_cool -t $ncore -g $gt -p $pair -f $genome $cool $odir $cell

Analysis from a .hic file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can apply the ``custardpy_process_hic`` command to the `.hic` file (see :ref:`process_hic`).

.. code-block:: bash

    # Downstream analysis using .hic
    odir=CustardPyResults/Cooler_$build/$cell
    hic=$odir/hic/contact_map.q30.hic
    norm=SCALE

    custardpy_process_hic -p $ncore -n $norm -g $gt -a $gene $hic $odir

Output data
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The ``custardpy_process_cool`` commands create the directories ``Eigen/``, ``InsulationScore/``, ``Matrix/``, ``TAD/``, and ``loops/`` in ``$odir`` and save the results there. 

- The ``Eigen/`` directory contains the eigenvectors and compartment structures for each chromosome.
- The ``InsulationScore/`` directory contains the insulation scores and TAD boundaries for each chromosome.
- The ``Matrix/`` directory contains the normalized contact matrices (observed and O/E).
- The ``TAD/`` directory contains the TADs called by `ArrowHead`.
- The ``loops/`` directory contains the chromatin loops called by `HiCCUPS`. If the GPU is unavailable, the chromatin loops will not be calculated and ``loops/`` directory will be blank.

The ``custardpy_juicer`` command also creates the ``distance/`` directory and saves the distance decay curve plot there.

The ``fastq/``, ``aligned/``, and ``splits/`` directories are direct outputs from **Juicer**, while several large intermediate files are automatically gzipped.


Hi-C/Micro-C analysis using HiC-Pro
---------------------------------------------

CustardPy also allows the Hi-C  analysis by `HiC-Pro <https://nservant.github.io/HiC-Pro/>`_.

First, edit ``config-hicpro.txt``. 
Set the ``BOWTIE2_IDX_PATH`` and ``GENOME_SIZE`` variables to absolute paths that match your environment.

The CustardPy docker image stores the restriction files for HiC-Pro in ``/opt/restriction_sites/HiCPro-restriction_sites/``. 
You can set ``GENOME_FRAGMENT`` accordingly without preparing them. An example is shown below.

.. code-block:: bash

    BOWTIE2_IDX_PATH = /home/youraccount/CustardPy/tutorial/Hi-C/03.HiC-Pro/bowtie2-indexes
    REFERENCE_GENOME = hg38
    GENOME_SIZE = /home/youraccount/CustardPy/tutorial/Hi-C/03.HiC-Pro/genometable.hg38.txt
    GENOME_FRAGMENT = /opt/restriction_sites/HiCPro-restriction_sites/MboI_resfrag_hg38.bed


Then, run the command below to perform the Hi-C analysis from FASTQ files using HiC-Pro.
CustardPy provices the ``hicpro`` command to run HiC-Pro.

.. code-block:: bash

    fqdir=fastq/
    hicpro -i $fqdir -o HiCProResults -c config-hicpro.txt

Analysis will be run for all subdirectories in ``$fqdir``. 
The output files are the same as that of HiC-Pro. For details, please refer to the official website.