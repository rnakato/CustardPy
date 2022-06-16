MicroC analysis
=====================

**CustardPy** docker image contains various tools for Hi-C/Micro-C analysis in addition to **CustardPy** itself, including:

- Cooler version 0.8.6
- cooltools version 0.5.1
- HiCExplorer version 3.5.1

See the original website for the full description about each tool.

microc_bwa
----------------------------------------------------------------

``microc_bwa`` generates ``.hic`` file from fastq using BWA.

The input fastq files can be gzipped (.fastq.gz).
The results including the ``.hic`` file is outputted in ``$odir``.
The BWA index file is needed.

.. code-block:: bash

  microc_bwa [options] <prefix> <outputdir> <fq1> <fq2> <bwaindex> <genometable>
    prefix: prefix of output files
    outputdir: directory where the output files are generaged
    fq1, fq2: input fastq files (paired-end, .fastq(.gz))
    bwaindex: index file for bwa
    genometable: genome table file (describing the chromosome length)

  Options:
    -p ncore : number of CPUs (default: 4)
    
  Example:
      microc_bwa -p 24 Hap1-A Results_bwa/Hap1-A Hap1A_1.fastq.gz Hap1A_2.fastq.gz index_bwa/UCSC-hg38 genometable.hg38.txt

- Output

    - 4-minus.cool
    - 4-minus.multires.cool
    - bam/
    - hic/
    - log/
    - loops/
    - pairs/
    - pairs.stats.txt
    - qc_report/


microc_chromap
-----------------------------------------------------------------

``microc_chromap`` generates ``.hic`` file from fastq using chromap.

The input fastq files can be gzipped (.fastq.gz).
The results including the ``.hic`` file is outputted in ``$odir``.
The BWA index file is needed.

.. code-block:: bash

  microc_chromap [options] <prefix> <outputdir> <fq1> <fq2> <chromapindex> <genometable> <genome>
    prefix: prefix of output files
    outputdir: directory where the output files are generaged
    fq1, fq2: input fastq files (paired-end, .fastq(.gz))
    chromapindex: index file for chromap
    genometable: genome table file (describing the chromosome length)
    genome: genome file (.fa)

  Options:
    -p ncore : number of CPUs (default: 4)
    
  Example:
      microc_chromap -p 24 Hap1-A Results_chromap/Hap1-A Hap1A_1.fastq.gz Hap1A_2.fastq.gz index_chromap/genometable.hg38.txt genome.hg38.fa

- Output

    - 4-minus.cool
    - 4-minus.multires.cool
    - hic/
    - log/
    - loops/
    - pairs/
    - qc_report/