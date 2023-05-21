Utility tools
=====================

distance_vs_count.Juicer
---------------------------

Count genomic distance of read pairs in the input file (supposing ``align/merged_nodups.txt.gz`` in Juicer outputs.)

.. code-block:: bash

    distance_vs_count.Juicer <file> <winsize> <MAPQ>
        <file>:    Input file  (merged_nodups.txt.gz)
        <winsize>: window size (default: 10000)
        <MAPQ>:    MAPQ threshold (default: 30)


convert_JuicerDump_to_dense.py
------------------------------------------------------

Convert interaction frequency file dumped by Juicer to dense matrix

.. code-block:: bash

    convert_JuicerDump_to_dense.py <inputfile> <outputfile> <genometable> <chr> <resolution> [--help]

    Example:
    convert_JuicerDump_to_dense.py \
        Matrix.observed.VC_SQRT.chrX.txt \
        Matrix.observed.VC_SQRT.chrX.matrix.gz \
        genome_table.txt \
        chrX \
        25000

DirectionalFreqRatio.py
------------------------------------------------------

Output  DirectionalFreqRatio as bedGraph

.. code-block:: bash

    DirectionalFreqRatio.py [-h] [--dfr_right] [--dfr_left] input control output chr resolution

    Example:
        DirectionalFreqRatio.py \
            Ctrl/matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matri x.gz \
            siCTCF/matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
            CTCF.DFR.chr21 \
            chr21 \
            25000
