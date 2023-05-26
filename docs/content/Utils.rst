Utility tools
=====================

This page describes the utility scripts included in **CustardPy**.

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

Convert interaction frequency file dumped by Juicertools to dense (two-dimensional) matrix.

.. code-block:: bash

    convert_JuicerDump_to_dense.py <inputfile> <outputfile> <genometable> <chr> <resolution> [--help]

    Example:
    convert_JuicerDump_to_dense.py \
        Matrix.observed.VC_SQRT.chrX.txt \
        Matrix.observed.VC_SQRT.chrX.matrix.gz \
        genome_table.txt \
        chrX \
        25000
