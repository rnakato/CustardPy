Quality check, reproducibility, and comparative analysis
============================================================

**CustardPy** provides the script to run `3DChromatin_ReplicateQC <https://github.com/kundajelab/3DChromatin_ReplicateQC>`_, which checks for quality and reproducibility. 
This reproducibility analysis can also be used to compare the overall similarity among Hi-C samples.

``3DChromatin_ReplicateQC`` runs:

- `QuASAR <http://github.com/bxlab/hifive>`_, 
- `HiCRep <http://github.com/qunhualilab/hicrep>`_,
- `GenomeDISCO <http://github.com/kundajelab/genomedisco>`_, and
- `HiC-Spector <http://github.com/gersteinlab/HiC-spector>`_

for quality check and similarity calculation.
See the original website for the detailed usage.

Here is a tutorial showing how to use 3DChromatin_ReplicateQC in CustardPy.
The data is the same as :doc:`StepbyStep`.
You can also use the ``07.QualityCheck.sh`` script in the `tutorial <https://github.com/rnakato/CustardPy/tree/main/tutorial>`_ on GitHub.


Step-by-step tutorial to run 3DChromatin_ReplicateQC in CustardPy
-------------------------------------------------------------------

Set the parameters
+++++++++++++++++++++++++++

Here we compare three Hi-C samples (Control, siCTCF, siRad21).
To save time, this tutorial only uses chromosomes 21 and 22 for calculation.
All results will be output to ``$outputdir``.

.. code-block:: bash

  build=hg38
  gt=genometable.$build.txt

  outputdir=3DChromatin_ReplicateQC  # output directory
  mkdir -p $outputdir

  samples="Control siCTCF siRad21" # Samples to be compared

  chrs="chr21 chr22" # chromosomes to be considered
  resolution=50000
  norm=SCALE

Prepare the metadata and input data
++++++++++++++++++++++++++++++++++++++++++

Create ``metadata.pairs``.

.. code-block:: bash

  pairlist=$outputdir/metadata.pairs
  rm -rf $pairlist
  for sample1 in $samples; do
      for sample2 in $samples; do
          echo -e $sample1"\t"$sample2 >> $pairlist
      done
  done

Generate contact data for all samples.

.. code-block:: bash

  rm -rf $outputdir/data
  mkdir -p $outputdir/data
  for cell in $samples
  do
      hic=CustardPyResults_Hi-C/Juicer_hg38/$cell/aligned/inter_30.hic

      echo "preparing $cell..."
      for chr in $chrs; do
          $sing juicertools.sh dump observed $norm $hic $chr $chr BP $resolution \
          | awk -v chr=$chr 'OFS="\t" {printf("%s\t%d\t%s\t%d\t%d\n", chr, $1, chr, $2, $3)}' \
          | grep -v NaN > $outputdir/data/$cell.$chr.txt
      done

      cat $outputdir/data/$cell.*.txt > $outputdir/data/$cell.res$resolution
      $sing pigz $outputdir/data/$cell.res$resolution
      rm $outputdir/data/$cell.*.txt
  done

Generate ``metadata.samples``.

.. code-block:: bash

  samplelist=$outputdir/metadata.samples
  rm -rf $samplelist
  for cell in $samples; do
      echo -e "$cell\t$(pwd)/$outputdir/data/$cell.res$resolution" >> $samplelist
  done

Generate the Bin list.

.. code-block:: bash

  binlist=$outputdir/data/Bins.$resolution.bed
  rm -rf $binlist
  for chr in $chrs; do
      $sing generate_binlist_from_gtfile.py $gt $chr $resolution >> $binlist
  done
  gzip -f $binlist


Run 3DChromatin_ReplicateQC
+++++++++++++++++++++++++++++++++++++++

``run_3DChromatin_ReplicateQC.sh run_all`` run all tools and output the results in ``$outputdir/output``.

.. code-block:: bash

  run_3DChromatin_ReplicateQC.sh run_all \
      --metadata_samples $samplelist --bins $binlist.gz --metadata_pairs $pairlist --outdir $outputdir/output

Plot figures from the output
+++++++++++++++++++++++++++++++++++++++

``visualize_QC.py`` plots figures for each tool. The pdf files are output to ``3DChromatin_ReplicateQC/pdf``.

.. code-block:: bash

  visualize_QC.py 3DChromatin_ReplicateQC/