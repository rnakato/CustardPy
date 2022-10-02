Visualization
===============================

plotHiCMatrix
----------------------------------------------------------------

Visualize Hi-C map as a simple square heatmap. The input data is a dense matrix output from `makeMatrix_intra.sh`.
The contact level is normalized by the total number of mapped reads for the chromosome.

.. code-block:: bash

    plotHiCMatrix <matrix> <output name (png)> <start> <end> <title in figure>

Example:

.. code-block:: bash
     chr=chr7
     start=15000000
     end=40000000
     norm=SCALE
     cell=Hap1-A
     binsize=100000
     matrix=JuicerResults_hg38/$cell/Matrix/intrachromosomal/$binsize/observed.$norm.$chr.matrix.gz

.. figure:: img/plotHiCMatrix.png
   :width: 500px
   :align: center
   :alt: Alternate

   plotHiCMatrix 


drawSquareMulti
------------------------------------------------------

Visualize contact map of multiple Hi-C samples as triangle heatmaps.
The input data is a dense matrix output from `makeMatrix_intra.sh`.

.. code-block:: bash

     # linear scale
     drawSquareMulti Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL SquareMulti.chr9 \
          chr9 --start 1000000 --end 38000000  --type VC_SQRT --vmax 20

The black dashed lines and blue circles indicate TADs and loops, respectively. 

.. figure:: img/SquareMulti.jpg
   :width: 500px
   :align: center
   :alt: Alternate

   SquareMulti

Add ``--log`` option to change the contact level to log scale.

.. code-block:: bash

     # log scale
     drawTriangleMulti Control_1:Control CTCFKD_1:siCTCF drawTriangleMulti.log.chr9 \
          chr9 --start 1000000 --end 38000000 --type VC_SQRT -d 5000000 --log

.. figure:: img/SquareMulti.log.jpg
   :width: 500px
   :align: center
   :alt: Alternate

   SquareMulti (log scale)


drawSquarePair
------------------------------------------------------

The first and second samples are visualzed in the upper and bottom triagles, respectively.

.. code-block:: bash

     drawSquarePair \
         Control/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         Rad21KD_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
         drawSquarePair.chr21 --start 24000000 --end 32000000

.. figure:: img/drawSquarePair.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   drawSquarePair

drawSquareRatioPair
------------------------------------------------------

This command visualize the log-scale frequency of ``sample2/sample1`` and ``sample4/sample3``.

.. code-block:: bash

     drawSquareRatioPair \
          Control_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
          CTCFKD_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
          Control_2/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
          Rad21KD_1/Matrix/intrachromosomal/25000/observed.VC_SQRT.chr21.matrix.gz \
          drawSquareRatioPair.chr21 --start 24000000 --end 32000000

.. figure:: img/drawSquareRatioPair.jpg
   :width: 400px
   :align: center
   :alt: Alternate

   drawSquareRatioPair

In this case, CTCFKD_1/Control_1 and Rad21KD_1/Control_2 are visualized in the upper and bottom triagles, respectively.


drawSquareRatioMulti
------------------------------------------------------

Visualize a relative contact frequency (log scale) of 2nd to the last samples against the first sample.
The input data is a dense matrix output from `makeMatrix_intra.sh`.

.. code-block:: bash

     drawSquareRatioMulti Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL drawSquareRatioMulti.chr9 \
          chr9 --start 1000000 --end 38000000 --type VC_SQRT 

.. figure:: img/drawSquareRatioMulti.jpg
   :width: 500px
   :align: center
   :alt: Alternate

   drawSquareRatioMulti


The bottom line plots are Directional frequency ratio.


drawTriangleMulti
------------------------------------------------------

Visualize contact map of multiple Hi-C samples as triangle heatmaps.
The input data is a dense matrix output from `makeMatrix_intra.sh`.

.. code-block:: bash

     # linear scale
     drawTriangleMulti Control_1:Control CTCFKD_1:siCTCF drawTriangleMulti.chr9 \
          chr9 --start 1000000 --end 38000000 --type VC_SQRT -d 5000000   

The black dashed lines and blue circles indicate TADs and loops, respectively. 

.. figure:: img/drawTriangleMulti.jpg
   :width: 600px
   :align: center
   :alt: Alternate

   drawTriangleMulti

drawTrianglePair
------------------------------------------------------

Visualize a contact frequency of the first and second sample in upper and lower triangles, respectively. 

.. code-block:: bash

     drawTrianglePair  Control_1:Control CTCFKD_1:siCTCF drawTriangleRatioMulti.chr9 \
          chr9 --start 1000000 --end 38000000 --type VC_SQRT -d 8000000

.. figure:: img/drawTrianglePair.jpg
   :width: 500px
   :align: center
   :alt: Alternate

   drawTrianglePair
   
The black dashed lines and blue circles indicate TADs and loops, respectively. 

drawTriangleRatioMulti
------------------------------------------------------

Visualize a relative contact frequency (log scale) of 2nd to the last samples against the first sample.
The input data is a dense matrix output from `makeMatrix_intra.sh`.

.. code-block:: bash

     drawTriangleRatioMulti Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL drawTriangleRatioMulti.chr9 \
          chr9 --start 1000000 --end 38000000 --type VC_SQRT -d 5000000

.. figure:: img/drawTriangleRatioMulti.jpg
   :width: 600px
   :align: center
   :alt: Alternate

   drawTriangleRatioMulti


The bottom line plots are Directional frequency ratio.


plotHiCfeature
------------------------------------------------------

Draw heatmap and line graphs for various features values of multiple Hi-C samples.

.. code-block:: bash

     plotHiCfeature [-h] [--type TYPE] [--distance DISTANCE]
                         [-r RESOLUTION] [-s START] [-e END] [--multi]
                         [--compartment] [--di] [--dfr] [--vmax VMAX]
                         [--vmin VMIN] [-d VIZDISTANCEMAX] [--xsize XSIZE]
                         [input [input ...]] output chr

``Input`` should be "<sample directory>:<label>", for instance:

.. code-block:: bash

     plotHiCfeature --type SCALE --start 1000000 --end 38000000 \
         Ctrl:Control CTCF:siCTCF InsulationScore.chr9.1M-38M 

``<sample directory>`` is the output directory by ``custardpy_juicer``.
In default, ``plotHiCfeature`` outputs a single insulation score (500 kbp distance).
``type`` is the normalization type defined by Juicer (SCALE/KR/VC_SQRT/NONE).

``plotHiCfeature`` can also output a multi-scale insulation score ranging 100 kbp to 1 Mbp by supplying ``--multi `` option.

.. code-block:: bash

     plotHiCfeature JuicerResults_hg38/Hap1-A JuicerResults_hg38/WaplKO_3.3-A \
          MultiIS.chr9.1M-38M --start 1000000 --end 38000000 \
          --multi --type SCALE -d 5000000

Other examples::

     # PC1 for compartment
     plotHiCfeature JuicerResults_hg38/Hap1-A JuicerResults_hg38/WaplKO_3.3-A \
          Compartment.chr9.1M-38M chr9 --start 1000000 --end 38000000 \
          --compartment --type SCALE -d 5000000

     # Directionality index
     plotHiCfeature JuicerResults_hg38/Hap1-A JuicerResults_hg38/WaplKO_3.3-A \
          DI.chr9.1M-38M chr9 --start 1000000 --end 38000000 \
          --di --type SCALE -d 5000000

     # Directional frequency ratio
     plotHiCfeature JuicerResults_hg38/Hap1-A JuicerResults_hg38/WaplKO_3.3-A \
          DFR.chr9.1M-38M chr9 --start 1000000 --end 38000000 \
          --dfr --type SCALE -d 5000000


plotCompartmentGenome
------------------------------------------------------

Plot a PC1 value of multiple samples for the whole genome. 

.. code-block:: bash

     plotCompartmentGenome [-h] [--type TYPE] [-r RESOLUTION] [--heatmap]
                       [input [input ...]] output
     Example:
        plotCompartmentGenome Control_1:Control CTCFKD_1:siCTCF NIPBLKD_1:siNIBPL \ 
               CompartmentGenome -r 25000 --type VC_SQRT

.. figure:: img/plotCompartmentGenome.jpg
   :width: 700px
   :align: center
   :alt: Alternate

   plotCompartmentGenome


plotInsulationScore
------------------------------------------------------

Plot a line graph of insulation score. The input data is a dense matrix output from `makeMatrix_intra.sh`.

.. code-block:: bash

     plotInsulationScore [-h] [--num4norm NUM4NORM] [--distance DISTANCE]
                                 [--sizex SIZEX] [--sizey SIZEY]
                                 matrix output resolution
     Example:
        plotInsulationScore WT/intrachromosomal/25000/observed.KR.chr7.matrix.gz InsulationScore_WT.chr7.png 25000

.. figure:: img/InsulationScore.png
   :width: 700px
   :align: center
   :alt: Alternate

   InsulationScore


plotMultiScaleInsulationScore
------------------------------------------------------

Plot multi-scale insulation scores from Juicer matrix

.. code-block:: bash

     plotMultiScaleInsulationScore [-h] [--num4norm NUM4NORM]
                                   [--sizex SIZEX] [--sizey SIZEY]
                                   matrix output resolution
     Example:
        plotInsulationScore WT/intrachromosomal/25000/observed.KR.chr7.matrix.gz MultiInsulationScore_WT.chr7.png 25000


