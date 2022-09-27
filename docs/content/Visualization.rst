Visualization
===============================

plotHiCMatrix
----------------------------------------------------------------

Visualize Hi-C map as a simple square heatmap. The input data is a dense matrix output from `makeMatrix_intra.sh`.
The contact level is normalized by the total number of mapped reads for the chromosome.

.. code-block:: bash

    plotHiCMatrix <matrix> <output name (png)> <start> <end> <title in figure>

    Example:
        plotHiCMatrix WT/intrachromosomal/25000/observed.KR.chr7.matrix.gz hetamap_WT.chr7.25000000-31000000.png 25000000 31000000 WT


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

    plotHiCfeature Ctrl:Control CTCF:siCTCF MultiIS.chr9.1M-38M --multi \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000

Other examples::

     # PC1 for compartment
     plotHiCfeature Ctrl:Control CTCF:siCTCF Compartment.chr9.1M-38M --compartment \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000

     # Directionality index
     plotHiCfeature Ctrl:Control CTCF:siCTCF DI.chr9.1M-38M --di \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000

     # Directional frequency ratio
     plotHiCfeature Ctrl:Control CTCF:siCTCF DFR.chr9.1M-38M --dfr \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000

     # DirectionalFreqRatio (right)
     plotHiCfeature Ctrl:Control CTCF:siCTCF DFRright.chr9.1M-38M --dfr_right \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000

     # DirectionalFreqRatio (left)
     plotHiCfeature Ctrl:Control CTCF:siCTCF DFRleft.chr9.1M-38M --dfr_left \
          chr9 --start 1000000 --end 38000000 --type SCALE -d 5000000


plotInsulationScore.py
------------------------------------------------------

Plot insulation score from Juicer matrix

.. code-block:: bash

     plotInsulationScore.py [-h] [--num4norm NUM4NORM] [--distance DISTANCE]
                                 [--sizex SIZEX] [--sizey SIZEY]
                                 matrix output resolution

plotMultiScaleInsulationScore.py
------------------------------------------------------

Plot multi-scale insulation scores from Juicer matrix

.. code-block:: bash

     plotMultiScaleInsulationScore.py [-h] [--num4norm NUM4NORM]
                                           [--sizex SIZEX] [--sizey SIZEY]
                                           matrix output resolution

drawTriangleMulti
------------------------------------------------------

Plot Interaction matrix from Juicer matrix

.. code-block:: bash

     drawTriangleMulti [input [input ...]] output region --type $type        # linear scale
     drawTriangleMulti [input [input ...]] output region --type $type --log  # log scale

drawTriangleRatioMulti
------------------------------------------------------

Plot Interaction ratio (from 2nd to the last samples divided by 1st sample).

.. code-block:: bash

     drawTriangleRatioMulti [input [input ...]] output region --type $type    # logratio

