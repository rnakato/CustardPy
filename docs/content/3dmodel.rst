3D/4D modeling form Hi-C data
====================================

**CustardPy** has an advanced command for 3D modeling and 4D simulatuion using the Hi-C matrix.


run-pastis.sh
---------------------------

``run-pastis.sh`` implements `PASTIS <https://members.cbio.mines-paristech.fr/~nvaroquaux/pastis/>`_ and estimates 3D models from the Hi-C contact matrix.

.. code-block:: bash

    odir=JuicerResults/Sample
    # genomic region to be simulated
    chr=chr21
    s=24000000
    e=32000000
    resolution=100000
    norm=SCALE  # normalization type
    
    run-pastis.sh $odir $chr $s $e $resolution $norm

phic_processing
------------------------------------------------------

``phic_processing`` implements `PHi-C2 <https://github.com/soyashinkai/PHi-C2>`_ and generates biophysical polymer model.

.. code-block:: bash

    hic=$odir/aligned/inter_30.hic
    odir=phic_result

    # genomic region to be simulated
    chr=chr21
    start=24000000
    end=32000000
    resolution=100000
    norm=SCALE # normalization type

    phic_processing $odir $hic $chr $start $end $resolution $norm
    