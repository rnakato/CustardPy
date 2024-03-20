Installation
================

CustardPy Docker image
---------------------------------

.. note::

    From version 1, the `CustardPy <https://hub.docker.com/r/rnakato/custardpy>`_ docker image supports all analyses previously offered by CustardPy and `CustardPy_Juicer <https://hub.docker.com/r/rnakato/custardpy_juicer>`_ images, rendering the latter unnecessary.

Docker image of **CustardPy** is available at `DockerHub <https://hub.docker.com/r/rnakato/custardpy>`_.
This image contains various tools for Hi-C/Micro-C analysis in addition to **CustardPy** core components as below:

- Mapping

    - `BWA <http://bio-bwa.sourceforge.net/>`_ v0.7.17
    - `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ v2.4.5
    - `chromap <https://github.com/haowenz/chromap>`_ v0.2.6

- Hi-C/Micro-C analysis

    - `Juicer <https://github.com/aidenlab/juicer/wiki>`_ v1.6
    - `Juicertools <https://github.com/aidenlab/juicer/wiki>`_ v2.13.07
    - `JuiceBox <https://github.com/aidenlab/Juicebox>`_ v2.13.07
    - `Cooler <https://cooler.readthedocs.io/en/stable/>`_ v0.9.1
    - `cooltools <https://cooltools.readthedocs.io/en/latest/>`_ v0.5.4
    - `Pairtools <https://pairtools.readthedocs.io/en/latest/>`_ v1.0.2
    - `coolpup.py <https://github.com/open2c/coolpuppy>`_ v1.0.0
    - `HiCExplorer <https://hicexplorer.readthedocs.io/en/latest/>`_ v3.5.1
    - `HOMER <http://homer.ucsd.edu/homer/interactions/index.html>`_
    - `FAN-C <https://fan-c.readthedocs.io/en/latest/index.html>`_ v0.9.25
    - `HiC-Pro <https://github.com/nservant/HiC-Pro>`_ v3.1.0
    - `HiCUP <https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/>`_ v0.9.2
    - `HiC1Dmetrics <https://h1d.readthedocs.io/en/latest/>`_ v0.2.9
    - `CALDER2 <https://github.com/CSOgroup/CALDER2>`_ v2.0

- Loop calling
    - `FitHiC <https://github.com/ay-lab/fithic>`_ v2.0.7
    - `CHiCAGO <https://bitbucket.org/chicagoTeam/chicago/src/master/>`_ v1.19.0

- Stripe analysis
    - `STRIPENN <https://github.com/VahediLab/stripenn>`_ v1.1.65.18

- Chromatin hub analysis
    - `FIREcaller <https://github.com/yycunc/FIREcaller>`_ v1.42

- Sample comparison
    - `GENOVA <https://github.com/robinweide/GENOVA>`_ v1.0.1
    - `CHESS <https://chess-hic.readthedocs.io/en/latest/index.html>`_ v0.3.7

- 3D/4D modeling
    - `PASTIS <https://github.com/hiclib/pastis>`_ v0.4.0
    - `PHi-C2 <https://github.com/soyashinkai/PHi-C2>`_ v2.0.10

- Quality check
    - `3DChromatin_ReplicateQC <https://github.com/kundajelab/3DChromatin_ReplicateQC>`_ v1.0.1
    - `MultiQC <https://multiqc.info/>`_ v

- Hi-ChIP
    - `FitHiChIP <https://ay-lab.github.io/FitHiChIP/html/index.html>`_ v11.0

- ChIA-PET
    - `Mango <https://github.com/dphansti/mango>`_
    - `ChIAPop <https://github.com/wh90999/ChIAPoP>`_

- Genome analysis
    - `MACS2 <https://github.com/macs3-project/MACS>`_ v2.2.9.1

- File processing
   - `SAMtools <http://www.htslib.org/>`_ v1.19.2
   - `BEDtools <https://bedtools.readthedocs.io/en/latest/>`_ v2.31.0

- Utility tools
   - `SRAtoolkit <https://github.com/ncbi/sra-tools>`_ v3.0.10
   - `genomepy <https://vanheeringen-lab.github.io/genomepy/>`_ v0.16.1

For a full description of each tool, visit the original website.


RUN
++++++++++++++

For Docker:

.. code-block:: bash

   # pull docker image
   docker pull rnakato/custardpy
   
   # container login
   docker run [--gpus all] --rm -it rnakato/custardpy /bin/bash

   # execute a command
   docker run [--gpus all] --rm -it rnakato/custardpy <command>

For Singularity:

.. code-block:: bash

   # build image
   singularity build custardpy.sif docker://rnakato/custardpy

   # execute a command
   singularity exec [--nv] custardpy.sif <command>

.. note::

    ``--gpus all`` for Docker and ``--nv`` option for Singularity allow using GPU. This option is needed only when calling loops by HiCCUPS.

CustardPy from PyPI
---------------------------------

Core components of **CustardPy** (e.g., commands for visualization) can by installed using pip (>= Python 3.7):

.. code-block:: bash

    pip3 install custardpy
