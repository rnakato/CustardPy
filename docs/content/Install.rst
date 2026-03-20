Installation
================

CustardPy Docker image
---------------------------------

Docker image of **CustardPy** is available at `DockerHub <https://hub.docker.com/r/rnakato/custardpy>`_.
This image contains various tools for 3D genome analysis in addition to **CustardPy** core components as below:

- Mapping
    - `BWA <http://bio-bwa.sourceforge.net/>`_
    - `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ 
    - `chromap <https://github.com/haowenz/chromap>`_ 

- Hi-C/Micro-C analysis
    - `Juicer <https://github.com/aidenlab/juicer/wiki>`_ (v1.6)
    - `Juicertools <https://github.com/aidenlab/juicer/wiki>`_ 
    - `JuiceBox <https://github.com/aidenlab/Juicebox>`_ 
    - `Cooler <https://cooler.readthedocs.io/en/stable/>`_
    - `cooltools <https://cooltools.readthedocs.io/en/latest/>`_
    - `Pairtools <https://pairtools.readthedocs.io/en/latest/>`_ 
    - `coolpup.py <https://github.com/open2c/coolpuppy>`_
    - `HiCExplorer <https://hicexplorer.readthedocs.io/en/latest/>`_
    - `HOMER <http://homer.ucsd.edu/homer/interactions/index.html>`_
    - `FAN-C <https://fan-c.readthedocs.io/en/latest/index.html>`_
    - `HiC-Pro <https://github.com/nservant/HiC-Pro>`_ 
    - `HiCUP <https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/>`_
    - `HiC1Dmetrics <https://h1d.readthedocs.io/en/latest/>`_ 
    - `CALDER2 <https://github.com/CSOgroup/CALDER2>`_

- Loop calling
    - `FitHiC <https://github.com/ay-lab/fithic>`_
    - `CHiCAGO <https://bitbucket.org/chicagoTeam/chicago/src/master/>`_
 
- Stripe analysis
    - `STRIPENN <https://github.com/VahediLab/stripenn>`_

- Chromatin hub analysis
    - `FIREcaller <https://github.com/yycunc/FIREcaller>`_

- Sample comparison
    - `GENOVA <https://github.com/robinweide/GENOVA>`_
    - `CHESS <https://chess-hic.readthedocs.io/en/latest/index.html>`_ 

- 3D/4D modeling
    - `PASTIS <https://github.com/hiclib/pastis>`_ 
    - `PHi-C2 <https://github.com/soyashinkai/PHi-C2>`_ 

- Quality check
    - `3DChromatin_ReplicateQC <https://github.com/kundajelab/3DChromatin_ReplicateQC>`_
    - `MultiQC <https://multiqc.info/>`_

- Hi-ChIP
    - `FitHiChIP <https://ay-lab.github.io/FitHiChIP/html/index.html>`_ 

- ChIA-PET
    - `Mango <https://github.com/dphansti/mango>`_
    - `ChIAPop <https://github.com/wh90999/ChIAPoP>`_

- Genome analysis
    - `MACS2 <https://github.com/macs3-project/MACS>`_ 

- File processing
   - `SAMtools <http://www.htslib.org/>`_
   - `BEDtools <https://bedtools.readthedocs.io/en/latest/>`_ 
   - `hictk <https://hictk.readthedocs.io/en/stable/>`_ 

- Utility tools
   - `SRAtoolkit <https://github.com/ncbi/sra-tools>`_
   - `genomepy <https://vanheeringen-lab.github.io/genomepy/>`_ 

For a full description of each tool, visit the original website.

Download and run the Docker image
+++++++++++++++++++++++++++++++++++++++++++++++

For Docker:

.. code-block:: bash

   # pull docker image
   docker pull rnakato/custardpy:<version>
   
   # container login
   docker run --rm -it rnakato/custardpy:<version> /bin/bash

   # execute a command
   docker run --rm -it rnakato/custardpy:<version> <command>

For Apptainer:  

.. code-block:: bash

   # build image
   apptainer build custardpy.<version>.sif docker://rnakato/custardpy:<version>

   # execute a command
   apptainer exec custardpy.<version>.sif <command>


CustardPy from PyPI
---------------------------------

Core components of **CustardPy** (e.g., commands for visualization) can by installed using pip (>= Python 3.7):

.. code-block:: bash

    pip3 install custardpy


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

