Installation
================

CustardPy from PyPI
---------------------------------

Core components of **CustardPy** can by installed using pip (>= Python 3.7):

.. code-block:: bash

    pip3 install custardpy


Docker images for CustardPy
---------------------------------

Most commands introduced in this manual are included in two Docker images for **CustardPy**.

- ``CustardPy``: https://hub.docker.com/r/rnakato/custardpy
    - An image that contains various tools for Hi-C/Micro-C analysis in addition to **CustardPy** itself, including:

        - Cooler version 0.8.6
        - cooltools version 0.5.1
        - HiCExplorer version 3.5.1
- ``CustardPy_Juicer``: https://hub.docker.com/r/rnakato/custardpy_juicer
    - An image for `Juicer <https://github.com/aidenlab/juicer/wiki>`_ analysis (because Juicer requires older environment: ``cuda:8.0-cudnn7-devel-ubuntu16.04``). 
    - This image internally implements the tools below:

        - Juicer version 1.6
        - Juicertools version 2.13.07
        - JuiceBox version 2.13.07
        - `BWA <http://bio-bwa.sourceforge.net/>`_ version 0.7.17

See the original website for the full description about each tool.

RUN
++++++++++++++

For Docker:

.. code-block:: bash

   # pull docker image
   docker pull rnakato/custardpy 
   docker pull rnakato/custardpy_juicer
   
   # container login
   docker run [--gpus all] --rm -it rnakato/custardpy_juicer /bin/bash
   # execute a command
   docker run --rm -it rnakato/custardpy <command>
   docker run --rm -it rnakato/custardpy_juicer <command>

For Singularity:

.. code-block:: bash

   # build image
   singularity build custardpy.sif docker://rnakato/custardpy
   singularity build custardpy_juicer.sif docker://rnakato/custardpy_juicer
   # execute a command
   singularity exec custardpy.sif <command>
   singularity exec [--nv] custardpy_juicer.sif <command>

.. note::

    ``--nv`` option allows using GPU in singularity. This option is needed only when calling loops by HiCCUPS. 