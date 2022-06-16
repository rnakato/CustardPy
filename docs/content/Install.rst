Installation
================

CustardPy from PyPI
---------------------------------

Core components of **CustardPy** can by installed using pip (>= Python 3.7):

.. code-block:: bash

    pip3 install custardpy


Docker images for CustardPy
---------------------------------

**CustardPy** has two Docker images for Hi-C/Micro-C analysis.

- ``CustardPy``: An image that contains additional scripts for Hi-C/Micro-C analysis.
- ``CustardPy_Juicer``: An image for `Juicer <https://github.com/aidenlab/juicer/wiki>`_ analysis (because Juicer requires older environment: ``cuda:8.0-cudnn7-devel-ubuntu16.04``).

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