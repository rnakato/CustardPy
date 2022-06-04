Installation
================

Docker images for CustardPy
---------------------------------

CustardPy contains two Docker images, ``CustardPy`` and ``CustardPy_Juicer``.

For Docker:

.. code-block:: bash

   # pull docker image
   docker pull rnakato/custardpy 
   docker pull rnakato/custardpy_juicer
   
   # container login
   docker run [--gpus all] --rm -it rnakato/custardpy_juicer /bin/bash
   # execute a command
   docker run [--gpus all] --rm -v (your directory):/opt/work rnakato/custardpy_juicer <command>

For Singularity:

.. code-block:: bash

   # build image
   singularity build custardpy.sif docker://rnakato/custardpy
   singularity build custardpy_juicer.sif docker://rnakato/custardpy_juicer
   # execute a command
   singularity exec [--nv] custardpy_juicer.sif <command>


CustardPy
---------------------------------

Stand-alone CustardPy can be installed by pip (>= Python 3.7):

.. code-block:: bash

    pip3 install custardpy

