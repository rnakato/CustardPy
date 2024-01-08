================================================================
CustardPy
================================================================

**CustardPy** is a 3D genome analysis tool written in Python3 and available using the `Docker system <https://www.docker.com/>`_.
It is mainly designed for multi-sample Hi-C analysis (e.g., comparison of depletion effects across multiple proteins) and provides various functions for analysis and visualization.

Since the docker image includes various existing tools in addition to the core **CustardPy** component (see :doc:`./content/Install`),
users can perform a variety of 3D genome analyses without having to install them individually.

- **Major Release! (version 1)**
   - Unified the Docker images for `CustardPy <https://hub.docker.com/r/rnakato/custardpy>`_ and `CustardPy_Juicer <https://hub.docker.com/r/rnakato/custardpy_juicer>`_. Version 1 of the `CustardPy <https://hub.docker.com/r/rnakato/custardpy>`_ docker image now supports all analyses previously offered by CustardPy and CustardPy_Juicer, rendering the latter unnecessary.

Contents:
---------------

.. toctree::
   :numbered:
   :glob:
   :maxdepth: 1

   content/Install
   content/QuickStart
   content/StepbyStep
   content/Visualization
   content/DEGanalysis
   content/Multisample
   content/3DChromatin_ReplicateQC
   content/3dmodel
   content/Command

Citation:
--------------

* Nakato R, Sakata T, Wang J, Nagai LAE, Nagaoka Y, Oba GM, Bando M, Shirahige K, Context-dependent perturbations in chromatin folding and the transcriptome by cohesin and related factors, *Nature Communications*, 2023. doi: 10.1038/s41467-023-41316-4
