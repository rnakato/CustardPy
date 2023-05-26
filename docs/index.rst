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
   content/Hi-C
   content/Visualization
   content/DEGanalysis
   content/Multisample
   content/3dmodel
   content/Utils

Citation:
--------------

* Nakato R, Sakata T, Wang J, Nagai LAE, Oba GM, Bando M, Shirahige K, Context-dependent 3D genome regulation by cohesin and related factors, *bioRxiv*, 2022. doi: 10.1101/2022.05.24.493188

Contact:
--------------

rnakato AT iqb.u-tokyo.ac.jp