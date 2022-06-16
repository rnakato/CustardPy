Python API
=====================

you can import the module from a Python console or script:

.. code-block:: Python3

    from custardpy.HiCmodule import *
    from custardpy.PlotModule import *

The following tutorials will assume that you have loaded the fanc module in this manner.

Loading datasets
-----------------------

Here we specify files of the contact matrix and eigenvector (PC1) that are stored in the directory ``JuicerResults/Hap1-A``.

.. code-block:: Python3

    from custardpy.HiCmodule import *
    from custardpy.PlotModule import *

    normalizetype = "RPM" # total read normalization
    resolution = 100000
    matrix = "JuicerResults/Hap1-A/Matrix/intrachromosomal/" + str(resolution) + "/observed.SCALE.chr21.matrix.gz"
    eigenfile = "JuicerResults/Hap1-A/Eigen/" + str(resolution) + "/eigen.SCALE.chr21.txt.gz"
    data = JuicerMatrix(normalizetype, matrix, resolution, eigenfile)

Then you can visualize the matrix:

.. code-block:: Python3

    s = 20000000
    e = 107000000

    mat = Jdata[sample].getmatrix().loc[s:e,s:e]
    plt.imshow(mat, clim=(0, 30), cmap=generate_cmap(['#FFFFFF', '#d10a3f']))
    plt.title(sample)


Quickstart
---------------------------

.. code-block:: Python3

    from custardpy.HiCmodule import *
    from custardpy.PlotModule import *

    def loadJuicerData(dirname, type, resolution):
        return JuicerMatrix(normalizetype, matrixfile, resolution, eigenfile)
      
    s = 20000000
    e = 107000000

    plt.figure(figsize=(12, 9))
    for i, sample in enumerate(samples):
        ax = plt.subplot(3, 4, i+1)
        mat = Jdata[sample].getmatrix().loc[s:e,s:e]
        ax.imshow(mat, clim=(0, 30), cmap=generate_cmap(['#FFFFFF', '#d10a3f']))
        ax.set_title(sample)
        axxticks(ax, 0, mat.shape[1], s, e, 5)
        ax.set_yticks([])
    plt.tight_layout()
