Python API
=====================

**CustardPy_Juicer** is a docker image for Juicer analysis in `CustardPy <https://github.com/rnakato/Custardpy>`_.
This is a wrapper of `Juicer <https://github.com/aidenlab/juicer/wiki>`_ and internally executes `juicertools <https://github.com/aidenlab/juicer/wiki/Feature-Annotation>`_.
See the original website for the full description about each command.

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
