.. _installation:

============
Installation
============

EnPT depends on some open source packages which are usually installed without problems by the automatic install
routine. However, for some projects, we strongly recommend resolving the dependency before the automatic installer
is run. This approach avoids problems with conflicting versions of the same software.
Using conda_, the recommended approach is:

.. _conda: https://conda.io/docs/

.. code-block:: bash

    # create virtual environment for enpt, this is optional but recommended
    conda create --name enpt python=3
    source activate enpt

    # avoid package incompatibilities
    conda config --set channel_priority strict

    # install some dependencies that cause trouble when installed via pip
    conda install -c conda-forge numpy pandas lxml

    # install not pip-installable deps of py_tools_ds / geoarray / sensormapgeo
    conda install -c conda-forge gdal libgdal scikit-image pyproj geopandas matplotlib basemap shapely pyresample

    # install not pip-installable deps of arosics
    conda install -c conda-forge pyfftw pykrige

    # install not pip-installable deps of sicor
    conda install -c conda-forge glymur cachetools pyhdf h5py pytables llvmlite numba scikit-learn

    # install enpt
    pip install enpt


This is the preferred method to install EnPT, as it will always install the most recent stable release. 

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
