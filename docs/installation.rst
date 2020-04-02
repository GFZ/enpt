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

    # clone the source code of enpt
    git clone https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT.git
    cd EnPT

    # avoid package incompatibilities
    - conda config --set channel_priority strict

    # install some enpt dependencies that may cause trouble when installed via pip
    conda install -c conda-forge scipy

    # install not pip-installable deps of geoarray
    conda install -c conda-forge numpy scikit-image matplotlib pandas gdal pyproj basemap shapely

    # install not pip-installable deps of sensormapgeo
    conda install -c conda-forge pyresample

    # install sicor
    conda install -c conda-forge pygrib h5py pytables pyfftw numba llvmlite scikit-learn
    rm -rf context/sicor
    git clone https://gitext.gfz-potsdam.de/EnMAP/sicor.git ./context/sicor
    cd ./context/sicor
    make install
    cd ../../

    # install enpt
    make install
    cd ..
    pwd
    ls


This is the preferred method to install EnPT, as it will always install the most recent stable release. 

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
