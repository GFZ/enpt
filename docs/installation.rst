.. _installation:

============
Installation
============


Using Anaconda or Miniconda (recommended)
-----------------------------------------

Using conda_ (latest version recommended), EnPT is installed as follows:


1. Create virtual environment for enpt (optional but recommended):

   .. code-block:: bash

    $ conda create -c conda-forge --name enpt python=3
    $ conda activate enpt


2. Then install EnPT itself:

   .. code-block:: bash

    $ conda install -c conda-forge enpt


This is the preferred method to install EnPT, as it always installs the most recent stable release and
automatically resolves all the dependencies.


Using pip (not recommended)
---------------------------

There is also a `pip`_ installer for EnPT. However, please note that EnPT depends on some
open source packages that may cause problems when installed with pip. Therefore, we strongly recommend
to resolve the following dependencies before the pip installer is run:

    * cachetools
    * cartopy
    * gdal
    * geopandas
    * glymur
    * h5py
    * numba
    * numpy
    * llvmlite
    * lxml
    * matplotlib
    * pandas
    * pygrib
    * pyhdf
    * pyproj >2.2.0
    * pyresample
    * pytables
    * scikit-image
    * scikit-learn
    * shapely


Then, the pip installer can be run by:

   .. code-block:: bash

    $ pip install enpt

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.



.. note::

    EnPT has been tested with Python 3.6+., i.e., should be fully compatible to all Python versions from 3.6 onwards.


.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
