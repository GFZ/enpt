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



Optional: Install ACwater for advanced atmospheric correction over water surfaces
---------------------------------------------------------------------------------

For atmospheric correction of water surfaces using the Polymer algorithm instead of SICOR_ (which is mainly
designed for land surfaces), the additional package ACwater_ (a Polymer wrapper developed by AWI Bremerhaven)
is required.

1. Using a previously created enpt conda environment (as described above), first install some dependencies:

   .. code-block:: bash

    $ conda activate enpt
    $ conda install -c conda-forge cdsapi cython gdal netcdf4 pygrib pyhdf xarray
    $ pip install ecmwf-api-client

2. Then register at the `HYGEOS support forum`_, download polymer_ from there, unpack it and
   run :code:`pip install .` from the root directory of Polymer.

3. Finally install ACwater:

   .. code-block:: bash

    $ pip install git+https://gitlab.awi.de/phytooptics/acwater.git


Further details about the installation of ACwater can be found in the `ACwater Polymer installation instructions`_.



.. note::

    EnPT has been tested with Python 3.6+., i.e., should be fully compatible to all Python versions from 3.6 onwards.


.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
.. _ACwater: https://gitlab.awi.de/phytooptics/acwater/
.. _`ACwater Polymer installation instructions`: https://gitlab.awi.de/phytooptics/acwater/-/blob/master/docs/installation.rst
.. _HYGEOS support forum: https://forum.hygeos.com
.. _`polymer`: https://forum.hygeos.com
.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
