.. _installation:

Installation
============

EnPT can be installed in two different ways:

- as a standalone package which can be used via the Python API or from the command line
- along with the EnMAP-Box (a QGIS plugin) which provides a GUI for EnPT.


Installing EnPT as a standalone package (backend code only)
***********************************************************


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
designed for land surfaces), the additional packages ACwater_ (a Polymer wrapper developed by AWI Bremerhaven)
and polymer_ are required.

1. Using a previously created enpt conda_ environment (as described above), first install some dependencies:

   .. code-block:: bash

    $ conda activate enpt
    $ conda install -c conda-forge cdsapi cython gdal netcdf4 pygrib pyhdf xarray
    $ pip install ecmwf-api-client

2. Then register at the `HYGEOS support forum`_, download polymer_ from there, unpack it and
   run the following commands from the unpacked root directory of polymer_:

   .. code-block:: bash

    $ make
    $ make auxdata_common
    $ make ancillary
    $ pip install -e .

  .. note::

    When using a conda_ environment on Linux or Mac OSX, the needed compilers to build polymer_
    should be already installed. On Windows, you need to install the `Microsoft build tools for visual studio`_
    including the C++ build tools, the latest versions of MSVCv142 - VS 2019 C++ x64/x86 build tools and Windows 10 SDK
    (see `here <https://wiki.python.org/moin/WindowsCompilers>`__ for details).
    However, polymer_ is currently *not Windows compatible* and will likely not run as expected.


Apart from that, you need to register at the `CDS registration page`_ and install a `CDS API key`_.
Further details are given `here <https://gitlab.awi.de/phytooptics/acwater/-/blob/master/docs/installation.rst>`__.

3. Finally install ACwater:

   .. code-block:: bash

    $ pip install git+https://gitlab.awi.de/phytooptics/acwater.git


Further details about the installation of ACwater can be found in the `ACwater Polymer installation instructions`_.


Installing EnPT along with QGIS and the EnMAP-Box (backend + GUI)
*****************************************************************

If you want to use EnPT including the GUI_ in the EnMAP-Box_, it is highly recommended to install QGIS_,
the EnMAP-Box_ requirements, the EnPT backend code and the EnPT GUI_ into a single conda_ environment.

To do so, run the following command on a conda_ command line:

   .. code-block:: bash

    $ conda env create -n enpt_full -f https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/raw/master/tests/gitlab_CI_docker/context/environment_enpt_full.yml

Then activate the newly created conda_ environment and start QGIS_:

   .. code-block:: bash

    $ conda activate enpt_full
    $ qgis

The EnMAP-Box_ QGIS_ plugin can then be installed via the QGIS_ Plugin manager and the EnPT GUI_ can be started
from within the EnMAP-Box_ as described
`here <https://enmap.git-pages.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp/doc/usage.html>`__.

If you want to use advanced atmospheric correction over water surfaces, please install the optional
requirement polymer_ as described above.




.. note::

    EnPT has been tested with Python 3.6+ on Linux, Windows and Mac OSX.


.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://conda.io/docs
.. _ACwater: https://gitlab.awi.de/phytooptics/acwater/
.. _`ACwater Polymer installation instructions`: https://gitlab.awi.de/phytooptics/acwater/-/blob/master/docs/installation.rst
.. _HYGEOS support forum: https://forum.hygeos.com
.. _polymer: https://forum.hygeos.com
.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
.. _GUI: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp
.. _EnMAP-Box: https://bitbucket.org/hu-geomatics/enmap-box
.. _QGIS: https://www.qgis.org
.. _CDS registration page: https://cds.climate.copernicus.eu/
.. _CDS API key: https://cds.climate.copernicus.eu/api-how-to
.. _Microsoft build tools for visual studio: https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16
