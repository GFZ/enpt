.. _installation:

Installation
============

EnPT can be installed in two different ways:

- as a standalone package which can be used via the Python API or from the command line
- along with the EnMAP-Box (a QGIS plugin) which provides a GUI for EnPT.


Installing EnPT as a standalone package (backend code only)
***********************************************************


Using Mambaforge (recommended)
------------------------------

This is the preferred way to install EnPT. It is the fastest one and it always installs the most
recent stable release and automatically resolves all the dependencies.

1. Install Mambaforge_.
2. Open a Mambaforge command line prompt and proceed there (e.g., on Windows you can find it in the start menu).
3. Install enpt into a separate environment and activate it:

   .. code-block:: bash

    $ mamba create --name enpt enpt
    $ mamba activate enpt


Using Anaconda or Miniconda (slower)
------------------------------------

Using conda_ (latest version recommended), EnPT is installed as follows:

.. code-block:: bash

    $ conda create --name enpt -c conda-forge enpt
    $ conda activate enpt


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
* packaging
* pandas
* pygrib
* pyhdf
* pyproj >3.4.0
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



Optional: Install Polymer for advanced atmospheric correction over water surfaces
---------------------------------------------------------------------------------

For atmospheric correction of water surfaces using the Polymer algorithm instead of SICOR_ (which is mainly
designed for land surfaces), the additional package polymer_ is required. The Polymer atmospheric correction is made
available in EnPT by the ACwater_ package, a wrapper developed by AWI Bremerhaven, which is already contained in the
default EnPT installation.

1. To install the optional package polymer_, first activate the previously created enpt conda_ environment:

   .. code-block:: bash

    $ mamba activate enpt

2. Then register at the `HYGEOS support forum`_, download polymer_ from there (EnPT was tested with polymer v4.16.1,
   later versions may fail), unpack it and run the following commands from the unpacked root directory of polymer_:

   .. code-block:: bash

    $ make
    $ make auxdata_common
    $ make ancillary
    $ mkdir -p ANCILLARY/ERA5
    $ pip install -e .

  .. note::

    When using a conda_ environment on Linux or Mac OSX, the needed compilers to build polymer_
    should be already installed. On Windows, you need to install the `Microsoft build tools for visual studio`_
    including the C++ build tools, the latest versions of MSVCv142 - VS 2019 C++ x64/x86 build tools and Windows 10 SDK
    (see `here <https://wiki.python.org/moin/WindowsCompilers>`__ for details).
    However, polymer_ is currently *not Windows compatible* and will likely not run as expected.


Apart from that, you need to register at the `CDS registration page`_ and install a `CDS API key`_. This is required
to automatically download atmospheric AUX data at runtime, which are needed to run Polymer. Further details are
given in the `ACwater Polymer installation instructions`_.


Installing EnPT along with QGIS and the EnMAP-Box (backend + GUI)
*****************************************************************

If you want to use EnPT including the GUI_ in the EnMAP-Box_, it is highly recommended to install QGIS_,
the EnMAP-Box_ requirements, the EnPT backend code and the EnPT GUI_ into a single conda_ environment
within Mambaforge_.

To do so, run the following command on a Mambaforge_ conda_ command line:

.. code-block:: bash

  $ mamba env create -n enpt_full -f https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/raw/main/tests/gitlab_CI_docker/context/environment_enpt_full.yml

Then activate the newly created conda_ environment and start QGIS_:

.. code-block:: bash

  $ mamba activate enpt_full
  $ qgis

The EnMAP-Box_ QGIS_ plugin can then be installed via the QGIS_ Plugin manager and the EnPT GUI_ can be started
from within the EnMAP-Box_ as described
`here <https://enmap.git-pages.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp/doc/usage.html>`__.

If you want to use advanced atmospheric correction over water surfaces, please install the optional
requirement polymer_ into the enpt_full environment as described above.


.. hint::

    **Contributors** of the EnPT source code or plugins may install EnPT along with all packages needed for development
    with:

    .. code-block:: bash

      $ mamba env create -n enpt_full -f https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/raw/main/tests/gitlab_CI_docker/context/environment_enpt_full_dev.yml


.. note::

    EnPT has been tested with Python 3.8+ on Linux, Windows and Mac OSX.


.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://docs.conda.io
.. _ACwater: https://gitlab.awi.de/phytooptics/acwater/
.. _`ACwater Polymer installation instructions`: https://gitlab.awi.de/phytooptics/acwater/-/blob/master/docs/installation.rst
.. _HYGEOS support forum: https://forum.hygeos.com
.. _polymer: https://forum.hygeos.com
.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
.. _GUI: https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp
.. _EnMAP-Box: https://bitbucket.org/hu-geomatics/enmap-box
.. _QGIS: https://www.qgis.org
.. _CDS registration page: https://cds.climate.copernicus.eu/
.. _CDS API key: https://cds.climate.copernicus.eu/api-how-to
.. _Microsoft build tools for visual studio: https://visualstudio.microsoft.com/de/thank-you-downloading-visual-studio/?sku=BuildTools&rel=16
