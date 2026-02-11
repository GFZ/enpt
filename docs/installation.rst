.. _installation:

Installation
============

EnPT can be installed in two different ways:

- as a standalone package which can be used via the Python API or from the command line
- along with the EnMAP-Box (a QGIS plugin) which provides a GUI for EnPT.


Installing EnPT as a standalone package (backend code only)
***********************************************************


Using Miniforge (recommended)
-----------------------------

This is the preferred way to install EnPT. It is the fastest one and it always installs the most
recent stable release and automatically resolves all the dependencies.

1. Install Miniforge_.
2. Open a Miniforge command line prompt and proceed there (e.g., on Windows you can find it in the start menu).
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
* isofit =3.4.3
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



Optional: Get an API key to download Polymer auxiliary data
-----------------------------------------------------------

If you intend to use `ACwater`_/`Polymer`_ for atmospheric correction over water surfaces,
you need to get CDS API key to download auxiliary data at runtime.
For this, please register at the `CDS registration page`_ and store your `CDS API key`_
in your home directory. Further details are given in the `ACwater Polymer installation instructions`_.

  .. note::

    As of EnPT v1.1.4 and `ACwater`_ v0.4.0, the installation of the latest version of `Polymer`_ is
    fully integrated into the installation of EnPT. A separate installation of Polymer is no longer required.


Installing EnPT along with QGIS and the EnMAP-Box (backend + GUI)
*****************************************************************

If you want to use EnPT including the GUI_ in the EnMAP-Box_, it is highly recommended to install QGIS_,
the EnMAP-Box_ requirements, the EnPT backend code and the EnPT GUI_ into a single conda_ environment
within Miniforge_.

To do so, run the following command on a Miniforge_ conda_ command line:

.. code-block:: bash

  $ conda env create -n enpt_full -f https://git.gfz.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/raw/main/tests/gitlab_CI_docker/context/environment_enpt_full.yml

Then activate the newly created conda_ environment and start QGIS_:

.. code-block:: bash

  $ conda activate enpt_full
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

      $ conda env create -n enpt_full -f https://git.gfz.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/raw/main/tests/gitlab_CI_docker/context/environment_enpt_full_dev.yml


.. note::

    EnPT has been tested with Python 3.9+ on Linux, Windows and Mac OSX.


.. _Miniforge: https://github.com/conda-forge/miniforge
.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _conda: https://docs.conda.io
.. _ACwater: https://gitlab.awi.de/phytooptics/acwater/
.. _`ACwater Polymer installation instructions`: https://gitlab.awi.de/phytooptics/acwater/-/blob/master/docs/installation.rst
.. _Polymer: https://github.com/hygeos/polymer
.. _SICOR: https://git.gfz.de/EnMAP/sicor
.. _ISOFIT: https://github.com/isofit/isofit
.. _GUI: https://git.gfz.de/EnMAP/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp
.. _EnMAP-Box: https://github.com/EnMAP-Box/enmap-box
.. _QGIS: https://www.qgis.org
.. _CDS registration page: https://cds.climate.copernicus.eu/
.. _CDS API key: https://cds.climate.copernicus.eu/how-to-api
