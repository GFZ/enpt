=======
History
=======


0.11.6 (2020-04-06)
-------------------

* Updated .zenodo.json.


0.11.5 (2020-04-06)
-------------------

* Fixed EnPT logo in README.rst.
* Updated README.rst text.
* Pages now expire after 10 days instead of 30 days.
* Added .zenodo.json.


0.11.4 (2020-04-06)
-------------------

* Removed deprecated channels from environment_enpt.yml.
* Simplified test_enpt_install.
* Added SICOR to setup.py requirements.
* Updated installation instructions.


0.11.3 (2020-04-03)
-------------------

* Fixed broken badge4.
* Replaced logo relative link in README.rst with URL.


0.11.2 (2020-04-02)
-------------------

* Updated setup.py and MANIFEST.in to exclude tests and examples directories from PyPI upload.


0.11.1 (2020-04-02)
-------------------

* Fixed invalid syntax for multiple authors and email addresses in setup.py.


0.11.0 (2020-04-02)
-------------------

New features / improvements:

* Added parameter 'vswir_overlap_algorithm' that provides 4 different algorithms how to deal with the VNIR/SWIR overlap.
* Revised orthorecifier module.
* Updated badges in README.rst.
* Added a GUI test mode to EnPTController.
* Added keywords to setup.py.
* Added 'deploy_pypi' CI job.
* Revised setup.py for a proper PyPI upload.
* Removed installation of 'icu=58.*' from installation.rst.

Bug fixes:

* Fixed issue 45 "Band central wavelength positions of L2A product cannot be read by ENVI."


0.10.0 (2020-03-03)
-------------------

New features / improvements:

* Added source code repository link to table of contents of documentation.
* Updated license notes, copyright info, contributor guidelines and logos.
* Updated author info.
* Revised package short description.
* Added arosics to requirements.
* SensorMapGeometryTransformer is now imported from new library sensormapgeo.
* Updated dependencies and added pip to environment_enpt.yml.
* Boolean values are now correctly passed from the command line interface to EnPT.
* Added a tutorial to the docs.
* Some code improvements.
* Added output validation to AC.
* The parameter 'disable_progressbars' is now correctly passed to SICOR.
* Added tqdm exception to license file and license headers.
* Adapted code to the current EnMAP format.

Bug fixes:

* Fixed "Encoding error: 'ascii' codec can't decode byte 0xc3 in position 320: ordinal not in range(128)".
* Fixed unexpected title exception during 'make docs'.
* Fixed broken badge. Removed ssh links.
* Fixed UTF-8 error when running setup.py. Updated installation instructions.
* Fix for wrong input parameter data types at 'enable_keystone_correction' and 'enable_vnir_swir_coreg'.
* Fixed scheme error: 'scale_factor_boa_ref must be of integer type'.
* Fix for not validating the input data for enmap_image_gapfill


0.9.0 (2019-10-18)
------------------

New features / improvements:

* added functionality to transform between EnMAP VNIR aand SWIR sensor geometry
  (improves accuracy of atmospheric correction and solves reflectance spikes within the VNIR / SWIR spectral overlap)


0.8.0 (2019-10-15)
------------------

New features / improvements:

* Fixed issue 29 (static TLS)
* Set DLR test data as default test data
* Enhanced logging in orthorectifier module
* Enhanced AC results due to updated SICOR implementation
  (currently dependent from SICOR branch "master")
* Fixed loggers failing to deserialize
* GitLab Pages are now working properly (documentation hosting)
* Fixed issue 28 (cutoff effect of orthorectification results)
* Fixed dead documentation links
* Updated DLR test data and revised DN/radiance conversion
  (fixes negative radiance and BOA reflectance values / saturation)
* AOT value is now read from metadata and passed to SICOR
* Added validation of EnMAP root directory
* Added documentation including algorithm descriptions, installation instructions usage examples and auto-generated docs
* Added license texts


0.7.0 (2019-01-21)
------------------

New features / improvements:

* Added a lot of software tests
* Added output writer for EnMAP Level-2 data
* Added metadata class for EnMAP Level-2 data
* Revised dead pixel correction (now 40-50 times faster; added spatial interpolation)
* Added support for dead pixel correction based on 3D dead pixel maps
* Added orthorectification module
* Added support for 3D (band-wise) geometry layers
* Added 3D geolayer generation based on band-wise RPC coefficients.
* Updated L1B reader to match DLR L1B format
* Added subsets of official DLR test data
* Improved DEM processor (added overlap and geographic datum check)


0.6.0 (2018-12-13)
-------------------

New features:

* Updated test datasets (bugfix for wrong corner coordinates)
* Added dem in map geometry to test data
* Added spatial_transform module to transform between sensor and map geometry
* Added first version of dem_preprocessor module for pre-processing elevation data
* Added tests for new modules
* Added parameters 'path_dem' and 'average_elevation' to config parameters


0.5.0 (2018-06-13)
------------------

New features:

* Added algorithm to automatically append a second EnMAP image to the main image in order to fill the along-track gap
* Updated test data (updated metadata header file, now 2 EnMAP subset scenes)
* Updated metadata reader
* Save extended image


0.4.0 (2018-06-01)
------------------
New features:

* Implemented dead pixel corrector
* Implemented SICOR atmospheric correction


0.3.0 (??)
----------

New features:

* TBD


0.2.0 (2017-08-24)
------------------

New features:

* Structure draft for all modules.
* First implementation of image and metadata classes.
* path_generator module
* Implemented Reader for EnMAP Level-1B products.


0.1.0 (2017-05)
---------------

* Initial development started.
