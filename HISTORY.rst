=======
History
=======

0.x.x (coming soon)
-------------------

New features:
* TBD


0.7.0 (coming soon)
-------------------

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
