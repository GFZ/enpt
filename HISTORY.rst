=======
History
=======

0.x.x (coming soon)
-------------------

New features:
* TBD


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
