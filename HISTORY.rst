=======
History
=======


0.18.6 (2021-12-20)
-------------------

* Fixed AttributeError in DEM_Processor related to new pyproj version.
* Fixed wrong FWHM metadata key in L2A data.
* Fixed slow conda environment creation.
* CI now uses Mambaforge.
* Migrated test calls from nosetests to pytest and implemented new test report (!60).


0.18.5 (2021-06-21)
-------------------

* Fixed incorrect gains in L2A output data.
* Band statistics in L2A metadata now exclude the image background.


0.18.4 (2021-06-21)
-------------------

* The 'enable_absolute_coreg' parameter is now automatically set to True if a reference image is provided via the CLI.
* Fixed that the 'progress' parameter was not correctly passed through to AROSICS.


0.18.3 (2021-06-18)
-------------------

* Fixed #79 (FileNotFoundError when using the official EnMAP test data downloaded from enmap.org).
* Removed required requirements_dev.txt which not needed anymore (also closes https://github.com/GFZ/enpt/pull/1).


0.18.2 (2021-06-17)
-------------------

* Updated installation instructions and GUI screenshots.
* Muted np.bool DeprecationWarning.


0.18.1 (2021-06-16)
-------------------

* Added method to EnPT_Controller to reveal Queue names to enpt_enmapboxapp
  (needed to fix missing log messages in EnPT GUI).


0.18.0 (2021-06-15)
-------------------

* Added config parameters to run EnPT in 3 AC modes: 'land', 'water', 'combined'.
* Added advanced atmospheric correction for water surfaces using ACwater/Polymer including tests and docs.


0.17.4 (2021-06-15)
-------------------

* Moved bin/enpt_cli.py to enpt/cli.py to fix missing argparser when installed from conda.


0.17.3 (2021-06-04)
-------------------

* Added documentation for pixel value 3 of land/water mask (#73).
* 'make lint' now directly prints errors instead of only logging them to logfiles.
* EnPTConfig and EnPT_Controller are now importable directly from the top level of EnPT.
* Added version pinnings to Python and QGIS.
* Replaced entry point 'enpt-cli.py' with 'enpt' which allows to call the argument parser with 'enpt <args>'.
* Don't assert existing file paths in config validation n case IS_ENPT_GUI_TEST==1.
* Default orthorectification algorithm is now 'gauss'.


0.17.2 (2021-03-04)
-------------------

* Added "if __name__ == '__main__'" blocks below test modules to ensure tests are also running properly on Windows.


0.17.1 (2021-02-23)
-------------------

* Fixed incorrectly pinned version of sicor.


0.17.0 (2021-02-22)
-------------------

* Moved sicor download from build_enpt_testsuite_image.sh to new before_script.sh
  and adjusted 'make gitlab_CI_docker' accordingly.
* Fixed wrong package name in environment_enpt.yml.
* Added sicor as conda environment to environment_enpt.yml which also makes before_script.sh obsolete.
* Coverage files are now completely deleted by running 'make clean-test'.
* Implemented new methods to transform raster arrays between sensor and map geometry.
* VNIR_SWIR_SensorGeometryTransformer now also allows to transform entire 3D data arrays.
* Added test_images_sensorgeo.py
* Removed deprecated config options 'enable_ice_retrieval' and 'auto_download_ecmwf'.
* EnPT is now compatible with the latest SICOR algorithm which fixes EnMAP/sicor#40
  ([EnMAP] SWIR parameter retrieval results are applied to wrong VNIR coordinates).


0.16.4 (2020-12-10)
-------------------

* Added URL checker CI job and fixed all dead URLs.


0.16.3 (2020-11-03)
-------------------

* Replaced deprecated 'source activate' by 'conda activate'.
* Updated installation instructions.
* Fixed DeprecationWarning (Using or importing the ABCs from 'collections' instead of from 'collections.abc' is
  deprecated since Python 3.3, and in 3.9 it will stop working)


0.16.2 (2020-10-12)
-------------------

* Added missing requirement 'scikit-image'.


0.16.1 (2020-10-12)
-------------------

* Added pygrib to pip installation instructions.
* Added test, doc and lint requirements as optional requirements to setup.py.


0.16.0 (2020-10-12)
-------------------

* Replaced deprecated function 'add_stylesheet' in conf.py.
* Use SPDX license identifier and set all files to GLP3+ to be consistent with license headers in the source files.
* EnPT is now on conda-forge! Updated the installation instructions accordingly. Added conda badge.
* Revised requirements and environment_enpt.yml.


0.15.0 (2020-09-21)
-------------------

* Added functionality to apply co-registration between an EnMAP image and a user-provided spatial reference dataset
  (still needs to be improved but can already be used). This includes: Spatial_Optimizer class, Test_Spatial_Optimizer
  class, updated config parameters, spatial reference test image.


0.14.1 (2020-09-01)
-------------------

* Updated installation instructions to include SICOR installation there.
* Some adaptions to recent changes in py_tools_ds and pyproj.


0.14.0 (2020-08-14)
-------------------

* Implemented new configuration option (parameter 'drop_bad_bands') to leave out certain water absorption bands from
  processing.
* Bugfix for allowing directories at filename config parameters.
* Total runtime is now printed.


0.13.2 (2020-08-11)
-------------------

* Fixed issue #37 (Empty FileNotFoundError in case no DEM is provided by the user). If no DEM is provided by the user,
  an average elevation (default=0) is used instead. Added corresponding warning and log messages.
* Fixed exception in case an average elevation is passed.
* Fixed a bug that caused 0 or 1 integers provided in the config parameters to be converted to booleans.
* User-provided file paths are now validated directly on startup.
* DEM_Processor instance can now be created from an average elevation (returns a flat DEM).
* RPC geolayer generators are now also compatible with an integer as input elevation.


0.13.1 (2020-08-07)
-------------------

* Updated the minimal version of sensormapgeo. This reduces the memory consumption during orthorectification when
  running in multiprocessing. Fixes issue #59 (High memory consumption of orthorectification causes deadlock or crash
  on smaller machines.)


0.13.0 (2020-05-18)
-------------------

* Renamed DEM for Arcachon test dataset.
* Fixed typo.
* Added config parameter 'output_format'.
* Implementated ENVI output format.
* Fixed log message.
* The output interleave is now also configurable via the parameter 'output_interleave'.
* Implemented 3 new config parameters: 'target_projection_type', 'target_epsg' and 'target_coord_grid'. This allows
  the user to choose between UTM and geographic L2A projection, to specify a custom L2A projection by providing an
  EPSG code or to specify a certain L2A coordinate grid.
* Added 'grid_res' and 'epsg' attributes to EnMAP_Metadata_L2A_MapGeo object.
* The L2A projection metadata is now correctly written to the XML file.


0.12.8 (2020-05-13)
-------------------

* L2A METADATA.XML now contains correct band characteristics, band statistics and information about the merged
  VNIR/SWIR cube (fixes issue #53).


0.12.7 (2020-05-12)
-------------------

* Revised RPC_3D_Geolayer_Generator such that the geolayer computation is now only executed for bands with differing
  RPC coefficients. Since we only have one set of RPC coefficients per detector by now, this speeds up the geolayer
  generation by the factor 20-50.
* Added Test_RPC_3D_Geolayer_Generator class.
* Added natsort to dependencies.
* Optimized RPC_3D_Geolayer_Generator algorithm to further reduce computation time. Its now 2-3x faster.
  Revised and added corresponding tests.


0.12.6 (2020-05-09)
-------------------

* DN2TOARadiance() now returns float32 radiance instead of float64.
* Bugfix in Orthorectifier._get_common_extent().
* Revised Geometry_Transformer and Geometry_Transformer_3D classes.
* Updated minimal version of sensormapgeo which makes the orthorectification much faster (factor 6-10)
  and fixes the deadlock within sensormapgeo.
* Bilinear orthorectification now uses 8 neighbours instead of 32 which makes it much faster.


0.12.5 (2020-05-04)
-------------------

* Dead-pixel correction is now called once more after AC to correct possible spectral spikes due to fringe effect.


0.12.4 (2020-05-04)
-------------------

* Revised computation of the common VNIR/SWIR extent within orthorectification (fixes issue #34). This computation now
  also respects deviations in per-band geolayers due to keystone or misregistrations.
* All pixels that have values in VNIR or SWIR only are not set to nodata in the L2A output (fixes issue #34).
* Nodata values of masks are now set.


0.12.3 (2020-04-21)
-------------------

* Fixed issue #50 (Quicklook images only contain noise).
* Fix for using the wrong bands for the SWIR quicklook image.


0.12.2 (2020-04-21)
-------------------

* L1B masks are now correctly written to the L2A output (fixes issue #48). However, the dead pixel map and the quality
  test flags are still missing.
* Silenced warning during closing of logging handler.


0.12.1 (2020-04-20)
-------------------

* Tests now use Arcachon test data as of 02/2020.
* Mask filenames are now correctly read from XML.
* Refactored filenames within metadata object to clean up the namespace.
* Disabled AC within tests.
* Converted type hints to Python 3.6 style.
* Dropped Python 3.5 support.
* Added filenames for masks to metadata.
* Added attribute 'epsg_ortho' to metadata.
* Revised _EnMAP_Image.generate_quicklook().
* __EnMAP_Image.paths is now correctly assigned.
* Split modules 'images' and 'metadata' into several sub-modules.
* Renamed image and metadata model modules for more clarity.
* Removed _EnMAP_Image properties 'mask_clouds_confidence', 'ac_errors' and 'ac_options'. Cleaned code duplicates.
* EnMAPL1Product_SensorGeo.transform_VNIR_to_SWIR_sensorgeo() now supports multiprocessing.
* Added mask attributes to sensor geometry image classes.
* Mask paths are now correctly set.
* Masks are now read from disk.
* Added subclasses EnMAP_VNIR_SensorGeo and EnMAP_SWIR_SensorGeo.
* Added functionality to set SWIR raster attributes with VNIR raster attributes tranformed to SWIR sensor geometry.
* The enmap_ImageL1 instance passed to SICOR now features a 'mask_water' attribute.
* Revised test_l1b_reader.py.
* Combined 'mask_water' and 'mask_land' attributes to 'mask_landwater'.
* Renamed metadata attribute 'filename_mask_deadpixel' to 'filename_deadpixelmap' for consistency.


0.12.0 (2020-04-09)
-------------------

* Added new L1B EnMAP test datasets for Arcachon, France (status 14.02.2020) + corresponding DEM.
* BSQ files now use Git LFS.


0.11.8 (2020-04-09)
-------------------

* Releases in the GitHub-Mirror-Repository are now created automatically
  (added create_release_from_gitlab_ci.sh and create_github_release CI job).
* Added GitHub issue template.


0.11.7 (2020-04-07)
-------------------

* Updated .zenodo.json.
* Added CITATION file.
* Added hint regarding citation to README.rst.


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
