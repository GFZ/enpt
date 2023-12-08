.. _algorithm_description:

Algorithm descriptions
======================

EnMAP Level 1B data reader
**************************

The EnPT reader module allows you to read in EnMAP Level-1B data in the official format as provided by the EnMAP Ground
Segment. After data extraction EnPT expects a folder that includes the following files:

    +-----------------------------------------------+------------------------------------------------------+
    | Filename                                      | Description                                          |
    +===============================================+======================================================+
    |ENMAP*L1B*-HISTORY.XML                         | EnMAP ground segment history file for L1B processing |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-LOG.XML                             | EnMAP ground segment log file                        |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-METADATA.XML                        | L1B metadata file                                    |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_PIXELMASK_SWIR.TIF               | Defective pixel mask SWIR                            |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_PIXELMASK_VNIR.TIF               | Defective pixel mask VNIR                            |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_CIRRUS.TIF               | Quality cirrus                                       |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_CLASSES.TIF              | Quality classes                                      |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUD.TIF                | Quality cloud                                        |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUDSHADOW.TIF          | Quality cloud shadow                                 |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_HAZE.TIF                 | Quality haze                                         |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_SNOW.TIF                 | Quality snow                                         |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_SWIR.TIF       | Quality test flags SWIR                              |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_VNIR.TIF       | Quality test flags VNIR                              |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_SWIR.TIF                         | SWIR quicklook image                                 |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-QL_VNIR.TIF                         | VNIR quicklook image                                 |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_SWIR.TIF             | EnMAP L1B SWIR spectral image                        |
    +-----------------------------------------------+------------------------------------------------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_VNIR.TIF             | EnMAP L1B VNIR spectral image                        |
    +-----------------------------------------------+------------------------------------------------------+

EnPT reads the raster files lazily / on-demand to save memory during execution. Metadata read from the
`ENMAP*-METADATA.XML` file are stored in memory and used at multiple stages of the pre-processing pipeline.

**Conversion to top-of-atmosphere radiance**

The image data is directly transformed from digital numbers (DNs, as provided by the DLR) to top-of-atmosphere radiance
in mW/m2/sr/nm.

**Filling the gap between adjacent EnMAP images**

There is a gap of around 20 pixels in along-track-direction between the VNIR and SWIR image of the official Level-1B
product and up to 1 pixel in across-track direction. Since the Level-2A processing in EnPT requires overlapping data,
the data is cut to the overlap area, which in turn leads to blank spaces between two adjacent Level-2A tiles. To create
overlapping Level-2A tiles in along-track direction for a further image mosaicking, EnPT offers the opportunity to
consider information from neighboring Level-1B image tiles which is then appended to the actual EnMAP image tile to be
processed.




Dead pixel correction
*********************

The EnPT dead pixel correction uses the pixel masks provided by the official Level-1B product and interpolates the
EnMAP image data at the indicated dead pixel positions. It supports two interpolation algorithms:

    1. spectral interpolation (default)
        * Interpolates the data in the spectral domain.
        * Points outside the data range are extrapolated.
        * possible interpolation methods: `linear`, `nearest`, `zero`, `slinear`, `quadratic`, `cubic`, etc.
    2. spatial interpolation
        * Interpolates the data spatially.
        * Remaining missing data positions (e.g., outermost columns) are spectrally interpolated.
        * possible interpolation methods: `linear`, `bilinear`, `cubic`, `spline`




Atmospheric correction
**********************

.. image:: https://git.gfz-potsdam.de/EnMAP/sicor/raw/main/docs/images/sicor_logo_lr.png
   :target: https://git.gfz-potsdam.de/EnMAP/sicor
   :width: 150px
   :alt: SICOR Logo

EnPT uses `SICOR`_ (Sensor Independent Atmosperic Correction of optical Earth observation data from multi- and
hyperspectral instruments) for atmospheric correction, i.e., for the conversion of TOA-(top-of-atmosphere) radiance
to BOA- (bottom-of-atmosphere / surface) reflectance. SICOR is a Python based open-source package developed at the
German Research Centre for Geosciences (GFZ) Potsdam. For details on the underlying algorithm, please refer to the
`documentation pages of SICOR`_.

Optionally, EnPT uses the Polymer_ algorithm for atmospheric correction over water (`Steinmetz et al. (2011)`_).
Polymer_ is a spectral matching algorithm in which atmospheric and oceanic signals are obtained simultaneously using
the fully available visible spectrum. The algorithm was developed by Hygeos (https://www.hygeos.com/); it is available
as a Python package and it has been largely applied to ocean colour sensors. The Polymer algorithm was integrated into
EnPT using the wrapper module ACwater_ developed at AWI Bremerhaven in cooperation with GFZ. In addition, Polymer_
was further adapted to process EnMAP L1B satellite data. For details on the underlying Polymer_ algorithm, please
refer to `Steinmetz et al. (2011)`_ and `Soppa et al. (2021)`_.

The following `ACwater`_/`Polymer`_ outputs can be included into the EnMAP Level 2A product generated by EnPT (added
in version 0.19.0):

- normalized water leaving reflectance
- chlorophyll-a concentration (logchl)
- fb coefficient that scales the backscattering coefficient of particles (logfb)
- reflectance of the sun glint (Rgli)
- TOA reflectance at 865 nm corrected for Rayleigh scattering (Rnir)
- quality flags (bitmask).

.. _Polymer: https://www.hygeos.com/polymer
.. _ACwater: https://gitlab.awi.de/phytooptics/acwater
.. _`Steinmetz et al. (2011)`: https://doi.org/10.1364/OE.19.009783
.. _`Soppa et al. (2021)`: https://doi.org/10.3390/s21124125

Spatial Co-Registration
***********************

.. image:: https://git.gfz-potsdam.de/danschef/arosics/raw/main/docs/images/arosics_logo.png
   :target: https://git.gfz-potsdam.de/danschef/arosics
   :width: 150px
   :alt: AROSICS Logo

EnPT enables the geometric adaptation of EnMAP data to a user-provided image scene (e.g. Sentinel-2). Spatial
misregistrations are detected using the open-source Python package `AROSICS`_ (An Automated and Robust Open-Source
Image Co-Registration Software for Multi-Sensor Satellite Data). It has been developed at the German Research Centre
for Geosciences (GFZ) Potsdam. For detailed algorithm description and use cases refer to the corresponding
(open-access) paper that can be found here:
`Scheffler D, Hollstein A, Diedrich H, Segl K, Hostert P. AROSICS: An Automated and Robust Open-Source Image
Co-Registration Software for Multi-Sensor Satellite Data. Remote Sensing. 2017; 9(7):676`__.

In EnPT, AROSICS is used to automatically compute thousands of tie points between a selected EnMAP band the
user-provided reference image. The computed shifts are considered in the orthorectification step.

__ https://www.mdpi.com/2072-4292/9/7/676



.. VNIR/SWIR coregistration estimation???
.. Keystone estimation???


Orthorectification
******************

EnMAP Level 1B data are provided in sensor geometry, i.e., the image data don't have map coordinates but only image
coordinates. For the ortho-rectification of the data EnPT uses a set of Rational Polynomial Coefficients (RPCs) provided
for each band of the two EnMAP subsystems (VNIR and SWIR). Together with a user provided digital elevation model these
RPC coefficients enable a highly accurate assignment of map coordinates to each pixel of the EnMAP Level-1B images.
The RPC coefficients already include the official information about detector coregistration and keystone. This way
image map coordinates are calculated internally for each pixel and band considering the spatial misregistrations
estimated by AROSICS on demand. Resampling is done using a fast KDTree gaussian weighting neighbour approach
implemented in the Python library
`pyresample`_ (find the documentation `here <https://pyresample.readthedocs.io/en/latest/>`__).

In this processing step, the EnMAP VNIR is merged with the SWIR subsystem and from now on stored in a single 3D array.




EnMAP Level 2A data writer
**************************

The EnPT writer module writes the computed EnMAP Level-2A data (orthorectified bottom-of-atmosphere reflectance for
land surfaces or normalized water-leaving reflectance for water surfaces if the atmospheric correction runs in `water`
or `combined` mode) to disk after finishing the processing pipeline. The data format produced by EnPT is based on the
official Level-2A format of the ground segment. However, due to differences in the underlying algorithms, EnPT also
produces a slightly different Level-2A data format. The current differences are summarized below:

    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    | Filename                                      | official L2A format | EnPT     | Description                                                                                       |
    +===============================================+=====================+==========+===================================================================================================+
    |ENMAP*L2A*.log                                 |         no          | yes      | EnPT log file                                                                                     |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-HISTORY.XML                         |         yes         | no       | EnMAP ground segment history file for L2A processing                                              |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-LOG.XML                             |         yes         | no       | EnMAP ground segment log file                                                                     |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-METADATA.XML                        |         yes         | yes      | L2A metadata file                                                                                 |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_PIXELMASK.TIF                    |         yes         | planned  | Defective pixel mask                                                                              |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_CIRRUS.TIF               |         yes         | yes      | Quality cirrus                                                                                    |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_CLASSES.TIF              |         yes         | yes      | Quality classes                                                                                   |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUD.TIF                |         yes         | yes      | Quality cloud                                                                                     |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUDSHADOW.TIF          |         yes         | yes      | Quality cloud shadow                                                                              |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_HAZE.TIF                 |         yes         | yes      | Quality haze                                                                                      |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_SNOW.TIF                 |         yes         | yes      | Quality snow                                                                                      |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_QUALITY_TESTFLAGS.TIF            |         yes         | no       | Quality test flags                                                                                |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_SWIR.TIF                         |         yes         | yes      | SWIR quicklook image                                                                              |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-QL_VNIR.TIF                         |         yes         | yes      | VNIR quicklook image                                                                              |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-SPECTRAL_IMAGE.TIF                  |         yes         | yes      | EnMAP L2A bottom-of-atmosphere reflectance (land) or normalized water leaving reflectance (water) |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-ACOUT_POLYMER_*RNIR.TIF             |         no          | optional | TOA reflectance at 863 nm corrected for Rayleigh scattering                                       |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-ACOUT_POLYMER_*RGLI.TIF             |         no          | optional | Reflectance of the sun glint predicted from ECMWF wind speed                                      |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-ACOUT_POLYMER_*LOGCHL.TIF           |         no          | optional | Chlorophyll-a concentration (mg/m3, in 10-based logarithm)                                        |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-ACOUT_POLYMER_*LOGFB.TIF            |         no          | optional | Particle scattering factor fb in `Park & Ruddick (2005)`_ (in 10-based logarithm)                 |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+
    |ENMAP*L2A*-ACOUT_POLYMER_*BITMASK.TIF          |         no          | optional | Polymer quality flags (more information below)                                                    |
    +-----------------------------------------------+---------------------+----------+---------------------------------------------------------------------------------------------------+

The **Polymer quality flags bitmask** represents a bit-encoded product with the following flag values:

    +--------------------+-------------+--------------------------------------------+
    | Flag name          | Flag value  | Description                                |
    +====================+=============+============================================+
    | LAND               | 1           | Land mask                                  |
    +--------------------+-------------+--------------------------------------------+
    | CLOUD_BASE         | 2           | Polymer's basic cloud mask                 |
    +--------------------+-------------+--------------------------------------------+
    | L1_INVALID         | 4           | Invalid level1 pixel                       |
    +--------------------+-------------+--------------------------------------------+
    | NEGATIVE_BB        | 8           | (deprecated flag)                          |
    +--------------------+-------------+--------------------------------------------+
    | OUT_OF_BOUNDS      | 16          | Retrieved marine parameters are outside    |
    |                    |             | valid bounds                               |
    +--------------------+-------------+--------------------------------------------+
    | EXCEPTION          | 32          | A processing error was encountered         |
    +--------------------+-------------+--------------------------------------------+
    | THICK_AEROSOL      | 64          | Thick aerosol flag                         |
    +--------------------+-------------+--------------------------------------------+
    | HIGH_AIR_MASS      | 128         | Air mass exceeds 5                         |
    +--------------------+-------------+--------------------------------------------+
    | EXTERNAL_MASK      | 512         | Pixel was masked using external mask       |
    +--------------------+-------------+--------------------------------------------+
    | CASE2              | 1024        | Pixel was processed in "case2" mode        |
    +--------------------+-------------+--------------------------------------------+
    | INCONSISTENCY      | 2048        | Inconsistent result was detected           |
    |                    |             | (atmospheric reflectance out of bounds)    |
    +--------------------+-------------+--------------------------------------------+
    | ANOMALY_RWMOD_BLUE | 4096        | Excessive difference was found at 412nm    |
    |                    |             | between Rw and Rwmod                       |
    +--------------------+-------------+--------------------------------------------+

Value 0 represents water (all fine, no flags), value -9999 represents no-data.


.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
.. _`documentation pages of SICOR`: https://enmap.git-pages.gfz-potsdam.de/sicor/doc/
.. _AROSICS: https://git.gfz-potsdam.de/danschef/arosics
.. _pyresample: https://github.com/pytroll/pyresample
.. _`Park & Ruddick (2005)`: https://opg.optica.org/ao/abstract.cfm?uri=ao-44-7-1236