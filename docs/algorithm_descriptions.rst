.. _algorithm_description:

Algorithm descriptions
======================

EnMAP Level 1B data reader
**************************

The EnPT reader module allows you to read in EnMAP Level-1B data in the official format as provided by the EnMAP Ground
Segment. After data extraction EnPT expects a folder that includes the following files:

    +-----------------------------------------------+-----------------+
    | Filename                                      | Description     |
    +===============================================+=================+
    |ENMAP*L1B*-HISTORY.XML                         |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-LOG.XML                             |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-METADATA.XML                        |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_PIXELMASK_SWIR.TIF               |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_PIXELMASK_VNIR.TIF               |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CIRRUS.TIF               |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLASSES.TIF              |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUD.TIF                |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUDSHADOW.TIF          |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_HAZE.TIF                 |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_SNOW.TIF                 |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_SWIR.TIF       |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_VNIR.TIF       |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_SWIR.TIF                         |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_VNIR.TIF                         |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_SWIR.TIF             |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_VNIR.TIF             |                 |
    +-----------------------------------------------+-----------------+

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

.. image:: https://git.gfz-potsdam.de/EnMAP/sicor/raw/master/docs/images/sicor_logo_lr.png
   :target: https://git.gfz-potsdam.de/EnMAP/sicor
   :width: 150px
   :alt: SICOR Logo

EnPT uses `SICOR`_ (Sensor Independent Atmosperic Correction of optical Earth observation data from multi- and
hyperspectral instruments) for atmospheric correction, i.e., for the conversion of TOA-(top-of-atmosphere) radiance
to BOA- (bottom-of-atmosphere / surface) reflectance. SICOR is a Python based open-source package developed at the
German Research Centre for Geosciences (GFZ) Potsdam. For details on the underlying algorithm, please refer to the
`documentation pages of SICOR`_.

Optionally, EnPT retrieves water reflectance above the surface using `ACwater Polymer`_.
ACwater Polymer is a "wrapper" package (developed at the Alfred-Wegener-Institute, Bremerhaven)
for the `Polymer`_ atmospheric correction (AC) algorithm (developed by Hygeos, Inc).
Polymer AC is based on an optimization technique that considers atmospheric and oceanic signals to retrieve
normalized spectral reflectance above water. For details regarding the Polymer algorithm,
users are referred to `Steinmetz F, Deschamps P-Y, Ramon R., Opt. Express. 2011; 19`__.

__ https://doi.org/10.1364/OE.19.009783

.. _`ACwater Polymer`: https://gitlab.awi.de/phytooptics/acwater
.. _Polymer: https://www.hygeos.com/polymer

Spatial Co-Registration
***********************

.. image:: https://git.gfz-potsdam.de/danschef/arosics/raw/master/docs/images/arosics_logo.png
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

__ http://www.mdpi.com/2072-4292/9/7/676


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

The EnPT writer module writes the computed EnMAP Level-2A data to disk after finishing the processing pipeline. The
data format produced by EnPT is based on the official Level-2A format. However, due to differences in the
underlying algorithms, EnPT also produces a slightly different Level-2A data format. The current differences are
summarized below:

    +-----------------------------------------------+---------------------+---------+-------------+
    | Filename                                      | official L2A format | EnPT    | Description |
    +===============================================+=====================+=========+=============+
    |ENMAP*L2A*-HISTORY.XML                         |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-LOG.XML                             |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-METADATA.XML                        |         yes         | yes     |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_PIXELMASK_SWIR.GEOTIFF           |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_PIXELMASK_VNIR.GEOTIFF           |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CIRRUS.GEOTIFF           |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLASSES.GEOTIFF          |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUD.GEOTIFF            |         yes         | yes     |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUDSHADOW.GEOTIFF      |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_HAZE.GEOTIFF             |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_SNOW.GEOTIFF             |         yes         | planned |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_SWIR.GEOTIFF                     |         yes         | yes     |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-QL_VNIR.GEOTIFF                     |         yes         | yes     |             |
    +-----------------------------------------------+---------------------+---------+-------------+
    |ENMAP*L2A*-SPECTRAL_IMAGE.GEOTIFF              |         yes         | yes     |             |
    +-----------------------------------------------+---------------------+---------+-------------+


.. _SICOR: https://git.gfz-potsdam.de/EnMAP/sicor
.. _`documentation pages of SICOR`: https://enmap.git-pages.gfz-potsdam.de/sicor/doc/
.. _AROSICS: https://git.gfz-potsdam.de/danschef/arosics
.. _pyresample: https://github.com/pytroll/pyresample