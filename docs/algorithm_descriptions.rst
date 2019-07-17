Algorithm descriptions
======================

EnMAP Level 1B data reader
**************************

The EnPT reader module allows you to read in EnMAP Level-1B data in the official format as provided by the DLR.
After extraction EnPT expects a folder with the following files:

    +-----------------------------------------------+-----------------+
    | Filename                                      | Description     |
    +===============================================+=================+
    |ENMAP*L1B*-HISTORY.XML                         |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-LOG.XML                             |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-METADATA.XML                        |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_PIXELMASK_SWIR.GEOTIFF           |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_PIXELMASK_VNIR.GEOTIFF           |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CIRRUS.GEOTIFF           |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLASSES.GEOTIFF          |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUD.GEOTIFF            |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_CLOUDSHADOW.GEOTIFF      |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_HAZE.GEOTIFF             |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_SNOW.GEOTIFF             |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_SWIR.GEOTIFF   |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_QUALITY_TESTFLAGS_VNIR.GEOTIFF   |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_SWIR.GEOTIFF                     |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-QL_VNIR.GEOTIFF                     |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_SWIR.GEOTIFF         |                 |
    +-----------------------------------------------+-----------------+
    |ENMAP*L1B*-SPECTRAL_IMAGE_VNIR.GEOTIFF         |                 |
    +-----------------------------------------------+-----------------+

EnPT reads the raster files lazily / on-demand to save memory during execution. Metadata read from the
`ENMAP*-METADATA.XML` file are stored in memory and used at multiple stages of the pre-processing pipeline.

The image data is directly transformed from digital numbers (DNs, as provided by the DLR) to top-of-atmosphere radiance
in mW/m2/sr/nm.


Dead pixel correction
*********************

The EnPT dead pixel correction uses the pixel masks provided by DLR and interpolates the EnMAP image data at
the indicated dead pixel positions. It supports two interpolation algorithms:

    1. spectral interpolation
        * Interpolates the data in the spectral domain.
        * Points outside the data range are extrapolated.
        * possible interpolation methods: `linear`, `nearest`, `zero`, `slinear`, `quadratic`, `cubic`, etc.
    2. spatial interpolation
        * Interpolates the data spatially.
        * Remaining missing data positions (e.g., outermost columns) are spectrally interpolated.
        * possible interpolation methods: `linear`, `bilinear`, `cubic`, `spline`

Import of an overlapping digital elevation model
************************************************

TBD

Atmospheric correction
**********************

EnPT uses `SICOR`_ (Sensor Independent Atmosperic Correction of optical Earth observation data from multi- and
hyperspectral instruments) for atmospheric correction. SICOR is a Python based open-source package developed at the
German Research Centre for Geosciences (GFZ) Potsdam. For details on the underlying algorithm, please refer to the
`documentation pages of SICOR`_.

.. _SICOR: https://gitext.gfz-potsdam.de/EnMAP/sicor
.. _`documentation pages of SICOR`: http://enmap.gitext.gfz-potsdam.de/sicor/doc/



Spatial Co-Registration
***********************

For the detection of spatial misregistrations with regard to a user-provided spatial reference EnPT builds on the
open-source Python package `AROSICS`_ (An Automated and Robust Open-Source Image Co-Registration Software for
Multi-Sensor Satellite Data). It has been developed at the German Research Centre for Geosciences (GFZ) Potsdam.
For detailed algorithm description and use cases refer to the corresponding (open-access) paper that can be found here:
`Scheffler D, Hollstein A, Diedrich H, Segl K, Hostert P. AROSICS: An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data. Remote Sensing. 2017; 9(7):676`_.

In EnPT, AROSICS is used to automacially compute thousands of tie points between a selected EnMAP band the
user-provided reference image. The computed shifts are later respected in the orthorectification step.

.. _AROSICS: https://gitext.gfz-potsdam.de/danschef/arosics
.. _`Scheffler D, Hollstein A, Diedrich H, Segl K, Hostert P. AROSICS: An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data. Remote Sensing. 2017; 9(7):676`: http://www.mdpi.com/2072-4292/9/7/676


Orthorectification
******************

EnMAP Level 1B data are provided in sensor geometry, i.e., the image data don't have map coordinates but only image
coordinates. For the geo-rectification of the data EnPT uses a set of Rational Polynomial Coefficients (RPCs) provided
for each band of the two EnMAP subsystems (VNIR and SWIR). Together with a user provided digital elevation model these
RPC coefficients enable a highly accurate assignment of map coordinates to each pixel of the EnMAP Level-1B images.
Resampling is done using a fast KDTree gaussian weighting neighbour approach implemented in the Python library
`pyresample`_ (`find the documentation here`_). The spatial shifts computed during the co-registration step are
respected here.

In this processing step, the EnMAP VNIR is merged with the SWIR subsystem and from now on stored in a single 3D array.

.. _pyresample: https://github.com/pytroll/pyresample
.. _find the documentation here: https://pyresample.readthedocs.io/en/latest/

EnMAP Level 2A data writer
**************************

The EnPT writer module writes the computed EnMAP Level-2A data to disk after finishing the processing pipeline. The
data format produced by EnPT is based on the official DLR Level-2A format. However, due to differences in the
underlying algorithms, EnPT also produces a slightly different Level-2A data format. The differences are summarized
below:

    +-----------------------------------------------+-----+---------+-------------+
    | Filename                                      | DLR | EnPT    | Description |
    +===============================================+=====+=========+=============+
    |ENMAP*L2A*-HISTORY.XML                         | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-LOG.XML                             | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-METADATA.XML                        | yes | yes     |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_PIXELMASK_SWIR.GEOTIFF           | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_PIXELMASK_VNIR.GEOTIFF           | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CIRRUS.GEOTIFF           | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLASSES.GEOTIFF          | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUD.GEOTIFF            | yes | yes     |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_CLOUDSHADOW.GEOTIFF      | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_HAZE.GEOTIFF             | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_QUALITY_SNOW.GEOTIFF             | yes | planned |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_SWIR.GEOTIFF                     | yes | yes     |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-QL_VNIR.GEOTIFF                     | yes | yes     |             |
    +-----------------------------------------------+-----+---------+-------------+
    |ENMAP*L2A*-SPECTRAL_IMAGE.GEOTIFF              | yes | yes     |             |
    +-----------------------------------------------+-----+---------+-------------+
