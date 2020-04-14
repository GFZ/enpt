# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz-potsdam.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz-potsdam.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz-potsdam.de)
#
# This software was developed within the context of the EnMAP project supported
# by the DLR Space Administration with funds of the German Federal Ministry of
# Economic Affairs and Energy (on the basis of a decision by the German Bundestag:
# 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version. Please note the following exception: `EnPT` depends on tqdm, which
# is distributed under the Mozilla Public Licence (MPL) v2.0 except for the files
# "tqdm/_tqdm.py", "setup.py", "README.rst", "MANIFEST.in" and ".gitignore".
# Details can be found here: https://github.com/tqdm/tqdm/blob/master/LICENCE.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""EnPT EnMAP object base classes."""

from types import SimpleNamespace
from typing import Tuple, Optional  # noqa: F401
import numpy as np

# noinspection PyPackageRequirements
from skimage import exposure  # contained in package requirements as scikit-image

from geoarray import GeoArray, CloudMask

from ...model.metadata import EnMAP_Metadata_L2A_MapGeo  # noqa: F401  # only used for type hint

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class _EnMAP_Image(object):
    """EnPT base class for all kinds of EnMAP images.

    NOTE:
        - Basic functionality that all EnMAP image objects have in common is to be implemented here.
        - All EnMAP image classes should (directly or indirectly) inherit from _EnMAP_image.

    Attributes:
        - to be listed here. Check help(_EnMAP_image) in the meanwhile!

    """

    def __init__(self):
        """Load and hold data needed for processing EnMAP Level-1B data to Level-2A.

        Intended usage:
            - only to be used as a base class to be inherited from!
            - Example:
                class EnMAP_Image_Subclass(_EnMAP_Image):
                    def __init__(self):
                        super(EnMAP_Image_Subclass, self).__init__()

        """
        # add EnMAP specific attributes
        self.paths = SimpleNamespace()

        # protected attributes
        self._data = None
        self._mask_nodata = None
        self._mask_clouds = None
        self._mask_clouds_confidence = None
        self._dem = None
        self._deadpixelmap = None
        self._ac_options = {}
        self._ac_errors = None
        self._subset = None  # FIXME(Stephane) how is _subset to be set?

        # defaults
        self.entity_ID = ''
        self.basename = ''

    @property
    def data(self) -> GeoArray:
        """Return the actual EnMAP image data.

        Bundled with all the corresponding metadata.

        Attributes and functions (most important; for a full list check help(self.data)!):
            - ALL attributes of numpy.ndarray!
            - is_inmem(bool):
                True if the image data are completely loaded into memory; False if GeoArray only holds a link to a file
                on disk.
            - arr: np.ndarray holding the pixel values (if is_mem is True)
            - rows(int)
            - cols(int)
            - bands(int)
            - shape(tuple)
            - gt(list):  GDAL geotransform: contains the geocoding
            - prj(str): WKT projection string
            - show(*args, **kwargs):  plot the image
            - show_map(*args, **kwargs):  plot a map of the image (based on Basemap library)
            - reproject_to_new_grid(*args, **kwargs)

        Usage (there will soon be detailed instructions on usage at https://gitext.gfz-potsdam.de/danschef/geoarray):

            - Use self.data like a normal numpy.ndarray!
                - NOTE: Operators like *, /, + , - will soon be implemented. In the meanwhile use:
                    result = self.data[:] *10

            - How to set self.data?
                - Link an image file to self.data -> all raster data is read into memory ON DEMAND:
                    self.data = '/path/to/image.tif'  # sets self.data to GeoArray('/path/to/image.tif')

                - Link a numpy.ndarray instance with self.data (remaining attributes like geocoding, projection, etc.
                    are copied from the previous self.data attribute.
                    self.data = numpy.array([[1,2,3],[4,5,6]])

                - Set self.data to an existing instance of GeoArray
                    (allows to set specific attributes of GeoArray by yourself)
                    self.data = GeoArray('/path/to/image.tif', geotransform=[...], projection='WKTString')

            - Delete self.data:
                del self.data or 'self.data = None'

        :return:    instance of geoarray.GeoArray
        """
        if self._subset is None:
            return self._data

        return GeoArray(self._data[self._subset[0]:self._subset[1], self._subset[2]:self._subset[3], :],
                        geotransform=self._data.gt, projection=self._data.prj)

    @data.setter
    def data(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            # TODO this must be able to handle subset inputs in tiled processing
            if isinstance(geoArr_initArgs[0], np.ndarray):
                self._data = GeoArray(geoArr_initArgs[0], geotransform=self._data.gt, projection=self._data.prj)
            else:
                self._data = GeoArray(*geoArr_initArgs)
        else:
            del self.data

    @data.deleter
    def data(self):
        self._data = None

    @property
    def mask_clouds(self) -> GeoArray:
        """Return the cloud mask.

        Bundled with all the corresponding metadata.

        For usage instructions and a list of attributes refer to help(self.data).
        self.mask_clouds works in the same way.

        :return: instance of geoarray.CloudMask
        """
        return self._mask_clouds

    @mask_clouds.setter
    def mask_clouds(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            cm = CloudMask(*geoArr_initArgs)
            if cm.shape[:2] != self.data.shape[:2]:
                raise ValueError("The 'mask_clouds' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(cm.shape))
            cm.nodata = 0
            cm.gt = self.data.gt
            cm.prj = self.data.prj
            self._mask_clouds = cm
        else:
            del self.mask_clouds

    @mask_clouds.deleter
    def mask_clouds(self):
        self._mask_clouds = None

    @property
    def mask_clouds_confidence(self) -> GeoArray:
        """Return pixelwise information on the cloud mask confidence.

        Bundled with all the corresponding metadata.

        For usage instructions and a list of attributes refer to help(self.data).
        self.mask_clouds_confidence works in the same way.

        :return: instance of geoarray.GeoArray
        """
        return self._mask_clouds_confidence

    @mask_clouds_confidence.setter
    def mask_clouds_confidence(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            cnfArr = GeoArray(*geoArr_initArgs)

            if not cnfArr.shape == self.data.shape[:2]:
                raise ValueError("The 'mask_clouds_confidence' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(cnfArr.shape))

            if cnfArr._nodata is None:
                cnfArr.nodata = 0  # DEF_D.get_outFillZeroSaturated(cnfArr.dtype)[0] # TODO
            cnfArr.gt = self.data.gt
            cnfArr.prj = self.data.prj

            self._mask_clouds_confidence = cnfArr
        else:
            del self.mask_clouds_confidence

    @mask_clouds_confidence.deleter
    def mask_clouds_confidence(self):
        self._mask_clouds_confidence = None

    @property
    def dem(self) -> GeoArray:
        """Return a DEM in the exact dimension and pixel grid of self.arr.

        :return: geoarray.GeoArray
        """
        if self._dem is None:
            raise NotImplementedError('An automatic DEM getter is not yet implemented.')
        return self._dem

    @dem.setter
    def dem(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dem = GeoArray(*geoArr_initArgs)

            if not dem.shape == self.data.shape[:2]:
                raise ValueError("The 'dem' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(dem.shape))
            self._dem = dem
            self._dem.nodata = 0  # FIXME
            self._dem.gt = self.data.gt
            self._dem.prj = self.data.prj
        else:
            del self.dem

    @dem.deleter
    def dem(self):
        self._dem = None

    @property
    def ac_errors(self) -> GeoArray:
        """Return error information calculated by the atmospheric correction.

        :return: geoarray.GeoArray
        """
        return self._ac_errors  # FIXME should give a warning if None

    @ac_errors.setter
    def ac_errors(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            errArr = GeoArray(*geoArr_initArgs)

            if errArr.shape != self.data.shape:
                raise ValueError("The 'ac_errors' GeoArray can only be instanced with an array of "
                                 "the same dimensions like _EnMAP_Image.arr. Got %s." % str(errArr.shape))

            if errArr._nodata is None:
                errArr.nodata = 0  # DEF_D.get_outFillZeroSaturated(errArr.dtype)[0] # TODO
            errArr.gt = self.data.gt
            errArr.prj = self.data.prj
            # errArr.bandnames = self.LBA2bandnames(self.LayerBandsAssignment)

            self._ac_errors = errArr
        else:
            del self.ac_errors

    @ac_errors.deleter
    def ac_errors(self):
        self._ac_errors = None

    @property
    def deadpixelmap(self) -> GeoArray:
        """Return the dead pixel map.

        Bundled with all the corresponding metadata. Dimensions: (bands x columns).

        For usage instructions and a list of attributes refer to help(self.data).
        self.mask_clouds_confidence works in the same way.

        :return: instance of geoarray.GeoArray
        """
        if self._deadpixelmap is not None:
            self._deadpixelmap.arr = self._deadpixelmap[:].astype(np.bool)  # ensure boolean map

        return self._deadpixelmap

    @deadpixelmap.setter
    def deadpixelmap(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dpm = GeoArray(*geoArr_initArgs)

            if dpm.ndim == 3 and dpm.shape != self.data.shape:
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with a 3D array with the same size "
                                 "like _EnMAP_Image.arr, i.e.: %s "
                                 "Received %s." % (str(self.data.shape), str(dpm.shape)))
            elif dpm.ndim == 2 and dpm.shape != (self.data.bands, self.data.cols):
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with an array with the size "
                                 "'bands x columns' of the GeoArray _EnMAP_Image.arr. "
                                 "Received %s. Expected %s" % (str(dpm.shape), str((self.data.bands, self.data.cols))))

            self._deadpixelmap = dpm
        else:
            del self.deadpixelmap

    @deadpixelmap.deleter
    def deadpixelmap(self):
        self._deadpixelmap = None

    def generate_quicklook(self, bands2use: Tuple[int, int, int]) -> GeoArray:
        """
        Generate image quicklook and save into a file

        :param bands2use:   (red, green, blue) band indices of self.data to be used for quicklook image, e.g., (3, 2, 1)
        :return: GeoArray
        """
        red, green, blue = [self.data[:, :, bandidx] for bandidx in bands2use]

        def rescale(bandarray):
            p2 = np.percentile(bandarray, 2)
            p98 = np.percentile(bandarray, 98)
            return exposure.rescale_intensity(bandarray, (p2, p98))

        pix = np.dstack((rescale(red), rescale(green), rescale(blue)))
        pix = np.uint8(pix * 255)

        return GeoArray(pix)
