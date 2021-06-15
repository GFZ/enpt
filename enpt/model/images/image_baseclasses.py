# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
from typing import List, Optional  # noqa: F401
import numpy as np

# noinspection PyPackageRequirements
from skimage import exposure  # contained in package requirements as scikit-image

from geoarray import GeoArray

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

        # private attributes
        self._data = None
        self._mask_landwater = None
        self._mask_clouds = None
        self._mask_cloudshadow = None
        self._mask_haze = None
        self._mask_snow = None
        self._mask_cirrus = None
        self._dem = None
        self._deadpixelmap = None
        self._subset = None  # FIXME how is _subset to be set?

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
            - show():  plot the image
            - show_map():  plot a map of the image (based on cartopy library)
            - reproject_to_new_grid()

        Usage (there will soon be detailed instructions on usage at https://git.gfz-potsdam.de/danschef/geoarray):

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
    def mask_landwater(self) -> GeoArray:
        """Return the land/water mask.

        pixel values:
        - 0: background within scene dimensions, e.g. due to missing values/errors
        - 1: no water
        - 2: water
        - 3: background outside the scene dimensions (artifact from resampling between map and sensor geometry)

        :return: geoarray.GeoArray
        """
        return self._mask_landwater

    @mask_landwater.setter
    def mask_landwater(self, *geoArr_initArgs):
        self._mask_landwater = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_landwater', nodataVal=0)

    @mask_landwater.deleter
    def mask_landwater(self):
        self._mask_landwater = None

    @property
    def mask_clouds(self) -> GeoArray:
        """Return the cloud mask.

        :return: geoarray.GeoArray
        """
        return self._mask_clouds

    @mask_clouds.setter
    def mask_clouds(self, *geoArr_initArgs):
        self._mask_clouds = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_clouds', nodataVal=0)

    @mask_clouds.deleter
    def mask_clouds(self):
        self._mask_clouds = None

    @property
    def mask_cloudshadow(self) -> GeoArray:
        """Return the cloud shadow mask (0=no cloud shadow, 1=cloud shadow)..

        :return: geoarray.GeoArray
        """
        return self._mask_cloudshadow

    @mask_cloudshadow.setter
    def mask_cloudshadow(self, *geoArr_initArgs):
        self._mask_cloudshadow = \
            self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_cloudshadow', nodataVal=0)

    @mask_cloudshadow.deleter
    def mask_cloudshadow(self):
        self._mask_cloudshadow = None

    @property
    def mask_haze(self) -> GeoArray:
        """Return the haze mask (0=no haze, 1=haze)..

        :return: geoarray.GeoArray
        """
        return self._mask_haze

    @mask_haze.setter
    def mask_haze(self, *geoArr_initArgs):
        self._mask_haze = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_haze', nodataVal=0)

    @mask_haze.deleter
    def mask_haze(self):
        self._mask_haze = None

    @property
    def mask_snow(self) -> GeoArray:
        """Return the snow mask (0=no snow, 1=snow)..

        :return: geoarray.GeoArray
        """
        return self._mask_snow

    @mask_snow.setter
    def mask_snow(self, *geoArr_initArgs):
        self._mask_snow = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_snow', nodataVal=0)

    @mask_snow.deleter
    def mask_snow(self):
        self._mask_snow = None

    @property
    def mask_cirrus(self) -> GeoArray:
        """Return the cirrus mask (0=none, 1=thin, 2=medium, 3=thick)..

        :return: geoarray.GeoArray
        """
        return self._mask_cirrus

    @mask_cirrus.setter
    def mask_cirrus(self, *geoArr_initArgs):
        self._mask_cirrus = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_cirrus', nodataVal=0)

    @mask_cirrus.deleter
    def mask_cirrus(self):
        self._mask_cirrus = None

    @property
    def dem(self) -> GeoArray:
        """Return a DEM in the exact dimension and pixel grid of self.data.

        :return: geoarray.GeoArray
        """
        if self._dem is None:
            raise NotImplementedError('DEM is missing. An automatic DEM getter is currently not implemented.')
        return self._dem

    @dem.setter
    def dem(self, *geoArr_initArgs):
        self._dem = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'dem', nodataVal=0)  # FIXME 0?

    @dem.deleter
    def dem(self):
        self._dem = None

    @property
    def deadpixelmap(self) -> GeoArray:
        """Return the dead pixel map.

        :return: geoarray.GeoArray
        """
        if self._deadpixelmap is not None:
            self._deadpixelmap.arr = self._deadpixelmap[:].astype(bool)  # ensure boolean map

        return self._deadpixelmap

    @deadpixelmap.setter
    def deadpixelmap(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dpm = GeoArray(*geoArr_initArgs)

            if dpm.ndim == 3 and dpm.shape != self.data.shape:
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with a 3D array with the same size "
                                 "like %s.data, i.e.: %s Received %s."
                                 % (self.__class__.__name__, str(self.data.shape), str(dpm.shape)))
            elif dpm.ndim == 2 and dpm.shape != (self.data.bands, self.data.cols):
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with an array with the size "
                                 "'bands x columns' of the GeoArray %s.data. Received %s. Expected %s"
                                 % (self.__class__.__name__, str(dpm.shape), str((self.data.bands, self.data.cols))))

            self._deadpixelmap = dpm
        else:
            del self.deadpixelmap

    @deadpixelmap.deleter
    def deadpixelmap(self):
        self._deadpixelmap = None

    def _get_geoarray_with_datalike_geometry(self,
                                             geoArr_initArgs: tuple,
                                             attrName: str,
                                             nodataVal: int = None,
                                             specialclass=None) -> Optional[GeoArray]:
        if geoArr_initArgs[0] is None:
            return None
        else:
            GeoArrayOrSubclass = GeoArray if not specialclass else specialclass
            gA = GeoArrayOrSubclass(*geoArr_initArgs)

            if gA.shape[:2] != self.data.shape[:2]:
                raise ValueError("The '%s' GeoArray can only be instanced with an array with "
                                 "the same X/Y dimensions like %s.data %s. Got %s." %
                                 (attrName, self.__class__.__name__, str(self.data.shape[:2]), str(gA.shape[:2])))

            # noinspection PyProtectedMember
            if gA._nodata is None and nodataVal is not None:
                gA.nodata = nodataVal
            gA.gt = self.data.gt
            gA.prj = self.data.prj

            return gA

    def generate_quicklook(self, bands2use: List[int]) -> GeoArray:
        """
        Generate image quicklook and save into a file.

        :param bands2use:   band indices of self.data to be used as (red, green, blue) for quicklook image,
                            e.g., [3, 2, 1]
        :return: GeoArray
        """
        red, green, blue = [self.data[:, :, bandidx] for bandidx in bands2use]

        def rescale(bandarray):
            pixvals = np.ma.masked_equal(bandarray, self.data.nodata).compressed() \
                if self.data.nodata is not None else bandarray
            p2, p98 = np.nanpercentile(pixvals, 2), np.nanpercentile(pixvals, 98)

            return exposure.rescale_intensity(bandarray, in_range=(p2, p98), out_range=(0, 255))

        pix = np.dstack((rescale(red), rescale(green), rescale(blue))).astype(np.uint8)

        return GeoArray(pix)
