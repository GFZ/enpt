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

"""EnPT images module. All EnMAP image objects are defined here."""

import logging
from types import SimpleNamespace
from typing import Tuple  # noqa: F401
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import numpy as np
from os import path, makedirs
from glob import glob
import utm
from scipy.interpolate import interp2d

# noinspection PyPackageRequirements
from skimage import exposure  # contained in package requirements as scikit-image

from geoarray import GeoArray, NoDataMask, CloudMask

from ..utils.logging import EnPT_Logger
from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo, EnMAP_Metadata_L1B_Detector_SensorGeo
from ..model.metadata import EnMAP_Metadata_L2A_MapGeo  # noqa: F401  # only used for type hint
from ..options.config import EnPTConfig
from ..processors.dead_pixel_correction import Dead_Pixel_Corrector
from ..processors.dem_preprocessing import DEM_Processor
from ..processors.spatial_transform import compute_mapCoords_within_sensorGeoDims

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


################
# BASE CLASSES #
################


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
        p2 = np.percentile(self.data[:, :, bands2use[0]], 2)
        p98 = np.percentile(self.data[:, :, bands2use[0]], 98)
        red_rescaled = exposure.rescale_intensity(self.data[:, :, bands2use[0]], (p2, p98))
        p2 = np.percentile(self.data[:, :, bands2use[1]], 2)
        p98 = np.percentile(self.data[:, :, bands2use[1]], 98)
        green_rescaled = exposure.rescale_intensity(self.data[:, :, bands2use[1]], (p2, p98))
        p2 = np.percentile(self.data[:, :, bands2use[2]], 2)
        p98 = np.percentile(self.data[:, :, bands2use[2]], 98)
        blue_rescaled = exposure.rescale_intensity(self.data[:, :, bands2use[2]], (p2, p98))
        pix = np.dstack((red_rescaled, green_rescaled, blue_rescaled))
        pix = np.uint8(pix * 255)

        return GeoArray(pix)

#######################################
# EnPT EnMAP objects in sensor geometry
#######################################


class EnMAP_Detector_SensorGeo(_EnMAP_Image):
    """Class representing a single detector of an EnMAP image (in sensor geometry).

    NOTE:
        - Inherits all attributes from _EnMAP_Image class.
        - All functionality that VNIR and SWIR detectors (sensor geometry) have in common is to be implemented here.

    Attributes:
        - to be listed here. Check help(_EnMAP_Detector_SensorGeo) in the meanwhile!

    """

    def __init__(self, detector_name: str, root_dir: str, config: EnPTConfig, logger=None, meta=None):
        """Get an instance of _EnMAP_Detector_SensorGeo.

        :param detector_name:   'VNIR' or 'SWIR'
        :param root_dir:
        :param logger:
        :param meta: import meta if already loaded
        """
        if detector_name not in ['VNIR', 'SWIR']:
            raise ValueError("'detector_name' must be 'VNIR' or 'SWIR'. Received %s!" % detector_name)

        self.cfg = config
        self._root_dir = root_dir
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(EnMAP_Detector_SensorGeo, self).__init__()

        # instance an empty metadata object
        if meta is None:
            self.detector_meta = \
                EnMAP_Metadata_L1B_Detector_SensorGeo(self.detector_name, config=self.cfg, logger=self.logger)
        else:
            self.detector_meta = meta  # type: EnMAP_Metadata_L1B_Detector_SensorGeo

    def get_paths(self):
        """
        Get all file paths associated with the current instance of EnMAP_Detector_SensorGeo
        These information are read from the detector_meta.
        :return: paths as SimpleNamespace
        """
        paths = SimpleNamespace()
        paths.root_dir = self._root_dir
        paths.data = path.join(self._root_dir, self.detector_meta.data_filename)

        paths.mask_clouds = path.join(self._root_dir, self.detector_meta.cloud_mask_filename) \
            if self.detector_meta.cloud_mask_filename else None
        paths.deadpixelmap = path.join(self._root_dir, self.detector_meta.dead_pixel_filename) \
            if self.detector_meta.dead_pixel_filename else None

        paths.quicklook = path.join(self._root_dir, self.detector_meta.quicklook_filename)

        return paths

    def correct_dead_pixels(self):
        """Correct dead pixels with respect to the dead pixel mask."""
        algo = self.cfg.deadpix_P_algorithm
        method_spectral, method_spatial = self.cfg.deadpix_P_interp_spectral, self.cfg.deadpix_P_interp_spatial
        self.logger.info("Correcting dead pixels of %s detector...\n"
                         "Used algorithm: %s interpolation in the %s domain"
                         % (self.detector_name, algo, method_spectral if algo == 'spectral' else method_spatial))

        self.data = \
            Dead_Pixel_Corrector(algorithm=algo,
                                 interp_spectral=method_spectral,
                                 interp_spatial=method_spatial,
                                 logger=self.logger)\
            .correct(self.data, self.deadpixelmap)

    def get_preprocessed_dem(self):
        if self.cfg.path_dem:
            self.logger.info('Pre-processing DEM for %s...' % self.detector_name)
            DP = DEM_Processor(self.cfg.path_dem, enmapIm_cornerCoords=tuple(zip(self.detector_meta.lon_UL_UR_LL_LR,
                                                                                 self.detector_meta.lat_UL_UR_LL_LR)),
                               CPUs=self.cfg.CPUs)
            DP.fill_gaps()  # FIXME this will also be needed at other places

            R, C = self.data.shape[:2]
            if DP.dem.is_map_geo:
                lons = self.detector_meta.lons
                lats = self.detector_meta.lats

                if not (lons.ndim == 2 and lats.ndim == 2) and not (lons.ndim == 3 and lats.ndim == 3):
                    raise ValueError((lons.ndim, lats.ndim), 'Geolayer must be either 2- or 3-dimensional.')

                msg_bandinfo = ''
                if lons.ndim == 3:
                    # 3D geolayer (the usual case for EnMAP data provided by DLR)
                    lons = lons[:, :, 0]
                    lats = lats[:, :, 0]
                    msg_bandinfo = ' (using first band of %s geolayer)' % self.detector_name
                else:
                    # 2D geolayer (GFZ-internal test case)
                    # FIXME replace linear interpolation by native geolayers
                    if lons.shape != self.data.shape:
                        lons = self.detector_meta.interpolate_corners(*self.detector_meta.lon_UL_UR_LL_LR, nx=C, ny=R)
                    if lats.shape != self.data.shape:
                        lats = self.detector_meta.interpolate_corners(*self.detector_meta.lat_UL_UR_LL_LR, nx=C, ny=R)

                self.logger.info(('Transforming DEM to %s sensor geometry%s...' % (self.detector_name, msg_bandinfo)))
                self.dem = DP.to_sensor_geometry(lons=lons, lats=lats)
            else:
                self.dem = DP.dem

        return self.dem

    def append_new_image(self, img2: 'EnMAP_Detector_SensorGeo', n_lines: int = None) -> None:
        # TODO convert method to function?
        """Check if a second image matches with the first image and if so, append the given number of lines below.

        In this version we assume that the image will be added below. If it is the case, the method will create
        temporary files that will be used in the following.

        :param img2:
        :param n_lines: number of lines to be added from the new image
        :return: None
        """
        if not self.cfg.is_dummy_dataformat:
            basename_img1 = self.detector_meta.data_filename.split('-SPECTRAL_IMAGE')[0] + '::%s' % self.detector_name
            basename_img2 = img2.detector_meta.data_filename.split('-SPECTRAL_IMAGE')[0] + '::%s' % img2.detector_name
        else:
            basename_img1 = path.basename(self._root_dir)
            basename_img2 = path.basename(img2._root_dir)

        self.logger.info("Check new image for %s: %s " % (self.detector_name, basename_img2))

        distance_min = 27.0
        distance_max = 34.0

        # compute distance between image1 LL and image2 UL
        x1, y1, _, _ = utm.from_latlon(self.detector_meta.lat_UL_UR_LL_LR[2], self.detector_meta.lon_UL_UR_LL_LR[2])
        x2, y2, _, _ = utm.from_latlon(img2.detector_meta.lat_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[0])
        distance_left = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        # compute distance between image1 LR and image2 UR
        x1, y1, _, _ = utm.from_latlon(self.detector_meta.lat_UL_UR_LL_LR[3], self.detector_meta.lon_UL_UR_LL_LR[3])
        x2, y2, _, _ = utm.from_latlon(img2.detector_meta.lat_UL_UR_LL_LR[1], img2.detector_meta.lon_UL_UR_LL_LR[1])
        distance_right = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        if distance_min < distance_left < distance_max and distance_min < distance_right < distance_max:
            self.logger.info("Append new image to %s: %s" % (self.detector_name, basename_img2))
        else:
            self.logger.warning("%s and %s don't fit to be appended." % (basename_img1, basename_img2))
            return

        # set new number of line
        n_lines = n_lines or img2.detector_meta.nrows

        if n_lines > img2.detector_meta.nrows:
            self.logger.warning("n_lines (%s) exceeds the total number of line of second image" % n_lines)
            self.logger.warning("Set to the image number of line")
            n_lines = img2.detector_meta.nrows

        if n_lines < 50:  # TODO: determine these values
            self.logger.warning("A minimum of 50 lines is required, only %s were selected" % n_lines)
            self.logger.warning("Set the number of line to 50")
            n_lines = 50

        self.detector_meta.nrows += n_lines

        # Compute new lower coordinates
        if not self.cfg.is_dummy_dataformat:
            img2_cornerCoords = tuple(zip(img2.detector_meta.lon_UL_UR_LL_LR,
                                          img2.detector_meta.lat_UL_UR_LL_LR))
            dem_validated = DEM_Processor(img2.cfg.path_dem,
                                          enmapIm_cornerCoords=img2_cornerCoords).dem
            LL, LR = compute_mapCoords_within_sensorGeoDims(
                sensorgeoCoords_YX=[(n_lines - 1, 0),  # LL
                                    (n_lines - 1, img2.detector_meta.ncols - 1)],  # LR
                rpc_coeffs=list(img2.detector_meta.rpc_coeffs.values())[0],  # RPC coeffs of first band of the detector
                dem=dem_validated,
                enmapIm_cornerCoords=img2_cornerCoords,
                enmapIm_dims_sensorgeo=(img2.detector_meta.nrows, img2.detector_meta.ncols)
            )

            self.detector_meta.lon_UL_UR_LL_LR[2], self.detector_meta.lat_UL_UR_LL_LR[2] = LL
            self.detector_meta.lon_UL_UR_LL_LR[3], self.detector_meta.lat_UL_UR_LL_LR[3] = LR
        else:
            # lats
            ff = interp2d(x=[0, 1],
                          y=[0, 1],
                          z=[[img2.detector_meta.lat_UL_UR_LL_LR[0], img2.detector_meta.lat_UL_UR_LL_LR[1]],
                             [img2.detector_meta.lat_UL_UR_LL_LR[2], img2.detector_meta.lat_UL_UR_LL_LR[3]]],
                          kind='linear')
            self.detector_meta.lat_UL_UR_LL_LR[2] = np.array(ff(0, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lat_UL_UR_LL_LR[3] = np.array(ff(1, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lats = self.detector_meta.interpolate_corners(*self.detector_meta.lat_UL_UR_LL_LR,
                                                                             self.detector_meta.ncols,
                                                                             self.detector_meta.nrows)
            # lons
            ff = interp2d(x=[0, 1],
                          y=[0, 1],
                          z=[[img2.detector_meta.lon_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[1]],
                             [img2.detector_meta.lon_UL_UR_LL_LR[2], img2.detector_meta.lon_UL_UR_LL_LR[3]]],
                          kind='linear')
            self.detector_meta.lon_UL_UR_LL_LR[2] = np.array(ff(0, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lon_UL_UR_LL_LR[3] = np.array(ff(1, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lons = self.detector_meta.interpolate_corners(*self.detector_meta.lon_UL_UR_LL_LR,
                                                                             self.detector_meta.ncols,
                                                                             self.detector_meta.nrows)

        # append the raster data
        self.data = np.append(self.data, img2.data[0:n_lines, :, :], axis=0)
        self.mask_clouds = np.append(self.mask_clouds, img2.mask_clouds[0:n_lines, :], axis=0)
        if not self.cfg.is_dummy_dataformat:
            self.deadpixelmap = np.append(self.deadpixelmap, img2.deadpixelmap[0:n_lines, :], axis=0)
        # TODO append remaining raster layers - additional cloud masks, ...

        # NOTE: We leave the quicklook out here because merging the quicklook of adjacent scenes might cause a
        #       brightness jump that can be avoided by recomputing the quicklook after DN/radiance conversion.

    def DN2TOARadiance(self):
        """Convert DNs to TOA radiance.

        Convert the radiometric unit of _EnMAP_Detector_SensorGeo.data from digital numbers to top-of-atmosphere
        radiance.

        :return: None
        """
        # TODO move to processors.radiometric_transform?
        if self.detector_meta.unitcode == 'DN':
            self.logger.info('Converting DN values to radiance [mW/m^2/sr/nm] for %s detector...' % self.detector_name)

            if self.detector_meta.l_min is not None and self.detector_meta.l_max is not None:
                # Lλ = (LMINλ + ((LMAXλ - LMINλ)/(QCALMAX-QCALMIN)) * (QCAL-QCALMIN))
                # FIXME this asserts LMIN and LMAX in mW m^-2 sr^-1 nm^-1

                QCALMIN = 1
                QCALMAX = 2 ** 16  # 65535 (16-bit maximum value)
                LMIN = self.detector_meta.l_min
                LMAX = self.detector_meta.l_max
                QCAL = self.data[:]

                self.data = ((LMAX - LMIN)/(QCALMAX - QCALMIN)) * (QCAL - QCALMIN) + LMIN

            elif self.detector_meta.gains is not None and self.detector_meta.offsets is not None:
                # Lλ = QCAL * GAIN + OFFSET
                # NOTE: - DLR provides gains between 2000 and 10000, so we have to DEVIDE by gains
                #       - DLR gains / offsets are provided in W/m2/sr/nm, so we have to multiply by 1000 to get
                #         mW/m2/sr/nm as needed later
                self.data = 1000 * (self.data[:] * self.detector_meta.gains + self.detector_meta.offsets)

            else:
                raise ValueError("Neighter 'l_min'/'l_max' nor 'gains'/'offsets' "
                                 "are available for radiance computation.")

            self.detector_meta.unit = "mW m^-2 sr^-1 nm^-1"
            self.detector_meta.unitcode = "TOARad"
        else:
            self.logger.warning(
                "No DN to Radiance conversion is performed because unitcode is not DN (found: {code}).".format(
                    code=self.detector_meta.unitcode)
            )


class EnMAPL1Product_SensorGeo(object):
    """Class for EnPT EnMAP object in sensor geometry.

    Attributes:
        - logger:
            - logging.Logger instance or subclass instance
        - vnir
            - instance of EnMAP_VNIR_SensorGeo class
        - swir
            - instance of EnMAP_SWIR_SensorGeo class
        - paths:
            - paths belonging to the EnMAP product
        - meta:
            - instance of EnMAP_Metadata_SensorGeo class
        - detector_attrNames:
            - list of attribute names for VNIR and SWIR detectors,

    """

    def __init__(self, root_dir: str, config: EnPTConfig, logger=None):
        """Get instance of EnPT EnMAP object in sensor geometry.

        :param root_dir:    Root directory of EnMAP Level-1B product
        :param config:      EnPT configuration object
        :param logger:      None or logging instance to be appended to EnMAPL1Product_SensorGeo instance
                            (If None, a default EnPT_Logger is used).
        """
        # protected attributes
        self._logger = None

        # populate attributes
        self.cfg = config
        if logger:
            self.logger = logger

        # Read metadata here in order to get all information needed by to get paths.
        if not self.cfg.is_dummy_dataformat:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*METADATA.XML"))[0],
                                                     config=self.cfg, logger=self.logger)
        else:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*_header.xml"))[0],
                                                     config=self.cfg, logger=self.logger)
        self.meta.read_metadata()

        # define VNIR and SWIR detector
        self.detector_attrNames = ['vnir', 'swir']
        self.vnir = EnMAP_Detector_SensorGeo('VNIR', root_dir, config=self.cfg, logger=self.logger, meta=self.meta.vnir)
        self.swir = EnMAP_Detector_SensorGeo('SWIR', root_dir, config=self.cfg, logger=self.logger, meta=self.meta.swir)

        # Get the paths according information delivered in the metadata
        self.paths = self.get_paths()

        # associate raster attributes with file links (raster data is read lazily / on demand)
        self.vnir.data = self.paths.vnir.data
        self.vnir.mask_clouds = self.paths.vnir.mask_clouds
        self.swir.data = self.paths.swir.data
        self.swir.mask_clouds = self.paths.swir.mask_clouds

        try:
            self.vnir.deadpixelmap = self.paths.vnir.deadpixelmap
            self.swir.deadpixelmap = self.paths.swir.deadpixelmap
        except ValueError:
            self.logger.warning("Unexpected dimensions of dead pixel mask. Setting all pixels to 'normal'.")
            self.vnir.deadpixelmap = np.zeros(self.vnir.data.shape)
            self.swir.deadpixelmap = np.zeros(self.swir.data.shape)

        # NOTE: We leave the quicklook out here because merging the quicklook of adjacent scenes might cause a
        #       brightness jump that can be avoided by recomputing the quicklook after DN/radiance conversion.

        # compute radiance
        self.DN2TOARadiance()

    def get_paths(self):
        """
        Get all file paths associated with the current instance of EnMAPL1Product_SensorGeo

        :return: paths.SimpleNamespace()
        """
        paths = SimpleNamespace()
        paths.vnir = self.vnir.get_paths()
        paths.swir = self.swir.get_paths()
        paths.root_dir = self.meta.rootdir
        paths.metaxml = self.meta.path_xml

        return paths

    @property
    def logger(self) -> EnPT_Logger:
        """Get a an instance of enpt.utils.logging.EnPT_Logger.

        NOTE:
            - The logging level will be set according to the user inputs of EnPT.
            - The path of the log file is directly derived from the attributes of the _EnMAP_Image instance.

        Usage:
            - get the logger:
                logger = self.logger
            - set the logger
                self.logger = logging.getLogger()  # NOTE: only instances of logging.Logger are allowed here
            - delete the logger:
                del self.logger  # or "self.logger = None"

        :return: EnPT_Logger instance
        """
        if self._logger and self._logger.handlers[:]:
            return self._logger
        else:
            basename = path.splitext(path.basename(self.cfg.path_l1b_enmap_image))[0]
            path_logfile = path.join(self.cfg.output_dir, basename + '.log') \
                if self.cfg.create_logfile and self.cfg.output_dir else ''
            self._logger = EnPT_Logger('log__' + basename, fmt_suffix=None, path_logfile=path_logfile,
                                       log_level=self.cfg.log_level, append=False)

            return self._logger

    @logger.setter
    def logger(self, logger: logging.Logger):
        assert isinstance(logger, logging.Logger) or logger in ['not set', None], \
            "%s.logger can not be set to %s." % (self.__class__.__name__, logger)

        # save prior logs
        # if logger is None and self._logger is not None:
        #     self.log += self.logger.captured_stream
        self._logger = logger

    @property
    def log(self) -> str:
        """Return a string of all logged messages until now.

        NOTE: self.log can also be set to a string.
        """
        return self.logger.captured_stream

    @log.setter
    def log(self, string: str):
        assert isinstance(string, str), "'log' can only be set to a string. Got %s." % type(string)
        self.logger.captured_stream = string

    # @classmethod
    # def from_L1B_provider_data(cls, path_enmap_image: str, config: EnPTConfig=None) -> EnMAPL1Product_SensorGeo:
    #     """
    #
    #     :param path_enmap_image:
    #     :param config:
    #     :return:
    #     """
    #     # input validation
    #     if not path.isdir(path_enmap_image) and \
    #        not (path.exists(path_enmap_image) and path_enmap_image.endswith('.zip')):
    #         raise ValueError("The parameter 'path_enmap_image' must be a directory or the path to an existing zip "
    #                          "archive.")
    #
    #     # extract L1B image archive if needed
    #     if path_enmap_image.endswith('.zip'):
    #         path_enmap_image = self.extract_zip_archive(path_enmap_image)
    #         if not path.isdir(path_enmap_image):
    #             raise NotADirectoryError(path_enmap_image)
    #
    #     # run the reader
    #     from ..io.reader import L1B_Reader
    #     RD = L1B_Reader(config=config)
    #     L1_obj = RD.read_inputdata(path_enmap_image, observation_time=datetime(2015, 12, 7, 10))
    #
    #     return L1_obj

    def DN2TOARadiance(self):
        """Convert self.vnir.data and self.swir.data from digital numbers to top-of-atmosphere radiance.

        :return: None
        """
        if self.vnir.detector_meta.unitcode != 'TOARad':
            self.vnir.DN2TOARadiance()
            self.meta.vnir.unitcode = self.vnir.detector_meta.unitcode  # FIXME possible duplicate?
            self.meta.vnir.unit = self.vnir.detector_meta.unit  # FIXME possible duplicate?
        if self.swir.detector_meta.unitcode != 'TOARad':
            self.swir.DN2TOARadiance()
            self.meta.swir.unitcode = self.swir.detector_meta.unitcode
            self.meta.swir.unit = self.swir.detector_meta.unit

    def append_new_image(self, root_dir, n_line_ext):
        """Append a second EnMAP Level-1B product below the last line of the current L1B product.

        NOTE:   We create new files that will be saved into a temporary directory and their path will be stored in the
                paths of self. We assume that the dead pixel map does not change between two adjacent images.

        :param root_dir:    root directory of EnMAP Level-1B product to be appended
        :param n_line_ext:  number of lines to be added from the new image
        :return:
        """
        l1b_ext_obj = EnMAPL1Product_SensorGeo(root_dir, config=self.cfg)

        self.vnir.append_new_image(l1b_ext_obj.vnir, n_line_ext)
        self.swir.append_new_image(l1b_ext_obj.swir, n_line_ext)

    def calc_snr_from_radiance(self):
        """Compute EnMAP SNR from radiance data.

        SNR equation:    SNR = p0 + p1 * LTOA + p2 * LTOA ^ 2   [W/(m^2 sr nm)].
        """
        with TemporaryDirectory(dir=self.cfg.working_dir) as tmpDir, \
                ZipFile(self.cfg.path_l1b_snr_model, "r") as zf:

            zf.extractall(tmpDir)

            if self.meta.vnir.unitcode == 'TOARad':
                self.vnir.detector_meta.calc_snr_from_radiance(rad_data=self.vnir.data, dir_snr_models=tmpDir)

            if self.meta.swir.unitcode == 'TOARad':
                self.swir.detector_meta.calc_snr_from_radiance(rad_data=self.swir.data, dir_snr_models=tmpDir)

    def correct_dead_pixels(self):
        """Correct dead pixels of both detectors."""
        self.vnir.correct_dead_pixels()
        self.swir.correct_dead_pixels()

    # def correct_VNIR_SWIR_shift(self):
    #     # use first geolayer bands for VNIR and SWIR
    #     Vlons, Vlats = self.vnir.detector_meta.lons[:, :, 0], self.vnir.detector_meta.lats[:, :, 0]
    #     Slons, Slats = self.swir.detector_meta.lons[:, :, 0], self.swir.detector_meta.lats[:, :, 0]
    #
    #     # get corner coordinates of VNIR and SWIR according to geolayer
    #     def get_coords(lons, lats):
    #         return tuple([(lons[Y, X], lats[Y, X]) for Y, X in [(0, 0), (0, -1), (-1, 0), (-1, -1)]])
    #
    #     VUL, VUR, VLL, VLR = get_coords(Vlons, Vlats)
    #     SUL, SUR, SLL, SLR = get_coords(Slons, Slats)
    #
    #     # get map coordinates of VNIR/SWIR overlap
    #     ovUL = max(VUL[0], SUL[0]), min(VUL[1], SUL[1])
    #     ovUR = min(VUR[0], SUR[0]), min(VUR[1], SUR[1])
    #     ovLL = max(VLL[0], SLL[0]), max(VLL[1], SLL[1])
    #     ovLR = min(VLR[0], SLR[0]), max(VLR[1], SLR[1])
    #
    #     # find nearest image positions for VNIR and SWIR to the map coordinates of the VNIR/SWIR overlap
    #     def nearest_imCoord(lons_arr, lats_arr, lon, lat):
    #         dists = np.sqrt((lons_arr - lon) ** 2 + (lats_arr - lat) ** 2)
    #         row, col = np.unravel_index(dists.argmin(), dists.shape)
    #
    #         return row, col
    #
    #     overlapImVNIR = tuple([nearest_imCoord(Vlons, Vlats, *ovCoords) for ovCoords in [ovUL, ovUR, ovLL, ovLR]])
    #     overlapImSWIR = tuple([nearest_imCoord(Slons, Slats, *ovCoords) for ovCoords in [ovUL, ovUR, ovLL, ovLR]])
    #
    #     raise NotImplementedError()  # FIXME
    #     self.vnir.data.get_mapPos()  # FIXME

    def get_preprocessed_dem(self):
        self.vnir.get_preprocessed_dem()
        self.swir.get_preprocessed_dem()

    def run_AC(self):
        from ..processors.atmospheric_correction import AtmosphericCorrector
        AC = AtmosphericCorrector(config=self.cfg)
        AC.run_ac(self)

    def save(self, outdir: str, suffix="") -> str:
        """Save the product to disk using almost the same input format.

        NOTE: Radiance is saved instead of DNs to avoid a radiometric jump between concatenated images.

        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        product_dir = path.join(path.abspath(outdir),
                                "{name}{suffix}".format(name=self.meta.scene_basename, suffix=suffix))

        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        # write the VNIR
        self.vnir.data.save(product_dir + path.sep + self.meta.vnir.data_filename, fmt="ENVI")
        self.vnir.mask_clouds.save(product_dir + path.sep + self.meta.vnir.cloud_mask_filename, fmt="GTiff")
        if self.vnir.deadpixelmap is not None:
            self.vnir.deadpixelmap.save(product_dir + path.sep + self.meta.vnir.dead_pixel_filename, fmt="GTiff")
        else:
            self.logger.warning('Could not save VNIR dead pixel map because there is no corresponding attribute.')

        # FIXME we could also write the quicklook included in DLR L1B format
        self.vnir.generate_quicklook(bands2use=self.vnir.detector_meta.preview_bands) \
            .save(path.join(product_dir, path.basename(self.meta.vnir.quicklook_filename) + '.png'), fmt='PNG')

        # write the SWIR
        self.swir.data.save(product_dir + path.sep + self.meta.swir.data_filename, fmt="ENVI")
        self.swir.mask_clouds.save(product_dir + path.sep + self.meta.swir.cloud_mask_filename, fmt="GTiff")
        if self.swir.deadpixelmap is not None:
            self.swir.deadpixelmap.save(product_dir + path.sep + self.meta.swir.dead_pixel_filename, fmt="GTiff")
        else:
            self.logger.warning('Could not save SWIR dead pixel map because there is no corresponding attribute.')
        self.swir.generate_quicklook(bands2use=self.swir.detector_meta.preview_bands) \
            .save(path.join(product_dir, path.basename(self.meta.swir.quicklook_filename) + '.png'), fmt='PNG')

        # Update the xml file
        metadata_string = self.meta.to_XML()
        with open(product_dir + path.sep + path.basename(self.paths.metaxml), "w") as xml_file:
            xml_file.write(metadata_string)

        self.logger.info("L1B product successfully written!")

        return product_dir


####################################
# EnPT EnMAP objects in map geometry
####################################


class EnMAP_Detector_MapGeo(_EnMAP_Image):
    """Base class representing a single detector of an EnMAP image (as map geometry).

    NOTE:
        - Inherits all attributes from _EnMAP_Image class.
        - All functionality that VNIR and SWIR detectors (map geometry) have in common is to be implemented here.
        - All EnMAP image subclasses representing a specific EnMAP detector (sensor geometry) should inherit from
          _EnMAP_Detector_SensorGeo.

    Attributes:
        - to be listed here. Check help(_EnMAP_Detector_SensorGeo) in the meanwhile!

    """

    def __init__(self, detector_name: str, logger=None):
        """Get an instance of _EnMAP_Detector_MapGeo.

        :param detector_name:   'VNIR' or 'SWIR'
        :param logger:
        """
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(EnMAP_Detector_MapGeo, self).__init__()

    @property
    def mask_nodata(self) -> GeoArray:
        """Return the no data mask.

        Bundled with all the corresponding metadata.

        For usage instructions and a list of attributes refer to help(self.data).
        self.mask_nodata works in the same way.

        :return: instance of geoarray.NoDataMask
        """
        if self._mask_nodata is None and isinstance(self.data, GeoArray):
            self.logger.info('Calculating nodata mask...')
            self._mask_nodata = self.data.mask_nodata  # calculates mask nodata if not already present

        return self._mask_nodata

    @mask_nodata.setter
    def mask_nodata(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            nd = NoDataMask(*geoArr_initArgs)
            if nd.shape[:2] != self.data.shape[:2]:
                raise ValueError("The 'mask_nodata' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(nd.shape))
            nd.nodata = False
            nd.gt = self.data.gt
            nd.prj = self.data.prj
            self._mask_nodata.prj = nd
        else:
            del self.mask_nodata

    @mask_nodata.deleter
    def mask_nodata(self):
        self._mask_nodata = None

    def calc_mask_nodata(self, fromBand=None, overwrite=False) -> GeoArray:
        """Calculate a no data mask with (values: 0=nodata; 1=data).

        :param fromBand:  <int> index of the band to be used (if None, all bands are used)
        :param overwrite: <bool> whether to overwrite existing nodata mask that has already been calculated
        :return:
        """
        self.logger.info('Calculating nodata mask...')

        if self._mask_nodata is None or overwrite:
            self.data.calc_mask_nodata(fromBand=fromBand, overwrite=overwrite)
            self.mask_nodata = self.data.mask_nodata
            return self.mask_nodata


class EnMAPL2Product_MapGeo(_EnMAP_Image):
    """Class for EnPT Level-2 EnMAP object in map geometry.

    Attributes:
        - logger:
            - logging.Logger instance or subclass instance
        - paths:
            - paths belonging to the EnMAP product
        - meta:
            - instance of EnMAP_Metadata_SensorGeo class
    """
    def __init__(self, config: EnPTConfig, logger=None):
        # protected attributes
        self._logger = None

        # populate attributes
        self.cfg = config
        if logger:
            self.logger = logger

        self.meta = None  # type: EnMAP_Metadata_L2A_MapGeo
        self.paths = None  # type: SimpleNamespace

        super(EnMAPL2Product_MapGeo, self).__init__()

    @property
    def logger(self) -> EnPT_Logger:
        """Get a an instance of enpt.utils.logging.EnPT_Logger.

        NOTE:
            - The logging level will be set according to the user inputs of EnPT.
            - The path of the log file is directly derived from the attributes of the _EnMAP_Image instance.

        Usage:
            - get the logger:
                logger = self.logger
            - set the logger
                self.logger = logging.getLogger()  # NOTE: only instances of logging.Logger are allowed here
            - delete the logger:
                del self.logger  # or "self.logger = None"

        :return: EnPT_Logger instance
        """
        if self._logger and self._logger.handlers[:]:
            return self._logger
        else:
            basename = path.splitext(path.basename(self.cfg.path_l1b_enmap_image))[0]
            path_logfile = path.join(self.cfg.output_dir, basename + '.log') \
                if self.cfg.create_logfile and self.cfg.output_dir else ''
            self._logger = EnPT_Logger('log__' + basename, fmt_suffix=None, path_logfile=path_logfile,
                                       log_level=self.cfg.log_level, append=False)

            return self._logger

    @logger.setter
    def logger(self, logger: logging.Logger):
        assert isinstance(logger, logging.Logger) or logger in ['not set', None], \
            "%s.logger can not be set to %s." % (self.__class__.__name__, logger)

        # save prior logs
        # if logger is None and self._logger is not None:
        #     self.log += self.logger.captured_stream
        self._logger = logger

    @property
    def log(self) -> str:
        """Return a string of all logged messages until now.

        NOTE: self.log can also be set to a string.
        """
        return self.logger.captured_stream

    @log.setter
    def log(self, string: str):
        assert isinstance(string, str), "'log' can only be set to a string. Got %s." % type(string)
        self.logger.captured_stream = string

    @classmethod
    def from_L1B_sensorgeo(cls, config: EnPTConfig, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        from ..processors.orthorectification import Orthorectifier
        L2_obj = Orthorectifier(config=config).run_transformation(enmap_ImageL1=enmap_ImageL1)

        return L2_obj

    def get_paths(self, l2a_outdir: str):
        """
        Get all file paths associated with the current instance of EnMAP_Detector_SensorGeo
        These information are read from the detector_meta.

        :param l2a_outdir:  output directory of EnMAP Level-2A dataset
        :return: paths as SimpleNamespace
        """
        paths = SimpleNamespace()
        paths.root_dir = l2a_outdir
        paths.metaxml = path.join(l2a_outdir, self.meta.metaxml_filename)
        paths.data = path.join(l2a_outdir, self.meta.data_filename)
        paths.mask_clouds = path.join(l2a_outdir, self.meta.cloud_mask_filename)
        paths.deadpixelmap_vnir = path.join(l2a_outdir, self.meta.dead_pixel_filename_vnir)
        paths.deadpixelmap_swir = path.join(l2a_outdir, self.meta.dead_pixel_filename_swir)
        paths.quicklook_vnir = path.join(l2a_outdir, self.meta.quicklook_filename_vnir)
        paths.quicklook_swir = path.join(l2a_outdir, self.meta.quicklook_filename_swir)

        return paths

    def save(self, outdir: str, suffix="") -> str:
        """
        Save the product to disk using almost the same input format
        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        # TODO optionally add more output formats
        product_dir = path.join(path.abspath(outdir),
                                "{name}{suffix}".format(name=self.meta.scene_basename, suffix=suffix))

        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        # define output paths
        outpath_data = path.join(product_dir, self.meta.data_filename)
        outpath_mask_clouds = path.join(product_dir, self.meta.cloud_mask_filename)
        outpath_quicklook_vnir = path.join(product_dir, self.meta.quicklook_filename_vnir)
        outpath_quicklook_swir = path.join(product_dir, self.meta.quicklook_filename_swir)
        outpath_meta = path.join(product_dir, self.meta.metaxml_filename)
        outpaths = [outpath_data, outpath_mask_clouds, outpath_quicklook_vnir, outpath_quicklook_swir, outpath_meta]

        # save raster data
        kwargs_save = dict(fmt='GTiff', creationOptions=["COMPRESS=LZW"])
        self.data.save(outpath_data, **kwargs_save)
        self.mask_clouds.save(outpath_mask_clouds, **kwargs_save)

        # TODO VNIR and SWIR
        # self.deadpixelmap.save(path.join(product_dir, self.meta.cloud_mask_filename), **kwargs_save)
        self.logger.warning('Currently, L2A dead pixel masks cannot be saved yet.')

        self.generate_quicklook(bands2use=self.meta.preview_bands_vnir).save(outpath_quicklook_vnir, **kwargs_save)
        self.generate_quicklook(bands2use=self.meta.preview_bands_swir).save(outpath_quicklook_swir, **kwargs_save)

        # TODO remove GDAL's *.aux.xml files?

        # save metadata
        self.meta.add_product_fileinformation(filepaths=outpaths)
        metadata_string = self.meta.to_XML()

        with open(outpath_meta, 'w') as metaF:
            self.logger.info("Writing metdata to %s" % outpath_meta)
            metaF.write(metadata_string)

        self.logger.info("L2A product successfully written!")

        return product_dir
