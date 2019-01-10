# -*- coding: utf-8 -*-
"""EnPT images module. All EnMAP image objects are defined here."""

import logging
from types import SimpleNamespace
import numpy as np
from os import path, makedirs
from xml.etree import ElementTree
from glob import glob
import utm
from scipy.interpolate import interp2d

# noinspection PyPackageRequirements
from skimage import exposure  # contained in package requirements as scikit-image

from geoarray import GeoArray, NoDataMask, CloudMask

from ..utils.logging import EnPT_Logger
from ..model.metadata import \
    EnMAP_Metadata_L1B_SensorGeo, \
    EnMAP_Metadata_L1B_Detector_SensorGeo, \
    EnMAP_Metadata_L2A_MapGeo, \
    L1B_product_props, \
    L1B_product_props_DLR
from ..options.config import EnPTConfig
from ..processors.dead_pixel_correction import Dead_Pixel_Corrector
from ..processors.dem_preprocessing import DEM_Processor
from ..processors.spatial_transform import compute_mapCoords_within_sensorGeoDims


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
        if self.cfg.is_dlr_dataformat:
            paths.metaxml = glob(path.join(self._root_dir, "*METADATA.XML"))[0]
        else:
            paths.metaxml = glob(path.join(self._root_dir, "*_header.xml"))[0]

        paths.data = path.join(self._root_dir, self.detector_meta.data_filename)

        paths.mask_clouds = path.join(self._root_dir, self.detector_meta.cloud_mask_filename) \
            if self.detector_meta.cloud_mask_filename else None
        paths.deadpixelmap = path.join(self._root_dir, self.detector_meta.dead_pixel_filename) \
            if self.detector_meta.dead_pixel_filename else None

        paths.quicklook = path.join(self._root_dir, self.detector_meta.quicklook_filename)
        return paths

    def correct_dead_pixels(self):
        """Correct dead pixels with respect to the dead pixel mask."""
        self.logger.info("Correcting dead pixels of %s detector..." % self.detector_name)

        self.data = \
            Dead_Pixel_Corrector(algorithm=self.cfg.deadpix_P_algorithm,
                                 interp=self.cfg.deadpix_P_interp,
                                 logger=self.logger)\
            .correct(self.data, self.deadpixelmap, progress=False if self.cfg.disable_progress_bars else True)

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
                    lons = lons[:, :, 0]
                    lats = lats[:, :, 0]
                    msg_bandinfo = ' (using first band of %s geolayer)' % self.detector_name
                else:
                    # 2D geolayer
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
        """
        Check if a second image could pass with the first image.
        In this version we assume that the image will be add below
        If it is the case, the method will create temporary files that will be used in the following.

        :param img2:
        :param n_lines: number of line to be added
        :return: None
        """
        if self.cfg.is_dlr_dataformat:
            basename_img1 = self.detector_meta.data_filename.split('-SPECTRAL_IMAGE')[0] + '::%s' % self.detector_name
            basename_img2 = img2.detector_meta.data_filename.split('-SPECTRAL_IMAGE')[0] + '::%s' % img2.detector_name
        else:
            basename_img1 = path.basename(self._root_dir)
            basename_img2 = path.basename(img2._root_dir)

        self.logger.info("Check new image for %s: %s " % (self.detector_name, basename_img2))

        distance_min = 27.0
        distance_max = 34.0

        # check bottom left
        x1, y1, _, _ = utm.from_latlon(self.detector_meta.lat_UL_UR_LL_LR[2], self.detector_meta.lon_UL_UR_LL_LR[2])
        x2, y2, _, _ = utm.from_latlon(img2.detector_meta.lat_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[0])
        distance_left = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        # check bottom right
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

        # Create new corner coordinate
        if self.cfg.is_dlr_dataformat:
            enmapIm_cornerCoords = tuple(zip(img2.detector_meta.lon_UL_UR_LL_LR,
                                             img2.detector_meta.lat_UL_UR_LL_LR))
            dem_validated = DEM_Processor(img2.cfg.path_dem,
                                          enmapIm_cornerCoords=enmapIm_cornerCoords).dem
            LL, LR = compute_mapCoords_within_sensorGeoDims(
                rpc_coeffs=list(img2.detector_meta.rpc_coeffs.values())[0],  # RPC coeffs of first band of the detector
                dem=dem_validated,
                enmapIm_cornerCoords=enmapIm_cornerCoords,
                enmapIm_dims_sensorgeo=(img2.detector_meta.nrows, img2.detector_meta.ncols),
                sensorgeoCoords_YX=[(n_lines - 1, 0),  # LL
                                    (n_lines - 1, img2.detector_meta.ncols - 1)]  # LR
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
            lon_lat_smpl = (15, 15)
            self.detector_meta.lats = self.detector_meta.interpolate_corners(*self.detector_meta.lat_UL_UR_LL_LR,
                                                                             *lon_lat_smpl)
            # lons
            ff = interp2d(x=[0, 1],
                          y=[0, 1],
                          z=[[img2.detector_meta.lon_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[1]],
                             [img2.detector_meta.lon_UL_UR_LL_LR[2], img2.detector_meta.lon_UL_UR_LL_LR[3]]],
                          kind='linear')
            self.detector_meta.lon_UL_UR_LL_LR[2] = np.array(ff(0, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lon_UL_UR_LL_LR[3] = np.array(ff(1, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lons = self.detector_meta.interpolate_corners(*self.detector_meta.lon_UL_UR_LL_LR,
                                                                             *lon_lat_smpl)

        # append the raster data
        self.data = np.append(self.data, img2.data[0:n_lines, :, :], axis=0)
        self.mask_clouds = np.append(self.mask_clouds, img2.mask_clouds[0:n_lines, :], axis=0)
        if self.cfg.is_dlr_dataformat:
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
                # Lλ = QCAL / GAIN + OFFSET
                # FIXME this asserts LMIN and LMAX in mW/cm2/sr/um
                # NOTE: - DLR provides gains between 2000 and 10000, so we have to DEVIDE by gains
                #       - DLR gains / offsets are provided in mW/cm2/sr/um, so we have to multiply by 10 to get
                #         mW m^-2 sr^-1 nm^-1 as needed later
                self.data = 10 * self.data[:] / self.detector_meta.gains + self.detector_meta.offsets

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

    def generate_quicklook(self) -> GeoArray:
        """
        Generate image quicklook and save into a file
        :return: GeoArray
        """
        p2 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[0]], 2)
        p98 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[0]], 98)
        red_rescaled = exposure.rescale_intensity(self.data[:, :, self.detector_meta.preview_bands[0]], (p2, p98))
        p2 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[1]], 2)
        p98 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[1]], 98)
        green_rescaled = exposure.rescale_intensity(self.data[:, :, self.detector_meta.preview_bands[1]], (p2, p98))
        p2 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[2]], 2)
        p98 = np.percentile(self.data[:, :, self.detector_meta.preview_bands[2]], 98)
        blue_rescaled = exposure.rescale_intensity(self.data[:, :, self.detector_meta.preview_bands[2]], (p2, p98))
        pix = np.dstack((red_rescaled, green_rescaled, blue_rescaled))
        pix = np.uint8(pix * 255)

        return GeoArray(pix)


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

    def __init__(self, root_dir: str, config: EnPTConfig, logger=None, lon_lat_smpl=None):
        """Get instance of EnPT EnMAP object in sensor geometry.

        :param root_dir:    Root directory of EnMAP Level-1B product
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
        if self.cfg.is_dlr_dataformat:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*METADATA.XML"))[0],
                                                     config=self.cfg, logger=self.logger)
        else:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*_header.xml"))[0],
                                                     config=self.cfg, logger=self.logger)
        self.meta.read_metadata(lon_lat_smpl)

        # define VNIR and SWIR detector
        self.vnir = EnMAP_Detector_SensorGeo('VNIR', root_dir, config=self.cfg, logger=self.logger, meta=self.meta.vnir)
        self.swir = EnMAP_Detector_SensorGeo('SWIR', root_dir, config=self.cfg, logger=self.logger, meta=self.meta.swir)

        # Get the paths according information delivered in the metadata
        self.paths = self.get_paths(root_dir)

        self.detector_attrNames = ['vnir', 'swir']

    def get_paths(self, root_dir):
        """
        Get all file paths associated with the current instance of EnMAPL1Product_SensorGeo
        :param root_dir: directory where the data are located
        :return: paths.SimpleNamespace()
        """
        if self.cfg.is_dlr_dataformat:
            paths = SimpleNamespace()
            paths.vnir = self.vnir.get_paths()
            paths.swir = self.swir.get_paths()
            paths.root_dir = root_dir
            paths.metaxml = glob(path.join(root_dir, "*METADATA.XML"))[0]
        else:
            paths = SimpleNamespace()
            paths.vnir = self.vnir.get_paths()
            paths.swir = self.swir.get_paths()
            paths.root_dir = root_dir
            paths.metaxml = glob(path.join(root_dir, "*_header.xml"))[0]
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

    def correct_dead_pixels(self):
        """Correct dead pixels of both detectors."""
        self.vnir.correct_dead_pixels()
        self.swir.correct_dead_pixels()

    def get_preprocessed_dem(self):
        self.vnir.get_preprocessed_dem()
        self.swir.get_preprocessed_dem()

    def run_AC(self):
        from ..processors.atmospheric_correction import AtmosphericCorrector
        AC = AtmosphericCorrector(config=self.cfg)
        AC.run_ac(self)

    # Define a new save to take into account the fact that 2 images might be appended
    # Here we save the radiance and not DN (otherwise there will be a problem with the concatened images)
    def save(self, outdir: str, suffix="") -> str:
        """
        Save the product to disk using almost the same input format
        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        product_dir = path.join(
            path.abspath(outdir), "{name}{suffix}".format(
                name=[ff for ff in self.paths.root_dir.split(path.sep) if ff != ''][-1],
                suffix=suffix)
        )
        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        # We can hardcode the detectors (?)
        # write the VNIR
        self.vnir.data.save(product_dir + path.sep + self.meta.vnir.data_filename, fmt="ENVI")
        self.vnir.mask_clouds.save(product_dir + path.sep + self.meta.vnir.cloud_mask_filename, fmt="GTiff")
        if self.vnir.deadpixelmap is not None:
            self.vnir.deadpixelmap.save(product_dir + path.sep + self.meta.vnir.dead_pixel_filename, fmt="GTiff")
        else:
            self.logger.warning('Could not save VNIR dead pixel map because there is no corresponding attribute.')

        # FIXME we could also write the quicklook included in DLR L1B format
        self.vnir.generate_quicklook() \
            .save(path.join(product_dir, path.basename(self.meta.vnir.quicklook_filename) + '.png'), fmt='PNG')

        # write the SWIR
        self.swir.data.save(product_dir + path.sep + self.meta.swir.data_filename, fmt="ENVI")
        self.swir.mask_clouds.save(product_dir + path.sep + self.meta.swir.cloud_mask_filename, fmt="GTiff")
        if self.swir.deadpixelmap is not None:
            self.swir.deadpixelmap.save(product_dir + path.sep + self.meta.swir.dead_pixel_filename, fmt="GTiff")
        else:
            self.logger.warning('Could not save SWIR dead pixel map because there is no corresponding attribute.')
        self.swir.generate_quicklook() \
            .save(path.join(product_dir, path.basename(self.meta.swir.quicklook_filename) + '.png'), fmt='PNG')

        # Update the xml file
        if self.cfg.is_dlr_dataformat:
            xml = ElementTree.parse(self.paths.metaxml).getroot()
            for detName, detMeta in zip(['VNIR', 'SWIR'], [self.meta.vnir, self.meta.swir]):
                lbl = L1B_product_props_DLR['xml_detector_label'][detName]
                xml.find("product/image/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
                xml.find("product/image/%s/dimension/columns" % lbl).text = str(detMeta.ncols)
                xml.find("product/quicklook/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
                xml.find("product/quicklook/%s/dimension/columns" % lbl).text = str(detMeta.ncols)

            with open(product_dir + path.sep + path.basename(self.paths.metaxml), "w") as xml_file:
                xml_file.write(ElementTree.tostring(xml, encoding='unicode'))
        else:
            xml = ElementTree.parse(self.paths.metaxml).getroot()
            for detName, detMeta in zip(['VNIR', 'SWIR'], [self.meta.vnir, self.meta.swir]):
                lbl = L1B_product_props['xml_detector_label'][detName]
                xml.find("ProductComponent/%s/Data/Size/NRows" % lbl).text = str(detMeta.nrows)
                xml.find("ProductComponent/%s/Data/Type/UnitCode" % lbl).text = detMeta.unitcode
                xml.find("ProductComponent/%s/Data/Type/Unit" % lbl).text = detMeta.unit
            with open(product_dir + path.sep + path.basename(self.paths.metaxml), "w") as xml_file:
                xml_file.write(ElementTree.tostring(xml, encoding='unicode'))

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

    def save(self, outdir: str, suffix="") -> str:
        """
        Save the product to disk using almost the same input format
        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        product_dir = path.join(
            path.abspath(outdir), "{name}{suffix}".format(
                name=[ff for ff in self.paths.root_dir.split(path.sep) if ff != ''][-1],
                suffix=suffix)
        )
        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        raise NotImplementedError()  # TODO implement save method for L2A data

        return product_dir
