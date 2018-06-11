# -*- coding: utf-8 -*-
"""EnPT images module. All EnMAP image objects are defined here."""

import logging
from types import SimpleNamespace
import numpy as np
from os import path, makedirs
from lxml import etree
from glob import glob
import utm
from scipy.interpolate import interp2d

# Use to generate preview
import imageio
from skimage import exposure

from geoarray import GeoArray, NoDataMask, CloudMask

from ..utils.logging import EnPT_Logger
from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo, EnMAP_Metadata_L1B_Detector_SensorGeo
from ..options.config import EnPTConfig
from ..processors.dead_pixel_correction import Dead_Pixel_Corrector


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
            del self._mask_clouds_confidence

    @mask_clouds_confidence.deleter
    def mask_clouds_confidence(self):
        self._mask_clouds_confidence = None

    @property
    def dem(self) -> GeoArray:
        """Return an SRTM DEM in the exact dimension an pixel grid of self.arr.

        :return: geoarray.GeoArray
        """
        if self._dem is None:
            raise NotImplementedError('An automatic DEM getter is not yet implemented.')
        return self._dem

    @dem.setter
    def dem(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dem = GeoArray(*geoArr_initArgs)
            assert self._dem.shape[:2] == self.data.shape[:2]

            self._dem = dem
            if not dem.shape == self.data.shape[:2]:
                raise ValueError("The 'dem' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(dem.shape))
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

            if dpm.shape != (self.data.bands, self.data.cols):
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with an array with the size "
                                 "'bands x columns' of the GeoArray _EnMAP_Image.arr. Got %s." % str(dpm.shape))

            self._deadpixelmap = dpm
        else:
            del self._deadpixelmap

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
            self.detector_meta = meta

    def get_paths(self):
        """
        Get all file paths associated with the current instance of EnMAP_Detector_SensorGeo
        These information are reading from the detector_meta
        :return: paths as SimpleNamespace
        """
        paths = SimpleNamespace()
        paths.root_dir = self._root_dir
        paths.metaxml = glob(path.join(self._root_dir, "*_header.xml"))[0]
        paths.data = path.join(self._root_dir, self.detector_meta.data_filename)
        paths.mask_clouds = path.join(self._root_dir, self.detector_meta.cloud_mask_filename)
        paths.deadpixelmap = path.join(self._root_dir, self.detector_meta.dead_pixel_filename)
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

    def DN2TOARadiance(self):
        """Convert DNs to TOA radiance.

        Convert the radiometric unit of _EnMAP_Detector_SensorGeo.data from digital numbers to top-of-atmosphere
        radiance.

        :return: None
        """
        # TODO move to processors.radiometric_transform?
        if self.detector_meta.unitcode == 'DN':
            self.logger.info('Converting DN values to radiance for %s detector...' % self.detector_name)
            self.data = (self.detector_meta.l_min + (self.detector_meta.l_max - self.detector_meta.l_min) /
                         (2 ** 16 - 1) * self.data[:])
            self.detector_meta.unit = "mW m^-2 sr^-1 nm^-1"
            self.detector_meta.unitcode = "TOARad"
        else:
            self.logger.warning(
                "No DN to Radiance conversion is performed because unitcode is not DN (found: {code}).".format(
                    code=self.detector_meta.unitcode)
            )

    def generate_quicklook(self, filename):
        """
        Generate image quicklook and save into a file
        :param filename: file path to store the image
        :return: None
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
        imageio.imwrite(filename, pix)


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
        self.vnir.deadpixelmap.save(product_dir + path.sep + self.meta.vnir.dead_pixel_filename, fmt="GTiff")
        self.vnir.generate_quicklook(product_dir + path.sep + self.meta.vnir.quicklook_filename)

        # write the SWIR
        self.swir.data.save(product_dir + path.sep + self.meta.swir.data_filename, fmt="ENVI")
        self.swir.mask_clouds.save(product_dir + path.sep + self.meta.swir.cloud_mask_filename, fmt="GTiff")
        self.swir.deadpixelmap.save(product_dir + path.sep + self.meta.swir.dead_pixel_filename, fmt="GTiff")
        self.swir.generate_quicklook(product_dir + path.sep + self.meta.swir.quicklook_filename)

        # Update the xml file
        xml = etree.parse(self.paths.metaxml)
        xml.findall("ProductComponent/VNIRDetector/Data/Size/NRows")[0].text = str(self.meta.vnir.nrows)
        xml.findall("ProductComponent/VNIRDetector/Data/Type/UnitCode")[0].text = self.meta.vnir.unitcode
        xml.findall("ProductComponent/VNIRDetector/Data/Type/Unit")[0].text = self.meta.vnir.unit
        xml.findall("ProductComponent/SWIRDetector/Data/Size/NRows")[0].text = str(self.meta.swir.nrows)
        xml.findall("ProductComponent/SWIRDetector/Data/Type/UnitCode")[0].text = self.meta.swir.unitcode
        xml.findall("ProductComponent/SWIRDetector/Data/Type/Unit")[0].text = self.meta.swir.unit
        xml_string = etree.tostring(xml, pretty_print=True, xml_declaration=True, encoding='UTF-8')
        with open(product_dir + path.sep + path.basename(self.paths.metaxml), "w") as xml_file:
            xml_file.write(xml_string.decode("utf-8"))

        return product_dir

    def append_new_image(self, img2, n_lines: int=None):
        """
        Check if a second image could pass with the first image.
        In this version we assume that the image will be add below
        If it is the case, the method will create temporary files that will be used in the following.
        :param img2:
        :param n_lines: number of line to be added
        :return: None
        """
        self.logger.info("Check new image %s" % img2.paths.root_dir)

        distance_min = 27.0
        distance_max = 34.0
        tag_vnir = False
        tag_swir = False

        # check vnir bottom left
        x1, y1, _, _ = utm.from_latlon(self.meta.vnir.lat_UL_UR_LL_LR[2], self.meta.vnir.lon_UL_UR_LL_LR[2])
        x2, y2, _, _ = utm.from_latlon(img2.meta.vnir.lat_UL_UR_LL_LR[0], img2.meta.vnir.lon_UL_UR_LL_LR[0])
        distance_left = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        # check vnir bottom right
        x1, y1, _, _ = utm.from_latlon(self.meta.vnir.lat_UL_UR_LL_LR[3], self.meta.vnir.lon_UL_UR_LL_LR[3])
        x2, y2, _, _ = utm.from_latlon(img2.meta.vnir.lat_UL_UR_LL_LR[1], img2.meta.vnir.lon_UL_UR_LL_LR[1])
        distance_right = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        if distance_min < distance_left < distance_max and distance_min < distance_right < distance_max:
            tag_vnir = True

        # check swir bottom left
        x1, y1, _, _ = utm.from_latlon(self.meta.swir.lat_UL_UR_LL_LR[2], self.meta.swir.lon_UL_UR_LL_LR[2])
        x2, y2, _, _ = utm.from_latlon(img2.meta.swir.lat_UL_UR_LL_LR[0], img2.meta.swir.lon_UL_UR_LL_LR[0])
        distance_left = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        # check swir bottom right
        x1, y1, _, _ = utm.from_latlon(self.meta.swir.lat_UL_UR_LL_LR[3], self.meta.swir.lon_UL_UR_LL_LR[3])
        x2, y2, _, _ = utm.from_latlon(img2.meta.swir.lat_UL_UR_LL_LR[1], img2.meta.swir.lon_UL_UR_LL_LR[1])
        distance_right = np.sqrt((x1-x2)**2 + (y1-y2)**2)

        if distance_min < distance_left < distance_max and distance_min < distance_right < distance_max:
            tag_swir = True

        if tag_vnir is False or tag_swir is False:
            self.logger.warning("%s and %s don't fit to be appended" % (self.paths.root_dir, img2.paths.root_dir))
            return

        self.logger.info("Append new image: %s" % img2.paths.root_dir)

        # set new number of line
        if n_lines is None:
            n_lines = img2.meta.vnir.nrows

        if n_lines > img2.meta.vnir.nrows:
            self.logger.warning("n_lines (%s) exceeds the total number of line of second image" % n_lines)
            self.logger.warning("Set to the image number of line")
            n_lines = img2.meta.vnir.nrows

        if n_lines < 50:  # TODO: determine these values
            self.logger.warning("A minimum of 50 lines is required, only %s were selected" % n_lines)
            self.logger.warning("Set the number of line to 50")
            n_lines = 50

        self.meta.vnir.nrows += n_lines

        # Create new corner coordinate - VNIR
        ff = interp2d(x=[0, 1], y=[0, 1], z=[[img2.meta.vnir.lat_UL_UR_LL_LR[0], img2.meta.vnir.lat_UL_UR_LL_LR[1]],
                                             [img2.meta.vnir.lat_UL_UR_LL_LR[2], img2.meta.vnir.lat_UL_UR_LL_LR[3]]],
                      kind='linear')
        self.meta.vnir.lat_UL_UR_LL_LR[2] = np.array(ff(0, n_lines/img2.meta.vnir.nrows))[0]
        self.meta.vnir.lat_UL_UR_LL_LR[3] = np.array(ff(1, n_lines/img2.meta.vnir.nrows))[0]
        lon_lat_smpl = (15, 15)
        self.meta.vnir.lats = self.meta.vnir.interpolate_corners(*self.meta.vnir.lat_UL_UR_LL_LR, *lon_lat_smpl)

        ff = interp2d(x=[0, 1], y=[0, 1], z=[[img2.meta.vnir.lon_UL_UR_LL_LR[0], img2.meta.vnir.lon_UL_UR_LL_LR[1]],
                                             [img2.meta.vnir.lon_UL_UR_LL_LR[2], img2.meta.vnir.lon_UL_UR_LL_LR[3]]],
                      kind='linear')
        self.meta.vnir.lon_UL_UR_LL_LR[2] = np.array(ff(0, n_lines/img2.meta.vnir.nrows))[0]
        self.meta.vnir.lon_UL_UR_LL_LR[3] = np.array(ff(1, n_lines/img2.meta.vnir.nrows))[0]
        self.meta.vnir.lons = self.meta.vnir.interpolate_corners(*self.meta.vnir.lon_UL_UR_LL_LR, *lon_lat_smpl)

        # Create new corner coordinate - SWIR
        ff = interp2d(x=[0, 1], y=[0, 1], z=[[img2.meta.swir.lat_UL_UR_LL_LR[0], img2.meta.swir.lat_UL_UR_LL_LR[1]],
                                             [img2.meta.swir.lat_UL_UR_LL_LR[2], img2.meta.swir.lat_UL_UR_LL_LR[3]]],
                      kind='linear')
        self.meta.swir.lat_UL_UR_LL_LR[2] = np.array(ff(0, n_lines/img2.meta.swir.nrows))[0]
        self.meta.swir.lat_UL_UR_LL_LR[3] = np.array(ff(1, n_lines/img2.meta.swir.nrows))[0]
        lon_lat_smpl = (15, 15)
        self.meta.swir.lats = self.meta.swir.interpolate_corners(*self.meta.swir.lat_UL_UR_LL_LR, *lon_lat_smpl)
        ff = interp2d(x=[0, 1], y=[0, 1], z=[[img2.meta.vnir.lon_UL_UR_LL_LR[0], img2.meta.vnir.lon_UL_UR_LL_LR[1]],
                                             [img2.meta.vnir.lon_UL_UR_LL_LR[2], img2.meta.vnir.lon_UL_UR_LL_LR[3]]],
                      kind='linear')
        self.meta.swir.lon_UL_UR_LL_LR[2] = np.array(ff(0, n_lines/img2.meta.swir.nrows))[0]
        self.meta.swir.lon_UL_UR_LL_LR[3] = np.array(ff(1, n_lines/img2.meta.swir.nrows))[0]
        self.meta.swir.lons = self.meta.swir.interpolate_corners(*self.meta.swir.lon_UL_UR_LL_LR, *lon_lat_smpl)

        # append the vnir/swir image
        img2.vnir.data = img2.paths.vnir.data
        img2.vnir.data = img2.vnir.data[0:n_lines, :, :]
        img2.swir.data = img2.paths.swir.data
        img2.swir.data = img2.swir.data[0:n_lines, :, :]
        img2.DN2TOARadiance()
        self.vnir.data = np.append(self.vnir.data, img2.vnir.data, axis=0)
        self.swir.data = np.append(self.swir.data, img2.swir.data, axis=0)

        # append the mask cloud
        self.vnir.mask_clouds = np.append(self.vnir.mask_clouds,
                                          GeoArray(img2.paths.vnir.mask_clouds)[0:n_lines, :], axis=0)
        self.swir.mask_clouds = np.append(self.swir.mask_clouds,
                                          GeoArray(img2.paths.swir.mask_clouds)[0:n_lines, :], axis=0)


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
