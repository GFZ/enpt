# -*- coding: utf-8 -*-
"""EnPT images module. All EnMAP image objects are defined here."""

import logging
from types import SimpleNamespace
import numpy as np
from os import path, sep, makedirs
from shutil import copyfile
from xml.etree import ElementTree

from geoarray import GeoArray, NoDataMask, CloudMask

from ..utils.path_generator import PathGenL1BProduct
from ..utils.logging import EnPT_Logger
from ..model.metadata import EnMAP_Metadata_L1B_SensorGeo, EnMAP_Metadata_L1B_Detector_SensorGeo
from ..options.config import EnPTConfig

##############
# BASE CLASSES
##############


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
        # get all attributes of base class "Dataset"
        super(_EnMAP_Image, self).__init__()

        # add EnMAP specific attributes
        self.paths = SimpleNamespace()

        # protected attributes
        self._logger = None
        self._log = ''
        self._arr = None
        self._mask_nodata = None
        self._mask_clouds = None
        self._mask_clouds_confidence = None
        self._dem = None
        self._deadpixelmap = None
        self._ac_options = {}
        self._ac_errors = None

        # defaults
        self.entity_ID = ''
        self.basename = ''

    @property
    def logger(self):
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
            self._logger = EnPT_Logger('log__' + self.basename, fmt_suffix=None, path_logfile='',
                                       log_level='INFO', append=True)  # TODO revise the logger

            return self._logger

    @logger.setter
    def logger(self, logger: logging.Logger):
        assert isinstance(logger, logging.Logger) or logger in ['not set', None], \
            "_EnMAP_Image.logger can not be set to %s." % logger

        # save prior logs
        # if logger is None and self._logger is not None:
        #     self.log += self.logger.captured_stream
        self._logger = logger

    @property  # FIXME does not work yet
    def log(self):
        """Return a string of all logged messages until now.

        NOTE: self.log can also be set to a string.
        """
        return self._log

    @log.setter
    def log(self, string):
        assert isinstance(string, str), "'log' can only be set to a string. Got %s." % type(string)
        self._log = string

    @property
    def data(self):
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
        # TODO this must return a subset if self.subset is not None
        return self._arr

    @data.setter
    def data(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            # TODO this must be able to handle subset inputs in tiled processing
            if self._arr and len(geoArr_initArgs[0]) and isinstance(geoArr_initArgs[0], np.ndarray):
                self._arr = GeoArray(geoArr_initArgs[0], geotransform=self._arr.gt, projection=self._arr.prj)
            else:
                self._arr = GeoArray(*geoArr_initArgs)
        else:
            del self.data

    @data.deleter
    def data(self):
        self._arr = None

    @property
    def mask_clouds(self):
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
    def mask_clouds_confidence(self):
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
    def dem(self):
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
    def ac_errors(self):
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
    def deadpixelmap(self):
        """Return the dead pixel map.

        Bundled with all the corresponding metadata. Dimensions: (bands x columns).

        For usage instructions and a list of attributes refer to help(self.data).
        self.mask_clouds_confidence works in the same way.

        :return: instance of geoarray.GeoArray
        """
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
    """Class representing a single detector of an EnMAP image (as sensor geometry).

    NOTE:
        - Inherits all attributes from _EnMAP_Image class.
        - All functionality that VNIR and SWIR detectors (sensor geometry) have in common is to be implemented here.

    Attributes:
        - to be listed here. Check help(_EnMAP_Detector_SensorGeo) in the meanwhile!

    """

    def __init__(self, detector_name: str, root_dir: str, config: EnPTConfig, logger=None):
        """Get an instance of _EnMAP_Detector_SensorGeo.

        :param detector_name:   'VNIR' or 'SWIR'
        :param root_dir:
        :param logger:
        """
        if detector_name not in ['VNIR', 'SWIR']:
            raise ValueError("'detector_name' must be 'VNIR' or 'SWIR'. Received %s!" % detector_name)

        self.cfg = config
        self._root_dir = root_dir
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(EnMAP_Detector_SensorGeo, self).__init__()
        self.paths = self.get_paths()
        # instance an empty metadata object
        self.detector_meta = EnMAP_Metadata_L1B_Detector_SensorGeo(self.detector_name, config=self.cfg, logger=logger)

    def get_paths(self):
        """Get all file paths associated with the current instance of _EnMAP_Detector_SensorGeo.

        :return: types.SimpleNamespace
        """
        pathGen = PathGenL1BProduct(self._root_dir, self.detector_name)
        paths = SimpleNamespace()
        paths.root_dir = self._root_dir
        paths.metaxml = pathGen.get_path_metaxml()
        paths.data = pathGen.get_path_data()
        paths.mask_clouds = pathGen.get_path_cloudmask()
        paths.deadpixelmap = pathGen.get_path_deadpixelmap()
        paths.quicklook = pathGen.get_path_quicklook()

        return paths

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
            self.logger.info(
                "No is DN to Radiance conversion is performed since unitcode is not DN (found: {code}).".format(
                    code=self.detector_meta.unitcode)
            )


class EnMAPL1Product_SensorGeo(object):
    """Class for EnPT EnMAP object in sensor geometry.

    Attributes:
        - logger:
            - logging.Logger instance or subclassed
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

        :param root_dir: Root directory of EnMAP Level-1B product
        :param logger: None or logging instance
        """
        self.cfg = config
        self.logger = logger or logging.getLogger(__name__)
        self.vnir = EnMAP_Detector_SensorGeo('VNIR', root_dir, config=self.cfg, logger=logger)
        self.swir = EnMAP_Detector_SensorGeo('SWIR', root_dir, config=self.cfg, logger=logger)
        self.paths = self.get_paths()
        self.meta = EnMAP_Metadata_L1B_SensorGeo(self.paths.metaxml, config=self.cfg, logger=logger)
        self.detector_attrNames = ['vnir', 'swir']

    def get_paths(self):
        """Get all file paths associated with the current instance of EnMAPL1Product_SensorGeo.

        :return: types.SimpleNamespace()
        """
        paths = SimpleNamespace()
        paths.vnir = self.vnir.get_paths()
        paths.swir = self.swir.get_paths()
        paths.root_dir = paths.vnir.root_dir
        paths.metaxml = paths.vnir.metaxml

        return paths

    def DN2TOARadiance(self):
        """Convert self.vnir.data and self.swir.data from digital numbers to top-of-atmosphere radiance.

        :return: None
        """
        self.vnir.DN2TOARadiance()
        self.swir.DN2TOARadiance()

    def save(self, outdir: str, suffix="") -> str:
        """Save this product to disk using almost the same format as for reading.

        :param outdir: Path to output directory
        :return: Root path of written product
        """
        product_dir = path.join(
            path.abspath(outdir),
            "{name}{suffix}".format(
                name=[ff for ff in self.paths.root_dir.split(sep) if ff != ''][-1],
                suffix=suffix)
        )
        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        for detector_name in self.detector_attrNames:
            detector = getattr(self, detector_name)
            detector_paths = getattr(self.paths, detector_name)

            for atts, fmt in ((("deadpixelmap", "mask_clouds"), "GTIff"),
                              (("data",), "ENVI")):
                for att in atts:
                    getattr(detector, att).save(
                        path.join(product_dir, path.basename(getattr(detector_paths, att))), fmt=fmt)

            copyfile(
                src=detector_paths.quicklook,
                dst=path.join(product_dir, path.basename(detector_paths.quicklook))
            )

        xml = ElementTree.parse(self.paths.metaxml)
        for xml_name, real_name in (("detector1", "vnir"), ("detector2", "swir")):
            ele = xml.getroot().find(xml_name)
            new_ele = ElementTree.Element("unitcode")
            new_ele.text = getattr(self.meta, real_name).unitcode
            ele.append(new_ele)
        xml.write(path.join(product_dir, path.basename(self.paths.metaxml)))

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
    def mask_nodata(self):
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

    def calc_mask_nodata(self, fromBand=None, overwrite=False):
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
