# -*- coding: utf-8 -*-

import logging
from types import SimpleNamespace
import numpy as np

from geoarray import GeoArray, NoDataMask, CloudMask

from ..utils.path_generator import PathGenL1BProduct
from ..model.metadata import EnMAP_Metadata_ImGeo, EnMAP_Metadata_VNIR_ImGeo, EnMAP_Metadata_SWIR_ImGeo


##############
# BASE CLASSES
##############


class _EnMAP_Image(object):
    def __init__(self):
        """This is the base class for all kinds of EnMAP images."""

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
        if self._logger and self._logger.handlers[:]:
            return self._logger
        else:
            self._logger = logging.getLogger()
            # DatasetLogger('log__' + self.baseN, fmt_suffix=self.scene_ID, path_logfile=self.path_logfile,
            #                          log_level='INFO', append=True)
            return self._logger

    @logger.setter
    def logger(self, logger):
        assert isinstance(logger, logging.Logger) or logger in ['not set', None], \
            "GMS_obj.logger can not be set to %s." % logger

        # save prior logs
        # if logger is None and self._logger is not None:
        #     self.log += self.logger.captured_stream
        self._logger = logger

    @property  # FIXME does not work yet
    def log(self):
        """Returns a string of all logged messages until now."""

        return self._log

    @log.setter
    def log(self, string):
        assert isinstance(string, str), "'log' can only be set to a string. Got %s." % type(string)
        self._log = string

    @property
    def arr(self):
        # TODO this must return a subset if self.subset is not None
        return self._arr

    @arr.setter
    def arr(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            # TODO this must be able to handle subset inputs in tiled processing
            if self._arr and len(geoArr_initArgs[0]) and isinstance(geoArr_initArgs[0], np.ndarray):
                self._arr = GeoArray(geoArr_initArgs[0], geotransform=self._arr.gt, projection=self._arr.prj)
            else:
                self._arr = GeoArray(*geoArr_initArgs)
        else:
            del self.arr

    @arr.deleter
    def arr(self):
        self._arr = None

    @property
    def mask_nodata(self):
        if self._mask_nodata is None and isinstance(self.arr, GeoArray):
            self.logger.info('Calculating nodata mask...')
            self._mask_nodata = self.arr.mask_nodata  # calculates mask nodata if not already present

        return self._mask_nodata

    @mask_nodata.setter
    def mask_nodata(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            nd = NoDataMask(*geoArr_initArgs)
            if nd.shape[:2] != self.arr.shape[:2]:
                raise ValueError("The 'mask_nodata' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(nd.shape))
            nd.nodata = False
            nd.gt = self.arr.gt
            nd.prj = self.arr.prj
            self._mask_nodata.prj = nd
        else:
            del self.mask_nodata

    @mask_nodata.deleter
    def mask_nodata(self):
        self._mask_nodata = None

    @property
    def mask_clouds(self):
        return self._mask_clouds

    @mask_clouds.setter
    def mask_clouds(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            cm = CloudMask(*geoArr_initArgs)
            if cm.shape[:2] != self.arr.shape[:2]:
                raise ValueError("The 'mask_clouds' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(cm.shape))
            cm.nodata = 0
            cm.gt = self.arr.gt
            cm.prj = self.arr.prj
            self._mask_clouds = cm
        else:
            del self.mask_clouds

    @mask_clouds.deleter
    def mask_clouds(self):
        self._mask_clouds = None

    @property
    def mask_clouds_confidence(self):
        return self._mask_clouds_confidence

    @mask_clouds_confidence.setter
    def mask_clouds_confidence(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            cnfArr = GeoArray(*geoArr_initArgs)

            if not cnfArr.shape == self.arr.shape[:2]:
                raise ValueError("The 'mask_clouds_confidence' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(cnfArr.shape))

            if cnfArr._nodata is None:
                cnfArr.nodata = 0  # DEF_D.get_outFillZeroSaturated(cnfArr.dtype)[0] # TODO
            cnfArr.gt = self.arr.gt
            cnfArr.prj = self.arr.prj

            self._mask_clouds_confidence = cnfArr
        else:
            del self._mask_clouds_confidence

    @mask_clouds_confidence.deleter
    def mask_clouds_confidence(self):
        self._mask_clouds_confidence = None

    @property
    def dem(self):
        """
        Returns an SRTM DEM in the exact dimension an pixel grid of self.arr as an instance of GeoArray.
        """

        if self._dem is None:
            raise NotImplementedError('An automatic DEM getter is not yet implemented.')
        return self._dem

    @dem.setter
    def dem(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dem = GeoArray(*geoArr_initArgs)
            assert self._dem.shape[:2] == self.arr.shape[:2]

            self._dem = dem
            if not dem.shape == self.arr.shape[:2]:
                raise ValueError("The 'dem' GeoArray can only be instanced with an array of the "
                                 "same dimensions like _EnMAP_Image.arr. Got %s." % str(dem.shape))
            self._dem.nodata = 0  # FIXME
            self._dem.gt = self.arr.gt
            self._dem.prj = self.arr.prj
        else:
            del self.dem

    @dem.deleter
    def dem(self):
        self._dem = None

    @property
    def ac_errors(self):
        """Returns an instance of GeoArray containing error information calculated by the atmospheric correction.

        :return:
        """

        return self._ac_errors  # FIXME should give a warning if None

    @ac_errors.setter
    def ac_errors(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            errArr = GeoArray(*geoArr_initArgs)

            if errArr.shape != self.arr.shape:
                raise ValueError("The 'ac_errors' GeoArray can only be instanced with an array of "
                                 "the same dimensions like _EnMAP_Image.arr. Got %s." % str(errArr.shape))

            if errArr._nodata is None:
                errArr.nodata = 0  # DEF_D.get_outFillZeroSaturated(errArr.dtype)[0] # TODO
            errArr.gt = self.arr.gt
            errArr.prj = self.arr.prj
            # errArr.bandnames = self.LBA2bandnames(self.LayerBandsAssignment)

            self._ac_errors = errArr
        else:
            del self.ac_errors

    @ac_errors.deleter
    def ac_errors(self):
        self._ac_errors = None

    @property
    def deadpixelmap(self):
        return self._deadpixelmap

    @deadpixelmap.setter
    def deadpixelmap(self, *geoArr_initArgs):
        if geoArr_initArgs[0] is not None:
            dpm = GeoArray(*geoArr_initArgs)

            if dpm.shape != (self.arr.bands, self.arr.cols):
                raise ValueError("The 'deadpixelmap' GeoArray can only be instanced with an array with the size "
                                 "'bands x columns' of the GeoArray _EnMAP_Image.arr. Got %s." % str(dpm.shape))

            self._deadpixelmap = dpm
        else:
            del self._deadpixelmap

    @deadpixelmap.deleter
    def deadpixelmap(self):
        self._deadpixelmap = None

    def calc_mask_nodata(self, fromBand=None, overwrite=False):
        """Calculates a no data mask with (values: 0=nodata; 1=data)

        :param fromBand:  <int> index of the band to be used (if None, all bands are used)
        :param overwrite: <bool> whether to overwrite existing nodata mask that has already been calculated
        :return:
        """

        self.logger.info('Calculating nodata mask...')

        if self._mask_nodata is None or overwrite:
            self.arr.calc_mask_nodata(fromBand=fromBand, overwrite=overwrite)
            self.mask_nodata = self.arr.mask_nodata
            return self.mask_nodata


class _EnMAP_Detector_ImGeo(_EnMAP_Image):
    def __init__(self, detector_name: str, root_dir: str, logger=None):
        """

        :param detector_name:   VNIR or SWIR
        :param root_dir:
        :param logger:
        """
        self._root_dir = root_dir
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(_EnMAP_Detector_ImGeo, self).__init__()
        self.paths = self.get_paths()
        # instance an empty metadata object
        self.meta = \
            EnMAP_Metadata_VNIR_ImGeo(self.paths.metaxml, logger=logger) if self.detector_name == 'VNIR' else \
            EnMAP_Metadata_SWIR_ImGeo(self.paths.metaxml, logger=logger)

    def get_paths(self):
        pathGen = PathGenL1BProduct(self._root_dir, self.detector_name)
        paths = SimpleNamespace()
        paths.root_dir = self._root_dir
        paths.metaxml = pathGen.get_path_metaxml()
        paths.imagedata = pathGen.get_path_imagedata()
        paths.mask_clouds = pathGen.get_path_cloudmask()
        paths.deadpixelmap = pathGen.get_path_deadpixelmap()
        paths.quicklook = pathGen.get_path_quicklook()

        return paths

    def DN2TOARadiance(self):
        # TODO move to processors.radiometric_transform?
        self.logger.info('Converting DN values to radiance for %s detector...' % self.detector_name)
        self.arr = (self.meta.l_min + (self.meta.l_max - self.meta.l_min) / (2 ** 16 - 1) * self.arr[:])
        self.meta.unit = "mW m^-2 sr^-1 nm^-1"
        self.meta.unitcode = "TOARad"


class _EnMAP_Detector_MapGeo(_EnMAP_Image):
    def __init__(self, detector_name: str, logger=None):
        """

        :param detector_name:   VNIR or SWIR
        :param logger:
        """
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(_EnMAP_Detector_MapGeo, self).__init__()


######################################
# EnPT EnMAP objects in image geometry
######################################


class EnMAP_VNIR_ImGeo(_EnMAP_Detector_ImGeo):
    def __init__(self, root_dir: str):
        """Get an instance of the VNIR of an EnMAP data Level-1B product.

        :param root_dir: Root directory of EnMAP Level-1B product
        """
        # get all attributes of base class "_EnMAP_Detector"
        super(EnMAP_VNIR_ImGeo, self).__init__('VNIR', root_dir)


class EnMAP_SWIR_ImGeo(_EnMAP_Detector_ImGeo):
    def __init__(self, root_dir: str):
        """Get an instance of the SWIR of an EnMAP data Level-1B product.

        :param root_dir: Root directory of EnMAP Level-1B product
        """
        # get all attributes of base class "_EnMAP_Detector"
        super(EnMAP_SWIR_ImGeo, self).__init__('SWIR', root_dir)


class EnMAPL1Product_ImGeo(object):
    """Class for EnPT EnMAP object in image geometry

    Attributes:
        - vnir
            - ...
        - swir
            - same as vor [vnir]
        - paths: paths belonging to the EnMAP product

    """
    def __init__(self, root_dir: str, logger=None):
        """Get instance of EnPT EnMAP object in image geometry.

        :param root_dir: Root directory of EnMAP Level-1B product
        :param logger: None or logging instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self.vnir = EnMAP_VNIR_ImGeo(root_dir)
        self.swir = EnMAP_SWIR_ImGeo(root_dir)
        self.paths = self.get_paths()
        self.meta = EnMAP_Metadata_ImGeo(self.paths.metaxml, logger=logger)

    def get_paths(self):
        """Get all the paths belonging to the EnMAP L1B product

        :return:
        """
        paths = SimpleNamespace()
        paths.vnir = self.vnir.get_paths()
        paths.swir = self.swir.get_paths()
        paths.root_dir = paths.vnir.root_dir
        paths.metaxml = paths.vnir.metaxml

        return paths

    def DN2TOARadiance(self):
        self.vnir.DN2TOARadiance()
        self.swir.DN2TOARadiance()


####################################
# EnPT EnMAP objects in map geometry
####################################


class EnMAP_VNIR_MapGeo(_EnMAP_Detector_MapGeo):
    """This class represents an EnPT EnMAP VNIR image in map geometry including all metadata and associated aux-data
    (masks, DEM, etc.).

    All attributes commonly used among different EnMAP images are inherited from the _EnMAP_Detector_MapGeo class.
    All VNIR_MapGeo specific modifications are to be implemented here."""
    def __init__(self, logger=None):
        super(EnMAP_VNIR_MapGeo, self).__init__('VNIR', logger=logger)


class EnMAP_SWIR_MapGeo(_EnMAP_Detector_MapGeo):
    """This class represents an EnPT EnMAP SWIR image in map geometry including all metadata and associated aux-data
    (masks, DEM, etc.).

    All attributes commonly used among different EnMAP images are inherited from the _EnMAP_Detector_MapGeo class.
    All SWIR_MapGeo specific modifications are to be implemented here."""
    def __init__(self, logger=None):
        super(EnMAP_SWIR_MapGeo, self).__init__('SWIR', logger=logger)
