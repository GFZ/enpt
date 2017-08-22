# -*- coding: utf-8 -*-

import logging
from types import SimpleNamespace

from geomultisens.model.dataset import Dataset

from ..utils.path_generator import PathGenL1BProduct
from ..model.metadata import EnMAP_Metadata_ImGeo, EnMAP_Metadata_VNIR_ImGeo, EnMAP_Metadata_SWIR_ImGeo


##############
# BASE CLASSES
##############


class _EnMAP_Image(Dataset):
    def __init__(self):
        """This is the base class for all kinds of EnMAP images."""

        # get all attributes of base class "Dataset"
        super(_EnMAP_Image, self).__init__()

        # add EnMAP specific attributes
        self.paths = SimpleNamespace()


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

    def DN2Radiance(self):
        self.logger.info('Converting DN values to radiance for %s detector...' % self.detector_name)
        self.arr = (self.meta.l_min + (self.meta.l_max - self.meta.l_min) / (2 ** 16 - 1) * self.arr[:])
        self.meta.unit = "radiance"


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

    def DN2Radiance(self):
        self.vnir.DN2Radiance()
        self.swir.DN2Radiance()


####################################
# EnPT EnMAP objects in map geometry
####################################


class EnMAP_VNIR_MapGeo(_EnMAP_Detector_MapGeo):
    def __init__(self, logger=None):
        super(EnMAP_VNIR_MapGeo, self).__init__('VNIR', logger=logger)


class EnMAP_SWIR_MapGeo(_EnMAP_Detector_MapGeo):
    def __init__(self, logger=None):
        super(EnMAP_SWIR_MapGeo, self).__init__('SWIR', logger=logger)


####################
# DEPRECATED CLASSES
####################


class EnMAP_L1B(_EnMAP_Image):
    """This class represents an EnMAP L1B image including all metadata and associated aux-data (masks, DEM, etc.).

    All attributes commonly used among different EnMAP images are inherited from the _EnMAP_Image class.
    L1B specific modifications are to be implemented here."""

    pass


class EnMAP_L2A(_EnMAP_Image):
    """This class represents an EnMAP L2A image including all metadata and associated aux-data (masks, DEM, etc.).

    All attributes commonly used among different EnMAP images are inherited from the _EnMAP_Image class.
    L2A specific modifications are to be implemented here."""

    pass
