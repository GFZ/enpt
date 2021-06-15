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

"""EnPT EnMAP objects in map geometry."""

import logging
from types import SimpleNamespace
from typing import Tuple, Optional  # noqa: F401
from os import path, makedirs

from geoarray import GeoArray, NoDataMask

from .image_baseclasses import _EnMAP_Image
from .images_sensorgeo import EnMAPL1Product_SensorGeo
from ...utils.logging import EnPT_Logger
from ...model.metadata import EnMAP_Metadata_L2A_MapGeo  # noqa: F401  # only used for type hint
from ...options.config import EnPTConfig

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


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

        # private attributes
        self._mask_nodata = None

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
        self._mask_nodata = self._get_geoarray_with_datalike_geometry(geoArr_initArgs, 'mask_nodata',
                                                                      nodataVal=False, specialclass=NoDataMask)

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

        self.meta: Optional[EnMAP_Metadata_L2A_MapGeo] = None
        self.paths: Optional[SimpleNamespace] = None

        super(EnMAPL2Product_MapGeo, self).__init__()

    @property
    def logger(self) -> EnPT_Logger:
        """Get an instance of enpt.utils.logging.EnPT_Logger.

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
        from ...processors.orthorectification import Orthorectifier
        L2_obj = Orthorectifier(config=config).run_transformation(enmap_ImageL1=enmap_ImageL1)

        return L2_obj

    def get_paths(self, l2a_outdir: str):
        """Get all file paths associated with the current instance of EnMAP_Detector_SensorGeo.

        NOTE: This information is read from the detector_meta.

        :param l2a_outdir:  output directory of EnMAP Level-2A dataset
        :return: paths as SimpleNamespace
        """
        paths = SimpleNamespace()
        paths.root_dir = l2a_outdir
        paths.metaxml = path.join(l2a_outdir, self.meta.filename_metaxml)
        paths.data = path.join(l2a_outdir, self.meta.filename_data)
        paths.mask_landwater = path.join(l2a_outdir, self.meta.filename_mask_landwater)
        paths.mask_clouds = path.join(l2a_outdir, self.meta.filename_mask_clouds)
        paths.mask_cloudshadow = path.join(l2a_outdir, self.meta.filename_mask_cloudshadow)
        paths.mask_haze = path.join(l2a_outdir, self.meta.filename_mask_haze)
        paths.mask_snow = path.join(l2a_outdir, self.meta.filename_mask_snow)
        paths.mask_cirrus = path.join(l2a_outdir, self.meta.filename_mask_cirrus)
        paths.deadpixelmap_vnir = path.join(l2a_outdir, self.meta.filename_deadpixelmap_vnir)
        paths.deadpixelmap_swir = path.join(l2a_outdir, self.meta.filename_deadpixelmap_swir)
        paths.quicklook_vnir = path.join(l2a_outdir, self.meta.filename_quicklook_vnir)
        paths.quicklook_swir = path.join(l2a_outdir, self.meta.filename_quicklook_swir)

        return paths

    def save(self, outdir: str, suffix="") -> str:
        """Save the product to disk using almost the same input format.

        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        # TODO optionally add more output formats
        product_dir = path.join(path.abspath(outdir),
                                "{name}{suffix}".format(name=self.meta.scene_basename, suffix=suffix))

        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        # save raster data
        kwargs_save = \
            dict(fmt='GTiff',
                 creationOptions=["COMPRESS=LZW",
                                  "NUM_THREADS=%d" % self.cfg.CPUs,
                                  "INTERLEAVE=%s" % ('BAND' if self.cfg.output_interleave == 'band' else 'PIXEL')]
                 ) if self.cfg.output_format == 'GTiff' else \
            dict(fmt='ENVI',
                 creationOptions=["INTERLEAVE=%s" % ("BSQ" if self.cfg.output_interleave == 'band' else
                                                     "BIL" if self.cfg.output_interleave == 'line' else
                                                     "BIP")])
        outpaths = dict(metaxml=path.join(product_dir, self.meta.filename_metaxml))

        for attrName in ['data', 'mask_landwater', 'mask_clouds', 'mask_cloudshadow', 'mask_haze', 'mask_snow',
                         'mask_cirrus', 'quicklook_vnir', 'quicklook_swir', 'deadpixelmap']:

            if attrName == 'deadpixelmap':
                # TODO VNIR and SWIR must be merged
                self.logger.warning('Currently, L2A dead pixel masks cannot be saved yet.')
                continue

            outpath = path.join(product_dir, getattr(self.meta, 'filename_%s' % attrName))

            attr_gA = \
                self.generate_quicklook(bands2use=self.meta.preview_bands_vnir) if attrName == 'quicklook_vnir' else \
                self.generate_quicklook(bands2use=self.meta.preview_bands_swir) if attrName == 'quicklook_swir' else \
                getattr(self, attrName)

            if attr_gA is not None:
                attr_gA.save(outpath, **kwargs_save)
                outpaths[attrName] = outpath
            else:
                self.logger.warning("The '%s' attribute cannot be saved because it does not exist in the current EnMAP "
                                    "image." % attrName)

        # TODO remove GDAL's *.aux.xml files?

        # save metadata
        self.meta.add_product_fileinformation(filepaths=list(outpaths.values()))
        metadata_string = self.meta.to_XML()

        with open(outpaths['metaxml'], 'w') as metaF:
            self.logger.info("Writing metadata to %s" % outpaths['metaxml'])
            metaF.write(metadata_string)

        self.logger.info("L2A product successfully written!")

        return product_dir
