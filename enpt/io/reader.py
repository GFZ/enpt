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

"""Reader module for reading all kinds of EnMAP images."""

import logging
import os
from fnmatch import filter
import numpy as np
from scipy.interpolate import interp1d

from ..model.images import EnMAPL1Product_SensorGeo
from ..options.config import EnPTConfig

__author__ = 'Daniel Scheffler'


class L1B_Reader(object):
    """Reader for EnMAP Level-1B products."""

    def __init__(self,
                 config: EnPTConfig,
                 logger: logging.Logger = None,
                 root_dir_main: str = None,
                 root_dir_ext: str = None,
                 n_line_ext: int = None):
        """Get an instance of L1B_Reader.

        :param config:  instance of EnPTConfig class
        :param logger:  instance of logging.Logger (NOTE: This logger is only used to log messages within L1B_Reader.
                                                          It is not appended to the read L1B EnMAP object).
        :param root_dir_main: Root directory of EnMAP Level-1B product (the main image)
        :param root_dir_ext:  Root directory of EnMAP Level-1B product (to extend the main image)
        :param n_line_ext:    [Optional] add number of line to be added to the main image from the extended image
        """
        self.cfg = config
        self.logger = logger or logging.getLogger(__name__)

        # read data if root_dir_main is given or not
        if root_dir_main is not None:
            self.read_inputdata(root_dir_main, root_dir_ext, n_line_ext)

    def read_inputdata(self,
                       root_dir_main,
                       root_dir_ext: str = None,
                       n_line_ext: int = None,
                       compute_snr: bool = True) -> EnMAPL1Product_SensorGeo:
        """Read L1B EnMAP data. Extend the image by adding a second image (entire, partial).

        :param root_dir_main: Root directory of the main EnMAP Level-1B product
        :param root_dir_ext:  Root directory of the extended EnMAP Level-1B product (optional)
        :param n_line_ext:    Number of lines to be added to the main image (if None, use the whole image)
        :param compute_snr:   whether to compute SNR or not (default: True)
        :return: instance of EnMAPL1Product_SensorGeo
        """
        self.validate_input(root_dir_main, root_dir_ext)
        self.logger.info("Reading Input Data")

        # Get a new instance of the EnMAPL1Product_SensorGeo for the main image (TOA radiance)
        l1b_main_obj = EnMAPL1Product_SensorGeo(root_dir_main, config=self.cfg, logger=self.logger)

        # append a second EnMAP L1B product below if given
        if root_dir_ext:
            l1b_main_obj.append_new_image(root_dir_ext, n_line_ext)

        # compute SNR
        if compute_snr:
            l1b_main_obj.calc_snr_from_radiance()

        # compute geolayer if not already done
        if l1b_main_obj.meta.vnir.lons is None or l1b_main_obj.meta.vnir.lats is None:
            l1b_main_obj.meta.vnir.lons, l1b_main_obj.meta.vnir.lats = \
                l1b_main_obj.meta.vnir.compute_geolayer_for_cube()

        if l1b_main_obj.meta.swir.lons is None or l1b_main_obj.meta.swir.lats is None:
            l1b_main_obj.meta.swir.lons, l1b_main_obj.meta.swir.lats = \
                l1b_main_obj.meta.swir.compute_geolayer_for_cube()

        # l1b_main_obj.correct_VNIR_SWIR_shift()

        # Validate and return the l1b_main_obj
        self.validate_output()
        return l1b_main_obj

    def validate_input(self, root_dir_main, root_dir_ext):
        """Validate user inputs."""
        if not self.cfg.is_dummy_dataformat:
            self._validate_enmap_l1b_rootdir(root_dir_main)
            if root_dir_ext:
                self._validate_enmap_l1b_rootdir(root_dir_ext)

    @staticmethod
    def _validate_enmap_l1b_rootdir(rootdir_l1b):
        """Check for valid EnMAP L1B root directory."""
        if not os.path.isdir(rootdir_l1b):
            raise NotADirectoryError(rootdir_l1b, 'EnMAP images have to be provided either as zip-archives or as '
                                                  'a directory containing all extracted files.')

        files = os.listdir(rootdir_l1b)

        if not files:
            raise RuntimeError("The root directory of the EnMAP image %s is empty." % rootdir_l1b)

        for pattern in [
            # '*-HISTORY.XML',  # only included in internal DLR test data, not in the zip archive on enmap.org
            # '*-LOG.XML',  # only included in internal DLR test data, not in the zip archive on enmap.org
            '*-METADATA.XML',
            '*-QL_PIXELMASK_SWIR.TIF',
            '*-QL_PIXELMASK_VNIR.TIF',
            '*-QL_QUALITY_CIRRUS.TIF',
            '*-QL_QUALITY_CLASSES.TIF',
            '*-QL_QUALITY_CLOUD.TIF',
            '*-QL_QUALITY_CLOUDSHADOW.TIF',
            '*-QL_QUALITY_HAZE.TIF',
            '*-QL_QUALITY_SNOW.TIF',
            '*-QL_QUALITY_TESTFLAGS_SWIR.TIF',
            '*-QL_QUALITY_TESTFLAGS_VNIR.TIF',
            '*-QL_SWIR.TIF',
            '*-QL_VNIR.TIF',
            '*-SPECTRAL_IMAGE_SWIR.TIF',
            '*-SPECTRAL_IMAGE_VNIR.TIF',
        ]:
            if not filter(files, pattern) and not filter(files, pattern.replace('.TIF', '.GEOTIFF')):
                raise FileNotFoundError('The root directory of the EnMAP image %s misses a file with the pattern %s.'
                                        % (rootdir_l1b, pattern))

    def validate_output(self):
        """Validate outputs of L1B_Reader."""
        pass


def Solar_Irradiance_reader(path_solar_irr_model: str, resol_nm: float = None, wvl_min_nm: float = None,
                            wvl_max_nm: float = None) -> np.ndarray:
    """Read the given solar irradiance file and return an array of irradiances.

    :param path_solar_irr_model:    file path to solar irradiance model

                                    -> must be arranged like that:
                                       col0 = Wavelength[nm]; col1 = Solar Irradiance [W/m2/µm])
    :param resol_nm:                spectral resolution for returned irradiances [nanometers]
    :param wvl_min_nm:              minimum wavelength of returned irradiances [nanometers]
    :param wvl_max_nm:              maximum wavelength of returned irradiances [nanometers]
    :return:
    """
    sol_irr = np.loadtxt(path_solar_irr_model, skiprows=1)
    if resol_nm is not None and isinstance(resol_nm, (float, int)):
        wvl_min = (np.min(sol_irr[:, 0]) if wvl_min_nm is None else wvl_min_nm)
        wvl_max = (np.max(sol_irr[:, 0]) if wvl_max_nm is None else wvl_max_nm)
        wvl_rsp = np.arange(wvl_min, wvl_max, resol_nm)
        sol_irr = interp1d(sol_irr[:, 0], sol_irr[:, 1], kind='linear')(wvl_rsp)
    return sol_irr
