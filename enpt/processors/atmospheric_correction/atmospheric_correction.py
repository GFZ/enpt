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

"""EnPT 'atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""
import pprint
import numpy as np
from multiprocessing import cpu_count
from os import path

import sicor
from sicor.sicor_enmap import sicor_ac_enmap
from sicor.options import get_options as get_ac_options

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig
from ...utils.path_generator import get_path_ac_options

__author__ = 'Daniel Scheffler'


class AtmosphericCorrector(object):
    """Class for performing atmospheric correction of EnMAP L1 images using SICOR."""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of AtmosphericCorrector."""
        self.cfg = config

    def get_ac_options(self, enmap_ImageL1: EnMAPL1Product_SensorGeo) -> dict:
        path_opts = get_path_ac_options()

        try:
            options = get_ac_options(path_opts, validation=False)  # FIXME validation is currently not implemented

            # adjust options
            # FIXME this path should be already known to sicor
            options["EnMAP"]["Retrieval"]["fn_LUT"] = \
                path.join(path.abspath(sicor.__path__[0]), 'tables', 'EnMAP_LUT_MOD5_formatted_1nm')
            # options["ECMWF"]["path_db"] = "./ecmwf"  # disbled as it is not needed at the moment

            if enmap_ImageL1.meta.aot is not None:
                options["EnMAP"]["FO_settings"]["aot"] = enmap_ImageL1.meta.aot

            # always use the fast implementation (the slow implementation was only a temporary solution)
            options["EnMAP"]["Retrieval"]["fast"] = True
            options["EnMAP"]["Retrieval"]["ice"] = self.cfg.enable_ice_retrieval
            options["EnMAP"]["Retrieval"]["cpu"] = self.cfg.CPUs or cpu_count()
            options["EnMAP"]["Retrieval"]["disable_progressbars"] = self.cfg.disable_progress_bars

            return options

        except FileNotFoundError:
            raise FileNotFoundError('Could not locate options file for atmospheric correction at %s.' % path_opts)

    def run_ac(self, enmap_ImageL1: EnMAPL1Product_SensorGeo) -> EnMAPL1Product_SensorGeo:
        options = self.get_ac_options(enmap_ImageL1)
        enmap_ImageL1.logger.debug('AC options: \n' + pprint.pformat(options))

        # run AC
        enmap_ImageL1.logger.info("Starting atmospheric correction for VNIR and SWIR detector. "
                                  "Source radiometric unit code is '%s'." % enmap_ImageL1.meta.vnir.unitcode)

        # run SICOR
        # NOTE: - enmap_l2a_vnir, enmap_l2a_swir: reflectance between 0 and 1
        #       - cwv_model, cwc_model, toa_model have the SWIR geometry
        #       - currently, the fast method is implemented,
        #           -> otherwise options["EnMAP"]["Retrieval"]["fast"] must be false
        #       - ice_model is None if self.cfg.enable_ice_retrieval is False
        enmap_l2a_vnir, enmap_l2a_swir, cwv_model, cwc_model, ice_model, toa_model, se, scem, srem = \
            sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=options, logger=enmap_ImageL1.logger)

        # validate results
        for detectordata, detectorname in zip([enmap_l2a_vnir, enmap_l2a_swir], ['VNIR', 'SWIR']):
            mean0, std0 = np.nanmean(detectordata[:, :, 0]), np.nanstd(detectordata[:, :, 0])
            if np.isnan(mean0) or mean0 == 0 or std0 == 0:
                enmap_ImageL1.logger.warning('The atmospheric correction returned empty %s bands!' % detectorname)

        # join results
        enmap_ImageL1.logger.info('Joining results of atmospheric correction.')

        for in_detector, out_detector in zip([enmap_ImageL1.vnir, enmap_ImageL1.swir],
                                             [enmap_l2a_vnir, enmap_l2a_swir]):
            in_detector.data = (out_detector * self.cfg.scale_factor_boa_ref).astype(np.int16)
            # NOTE: geotransform and projection are missing due to sensor geometry

            in_detector.detector_meta.unit = '0-%d' % self.cfg.scale_factor_boa_ref
            in_detector.detector_meta.unitcode = 'BOARef'

            # FIXME what about mask_clouds, mask_clouds_confidence, ac_errors?
            # FIXME use cwv_model, cwc_model, toa_model also for EnPT?

        return enmap_ImageL1
