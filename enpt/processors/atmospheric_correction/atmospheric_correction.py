# -*- coding: utf-8 -*-
"""EnPT 'atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""
import pprint
import numpy as np

from sicor.sicor_enmap import sicor_ac_enmap
from sicor.options import get_options as get_ac_options

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig
from ...utils.path_generator import get_path_ac_options


class AtmosphericCorrector(object):
    """Class for performing atmospheric correction of EnMAP L1 images using SICOR."""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of AtmosphericCorrector."""
        self.cfg = config

    def get_ac_options(self):
        path_opts = get_path_ac_options()

        try:
            options = get_ac_options(path_opts, validation=False)  # FIXME validation is currently not implemented

            # adjust options
            options["EnMAP"]["Retrieval"]["fn_LUT"] = \
                "/home/gfz-fe/scheffler/python_deployed/sicor_aux_data/EnMAP_LUT_MOD5_formatted_1nm"  # FIXME hardcoded
            # options["ECMWF"]["path_db"] = "./ecmwf"  # disbled as it is not needed at the moment
            # TODO disable_progress_bars?

            options["EnMAP"]["Retrieval"]["ice"] = self.cfg.enable_ice_retrieval

            return options

        except FileNotFoundError:
            raise FileNotFoundError('Could not locate options file for atmospheric correction at %s.' % path_opts)

    def run_ac(self, enmap_ImageL1: EnMAPL1Product_SensorGeo) -> EnMAPL1Product_SensorGeo:
        options = self.get_ac_options()
        enmap_ImageL1.logger.debug('AC options: \n' + pprint.pformat(options))

        # run AC
        enmap_ImageL1.logger.info("Starting atmospheric correction for VNIR and SWIR detector. "
                                  "Source radiometric unit code is '%s'." % enmap_ImageL1.meta.vnir.unitcode)

        # run SICOR
        # NOTE: - enmap_l2a_vnir, enmap_l2a_swir: reflectance between 0 and 1
        #       - cwv_model, cwc_model, toa_model have the SWIR geometry
        #       - currently, the slower method is implemented,
        #           -> otherwise options["EnMAP"]["Retrieval"]["fast"] must be true  # TODO
        if self.cfg.enable_ice_retrieval:
            enmap_l2a_vnir, enmap_l2a_swir, cwv_model, cwc_model, ice_model, toa_model = \
                sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=options, logger=enmap_ImageL1.logger)
        else:
            enmap_l2a_vnir, enmap_l2a_swir, cwv_model, cwc_model, toa_model = \
                sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=options, logger=enmap_ImageL1.logger)

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
