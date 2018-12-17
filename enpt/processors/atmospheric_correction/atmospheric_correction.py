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

    def get_ac_options(self, buffer_dir):
        path_opts = get_path_ac_options()

        try:
            options = get_ac_options(path_opts)

            # adjust options
            options['EnMAP']['buffer_dir'] = buffer_dir
            for vv in options["RTFO"].values():
                vv["hash_formats"] = dict(spr='%.0f',
                                          coz='%.0f,',
                                          cwv='%.0f,',
                                          tmp='%0f,',
                                          tau_a='%.2f,',
                                          vza='%.0f,')
                vv["disable_progress_bars"] = self.cfg.disable_progress_bars
            options["ECMWF"]["path_db"] = "./ecmwf"

            return options

        except FileNotFoundError:
            raise FileNotFoundError('Could not locate options file for atmospheric correction at %s.' % path_opts)

    def run_ac(self, enmap_ImageL1: EnMAPL1Product_SensorGeo) -> EnMAPL1Product_SensorGeo:
        options = self.get_ac_options(buffer_dir=self.cfg.sicor_cache_dir)
        enmap_ImageL1.logger.debug('AC options: \n' + pprint.pformat(options))

        # run AC
        enmap_ImageL1.logger.info("Starting atmospheric correction for VNIR and SWIR detector. "
                                  "Source radiometric unit code is '%s'." % enmap_ImageL1.meta.vnir.unitcode)
        enmap_l2a_sens_geo, state, cwv_map, ch4_map = sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=options,
                                                                     logger=enmap_ImageL1.logger)

        # join results
        enmap_ImageL1.logger.info('Joining results of atmospheric correction.')

        for in_detector, out_detector in zip([enmap_ImageL1.vnir, enmap_ImageL1.swir],
                                             [enmap_l2a_sens_geo.vnir, enmap_l2a_sens_geo.swir]):
            in_detector.data = (out_detector.data[:] * self.cfg.scale_factor_boa_ref).astype(np.int16)
            # NOTE: geotransform and projection are missing due to sensor geometry

            del in_detector.data_l2a  # FIXME sicor sets data_l2a to float array -> not needed
            del in_detector.unit  # FIXME sicor sets unit to '1' -> not needed

            in_detector.detector_meta.unit = '0-%d' % self.cfg.scale_factor_boa_ref
            in_detector.detector_meta.unitcode = 'BOARef'

            # FIXME what about mask_clouds, mask_clouds_confidence, ac_errors?

        return enmap_ImageL1
