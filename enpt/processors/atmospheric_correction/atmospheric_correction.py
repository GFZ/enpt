# -*- coding: utf-8 -*-
"""EnPT 'atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""

from sicor.sicor_enmap import sicor_ac_enmap
from sicor.options import get_options as get_ac_options

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig
from ...utils.path_generator import get_path_ac_options


class AtmosphericCorrector(object):
    """Class for performing atmospheric correction of EnMAP L1 images using SICOR."""

    def __init__(self, config: EnPTConfig=None):
        """Create an instance of AtmosphericCorrector."""
        self.cfg = config

    @property
    def options(self):
        path_opts = get_path_ac_options()

        try:
            return get_ac_options(path_opts)
        except FileNotFoundError:
            raise FileNotFoundError('Could not locate options file for atmospheric correction at %s.' % path_opts)

    def run_ac(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        enmap_ImageL1.logger.info("Starting atmospheric correction for VNIR and SWIR detector. "
                                  "Source radiometric unit code is '%s'." % enmap_ImageL1.meta.vnir.unitcode)

        # run AC
        enmap_l2a_sens_geo, state, fits = sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=self.options,
                                                         logger=enmap_ImageL1.logger, debug=True)

        return enmap_l2a_sens_geo
