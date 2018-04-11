# -*- coding: utf-8 -*-
"""EnPT 'atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""
import tempfile
import pprint

from sicor.sicor_enmap import sicor_ac_enmap
from sicor.options import get_options as get_ac_options

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig
from ...utils.path_generator import get_path_ac_options, get_path_ac_aerosol_table, get_path_ac_ch4_table


class AtmosphericCorrector(object):
    """Class for performing atmospheric correction of EnMAP L1 images using SICOR."""

    def __init__(self, config: EnPTConfig=None):
        """Create an instance of AtmosphericCorrector."""
        self.cfg = config

    @staticmethod
    def get_ac_options(buffer_dir):
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
            options["ECMWF"]["path_db"] = "./ecmwf"
            path_tbl_aerosol = get_path_ac_aerosol_table()
            path_tbl_ch4 = get_path_ac_ch4_table()
            for name, val in options["RTFO"].items():
                if "aerosol" in name:
                    val['atm_tables_fn'] = path_tbl_aerosol
                if "ch4" in name:
                    val['atm_tables_fn'] = path_tbl_ch4

            return options

        except FileNotFoundError:
            raise FileNotFoundError('Could not locate options file for atmospheric correction at %s.' % path_opts)

    def run_ac(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        with tempfile.TemporaryDirectory(self.cfg.working_dir) as buffer_dir:
            options = self.get_ac_options(buffer_dir)
            enmap_ImageL1.logger.debug('AC options: \n' + pprint.pformat(options))

            # run AC
            enmap_ImageL1.logger.info("Starting atmospheric correction for VNIR and SWIR detector. "
                                      "Source radiometric unit code is '%s'." % enmap_ImageL1.meta.vnir.unitcode)
            enmap_l2a_sens_geo, state, fits = sicor_ac_enmap(enmap_l1b=enmap_ImageL1, options=options,
                                                             logger=enmap_ImageL1.logger, debug=True)

        return enmap_l2a_sens_geo
