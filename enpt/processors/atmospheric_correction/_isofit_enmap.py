# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""EnPT atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""
from types import SimpleNamespace
from tempfile import TemporaryDirectory
from multiprocessing import cpu_count

from isofit.utils.apply_oe import apply_oe

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig

__author__ = 'Daniel Scheffler'





class IsofitEnMAP(object):
    """"""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of IsofitEnMAP."""
        self.cfg = config

    @staticmethod
    def _apply_oe(input_radiance: str,
                  input_loc: str,
                  input_obs: str,
                  working_directory: str,
                  surface_path: str,
                  wavelength_path: str,
                  log_file: str,
                  n_cores: int = cpu_count(),
                  lut_config_file: str = None
                  ):
        apply_oe(
            SimpleNamespace(**dict(
                input_radiance=input_radiance,
                input_loc=input_loc,
                input_obs=input_obs,
                working_directory=working_directory,
                sensor='enmap',
                surface_path=surface_path,
                wavelength_path=wavelength_path,
                n_cores=n_cores,
                presolve=True,
                emulator_base='',
                log_file=log_file,
                pressure_elevation=True,
                lut_config_file=lut_config_file,

                # defaults
                copy_input_files=False,
                modtran_path=None,
                surface_category='multicomponent_surface',
                aerosol_climatology_path=None,
                rdn_factors_path=None,
                atmosphere_type='ATM_MIDLAT_SUMMER',
                channelized_uncertainty_path=None,
                model_discrepancy_path=None,
                multiple_restarts=False,
                logging_level='INFO',
                num_cpus=1,
                memory_gb=-1,
                empirical_line=False,
                analytical_line=False,
                ray_temp_dir="/tmp/ray",
                segmentation_size=40,
                num_neighbors=None,
                atm_sigma=[2],
                prebuilt_lut=''
            ))
        )

    def run(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        with TemporaryDirectory() as td:
            breakpoint()

            self._apply_oe()
