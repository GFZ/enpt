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
import shutil
from types import SimpleNamespace
from tempfile import TemporaryDirectory
from typing import Tuple
from fnmatch import fnmatch
import os
from os.path import join as pjoin
import subprocess
from glob import glob
import json
from collections.abc import Mapping

from isofit.core.isofit import Isofit
from isofit.utils.apply_oe import apply_oe
import isofit
import ray

from ...model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo
from ...options.config import EnPTConfig, path_enptlib

__author__ = 'Daniel Scheffler'


class IsofitEnMAP(object):
    """"""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of IsofitEnMAP."""
        self.cfg = config

        os.environ['SIXS_DIR'] = "/home/gfz-fe/scheffler/6sV2.1"  # FIXME hardcoded
        # os.environ['EMULATOR_PATH'] = '/home/gfz-fe/scheffler/sRTMnet_v100/sRTMnet_v100'  # duplicate of emulator_base

        # make sure ISOFIT's extra-files are downloaded
        # FIXME: somehow ISOFIT expects the data and examples dirs at .../site-packages/, not at ../site-packages/isofit
        isofit_root = isofit.__path__[0]  # .../site-packages/isofit
        if not glob(pjoin(isofit_root, '..', 'data', '*')):
            subprocess.call('isofit download data', shell=True)
        if not glob(pjoin(isofit_root, '..', 'examples', '*')):
            subprocess.call('isofit download examples', shell=True)


    @staticmethod
    def _apply_oe(input_radiance: str,
                  input_loc: str,
                  input_obs: str,
                  working_directory: str,
                  sensor: str = 'enmap',
                  copy_input_files: bool = False,
                  modtran_path: str = None,
                  wavelength_path: str = None,
                  surface_category: str = 'multicomponent_surface',
                  aerosol_climatology_path: str = None,
                  rdn_factors_path: str = None,
                  surface_path: str = None,
                  atmosphere_type: str = 'ATM_MIDLAT_SUMMER',
                  channelized_uncertainty_path: str = None,
                  model_discrepancy_path: str = None,
                  lut_config_file: str = None,
                  multiple_restarts: bool = False,
                  logging_level: str = 'INFO',
                  log_file: str = None,
                  n_cores: int = 1,
                  num_cpus: int = 1,
                  memory_gb: int = -1,
                  presolve: bool = False,
                  empirical_line: bool = False,
                  analytical_line: bool = False,
                  ray_temp_dir: str = "/tmp/ray",
                  emulator_base: str = None,
                  segmentation_size: int = 40,
                  num_neighbors: Tuple[int] = (),
                  atm_sigma: Tuple[float] = (2.0,),
                  pressure_elevation: bool = False,
                  prebuilt_lut: str = None
                  ):
        params = {k: v for k, v in locals().items() if not k.startswith('__')}

        try:
            apply_oe(SimpleNamespace(**params))

        except FileNotFoundError as e:
            print('Attempt to run apply_oe() failed due to FileNotFoundError.')
            if fnmatch(os.path.dirname(e.filename), '*/site-packages/data'):
                raise FileNotFoundError(e.filename,
                                        'The ISOFIT extra-files (data directory) are not properly downloaded.')
            elif fnmatch(os.path.dirname(e.filename), '*/site-packages/examples'):
                raise FileNotFoundError(e.filename,
                                        'The ISOFIT extra-files (examples directory) are not properly downloaded.')
            else:
                raise
        
        finally:
            print('Stopping ray.')
            ray.shutdown()  # FIXME: This should be done by ISOFIT itself (calling ray stop --force is not sufficient)

    def apply_oe_on_sensor_geometry(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        with TemporaryDirectory() as td:
            breakpoint()

            self._apply_oe()

    def apply_oe_on_map_geometry(self, enmap_ImageL2: EnMAPL2Product_MapGeo):
        from geoarray import GeoArray

        with TemporaryDirectory() as td:
            toa_rad = enmap_ImageL2.data
            # TODO


            self._apply_oe()

    @staticmethod
    def run(path_toarad: str,
            path_loc: str,
            path_obs: str,
            path_outdir: str,
            path_workdir: str,
            path_enmap_wavelengths: str,
            path_emulator_basedir: str,
            path_surface_file: str
            ):
        path_isocfg_default = pjoin(path_enptlib, 'options', 'isofit_config_default.json')
        path_isocfg = pjoin(path_workdir, 'config', 'isofit_config.json')
        path_isofit_root = isofit.__path__[0]
        path_data = os.path.abspath(pjoin(path_isofit_root, '..', 'data'))
        path_examples = os.path.abspath(pjoin(path_isofit_root, '..', 'examples'))

        shutil.rmtree(path_workdir)

        for d in [
            path_workdir,
            pjoin(path_workdir, 'config'),
            pjoin(path_workdir, 'lut_full'),
            path_outdir
        ]:
            os.makedirs(d, exist_ok=True)

        enmap_timestamp = os.path.basename(path_toarad).split('____')[1].split('_')[1]

        with open(path_isocfg_default) as json_file:
            isocfg_default = json.load(json_file)

        updatedict = dict(
            ISOFIT_base=os.path.dirname(isofit.__path__[0]),
            forward_model=dict(
                instrument=dict(
                    wavelength_file=path_enmap_wavelengths,
                ),
                radiative_transfer=dict(
                    radiative_transfer_engines=dict(
                        vswir=dict(
                            aerosol_model_file=pjoin(path_data, 'aerosol_model.txt'),
                            aerosol_template_file=pjoin(path_data, 'aerosol_template.json'),
                            earth_sun_distance_file=pjoin(path_data, 'earth_sun_distance.txt'),
                            emulator_aux_file=pjoin(os.path.dirname(path_emulator_basedir), 'sRTMnet_v100_aux.npz'),
                            emulator_file=path_emulator_basedir,
                            engine_base_dir=os.environ['SIXS_DIR'],
                            interpolator_base_path=pjoin(path_workdir, 'lut_full', 'sRTMnet_v100_vi'),
                            irradiance_file=pjoin(path_examples, '20151026_SantaMonica/data/prism_optimized_irr.dat'),
                            lut_path=pjoin(path_workdir, 'lut_full', 'lut.nc'),
                            sim_path=pjoin(path_workdir, 'lut_full'),
                            template_file=pjoin(path_workdir, 'config', f'{enmap_timestamp}_modtran_tpl.json')
                        )
                    ),
                    # statevector=dict(
                    #     AOT550=dict(
                    #         init='TBD',  # TODO
                    #         prior_mean='TBD',  # TODO
                    #     ),
                    #     H2OSTR=dict(
                    #         init='TBD',  # TODO
                    #         prior_mean='TBD',  # TODO
                    #     )
                    # ),
                    surface=path_surface_file
                )
            ),
            implementation=dict(
                debug_mode=False,
                n_cores=30,
                ray_temp_dir='/tmp/ray',
            ),
            input=dict(
                measured_radiance_file=path_toarad,
                loc_file=path_loc,
                obs_file=path_obs
            ),
            output=dict(
                estimated_reflectance_file=pjoin(path_outdir, f'{enmap_timestamp}_estimated_reflectance.bsq'),
                estimated_state_file=pjoin(path_outdir, f'{enmap_timestamp}_estimated_state.bsq'),
                posterior_uncertainty_file=pjoin(path_outdir, f'{enmap_timestamp}_posterior_uncertainty.bsq'),
            )
        )

        def update_nested_dict(d, u):
            for k, v in u.items():
                if isinstance(v, Mapping):
                    d[k] = update_nested_dict(d.get(k, {}), v)
                else:
                    d[k] = v
            return d

        isocfg = update_nested_dict(isocfg_default, updatedict)

        with open(path_isocfg, 'w') as json_file:
            json.dump(isocfg, json_file, skipkeys=False, indent=4)

        Isofit(
            config_file=path_isocfg,
            level='INFO',
            logfile=pjoin(path_outdir, f'{enmap_timestamp}_isofit.log')
        ).run(row_column=None)