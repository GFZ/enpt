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
from datetime import datetime

import numpy as np
from pyproj.crs import CRS
import isofit
from isofit.core.isofit import Isofit
from isofit.utils.apply_oe import apply_oe
from isofit.utils.template_construction import (
    write_modtran_template,
    get_metadata_from_obs,
    get_metadata_from_loc,
    LUTConfig
)
import ray
from py_tools_ds.geo.coord_grid import get_coord_grid
from py_tools_ds.geo.coord_trafo import transform_coordArray
from geoarray import GeoArray

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
    def _build_modtran_template_file(path_emulator_basedir: str,
                                     path_obs: str,
                                     path_loc: str,
                                     path_workdir: str,
                                     enmap_timestamp: str):
        lut_params = LUTConfig(lut_config_file=None, emulator=path_emulator_basedir)
        (
            h_m_s,
            day_increment,
            mean_path_km,
            mean_to_sensor_azimuth,
            mean_to_sensor_zenith,
            mean_to_sun_zenith,
            mean_relative_azimuth,
            valid,
            to_sensor_zenith_lut_grid,
            to_sun_zenith_lut_grid,
            relative_azimuth_lut_grid,
        ) = get_metadata_from_obs(path_obs, lut_params)

        (
            mean_latitude,
            mean_longitude,
            mean_elevation_km,
            elevation_lut_grid,
        ) = get_metadata_from_loc(
            path_loc, lut_params, pressure_elevation=False
        )

        write_modtran_template(
            atmosphere_type='ATM_MIDLAT_SUMMER',
            fid=enmap_timestamp,
            altitude_km=mean_elevation_km + np.cos(np.deg2rad(180 - mean_to_sensor_zenith)) * mean_path_km,
            dayofyear=datetime.strptime(enmap_timestamp[:15], "%Y%m%dt%H%M%S").timetuple().tm_yday,
            to_sun_zenith=mean_to_sun_zenith,
            to_sensor_azimuth=mean_to_sensor_azimuth,
            to_sensor_zenith=mean_to_sensor_zenith,
            relative_azimuth=mean_relative_azimuth,
            gmtime=float(h_m_s[0] + h_m_s[1] / 60.0),
            elevation_km=mean_elevation_km,
            output_file=pjoin(path_workdir, 'config', f'{enmap_timestamp}_modtran_tpl.json'),
            ihaze_type="AER_NONE",
        )

    def generate_input_files(self, enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        self._generate_radiance_file(enmap_ImageL2, path_outdir)
        self._generate_loc_file(enmap_ImageL2, path_outdir)
        self._generate_obs_file(enmap_ImageL2, path_outdir)

    @staticmethod
    def _generate_radiance_file(enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        # ISOFIT expects radiance in uW/cm²/sr/nm, EnPT provides mW/m²/sr/nm
        # 1000 uW/10000 cm²/sr/nm corresponds to mW/m²/sr/nm
        radiance = enmap_ImageL2.data[:] / 10.0  # TODO consider nodata value

        timestamp = enmap_ImageL2.meta.scene_basename.split('____')[1].split('_')[1]
        gA = GeoArray(radiance, enmap_ImageL2.data.gt, enmap_ImageL2.data.prj)
        gA.meta.band_meta = enmap_ImageL2.data.meta.band_meta
        gA.save(os.path.join(path_outdir, f'{timestamp}_rdn.bsq'))

    @staticmethod
    def _generate_loc_file(enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        xmin, xmax, ymin, ymax = enmap_ImageL2.data.box.boundsMap
        xgsd, ygsd = (enmap_ImageL2.data.xgsd, enmap_ImageL2.data.ygsd)
        x_grid, ygrid = get_coord_grid((xmin, ymax), (xmax, ymin), (xgsd, -ygsd))
        if enmap_ImageL2.data.epsg == 4326:
            lons = x_grid
            lats = ygrid
        else:
            lons, lats = transform_coordArray(
                CRS(enmap_ImageL2.data.epsg).to_wkt(),
                CRS(4326).to_wkt(),
                x_grid, ygrid
            )

        elev = enmap_ImageL2.dem[:]  # FIXME nodata value 0

        timestamp = enmap_ImageL2.meta.scene_basename.split('____')[1].split('_')[1]
        GeoArray(np.dstack([lons, lats, elev]), enmap_ImageL2.data.gt, enmap_ImageL2.data.prj
                 ).save(os.path.join(path_outdir, f'{timestamp}_loc.bsq'))

    @staticmethod
    def _generate_obs_file(enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        path_length = 650000  # from ~650km EnMAP flight height
        vaa = None
        vza = None
        saa = None
        sza = None
        phase = None
        slope = 90
        aspect = 0
        cos_i = None
        utc = None
        earth_sun_dist = enmap_ImageL2.meta.earthSunDist

    def run(self,
            path_toarad: str,
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
                ),
                surface=dict(
                    surface_file=path_surface_file
                ),
            ),
            implementation=dict(
                debug_mode=False,
                n_cores=30,  # FIXME harcoded
                ray_temp_dir='/tmp/ray',  # TODO not Windows-compatible
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

        self._build_modtran_template_file(path_emulator_basedir, path_obs, path_loc, path_workdir, enmap_timestamp)

        Isofit(
            config_file=path_isocfg,
            level='INFO',
            logfile=pjoin(path_outdir, f'{enmap_timestamp}_isofit.log')
        ).run(row_column=None)