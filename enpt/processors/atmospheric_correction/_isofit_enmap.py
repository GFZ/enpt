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
from tempfile import TemporaryDirectory
from typing import Tuple, List
from fnmatch import fnmatch
import os
from os.path import join as pjoin
from pathlib import Path
from glob import glob
import json
from collections.abc import Mapping
from datetime import datetime
from zipfile import ZipFile
from multiprocessing import cpu_count

import numpy as np
from pyproj.crs import CRS
from pandas import DataFrame
import netCDF4 as nc
from scipy.interpolate import interp1d

from ...utils import EnvContextManager
with EnvContextManager(ISOFIT_DEBUG='0',
                       MKL_NUM_THREADS='1',
                       OMP_NUM_THREADS='1'):
    import isofit
    from isofit.core.isofit import Isofit
    from isofit.utils import surface_model, analytical_line, extractions, segment
    from isofit.data.cli.data import download as download_data
    from isofit.data.cli.examples import download as download_examples
    from isofit.utils.apply_oe import apply_oe, CHUNKSIZE
    from isofit.utils.template_construction import (
        write_modtran_template,
        get_metadata_from_obs,
        get_metadata_from_loc,
        LUTConfig,
        Pathnames
    )
from py_tools_ds.geo.coord_grid import get_coord_grid
from py_tools_ds.geo.coord_trafo import transform_coordArray
from geoarray import GeoArray

from ...model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo
from ...options.config import EnPTConfig, path_enptlib

__author__ = 'Daniel Scheffler'


class IsofitEnMAP(object):
    """"""

    def __init__(self,
                 config: EnPTConfig = None,
                 log_level: str = None
                 ) -> None:
        """Create an instance of IsofitEnMAP."""
        self.cfg = config
        self.log_level = log_level or config.log_level if config else 'INFO'

        os.environ['SIXS_DIR'] = pjoin(Path.home(), '.isofit', 'sixs')
        # os.environ['EMULATOR_PATH'] = '/home/gfz-fe/scheffler/srtmnet/sRTMnet_v120.h5'  # duplicate of emulator_base

        # make sure ISOFIT's extra-files are downloaded
        download_data(path=None, tag="latest")
        download_examples(path=None, tag="latest")

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
            mean_to_sun_azimuth,
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
            altitude_km=mean_elevation_km + np.cos(np.deg2rad(mean_to_sensor_zenith)) * mean_path_km,
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
        fp_rad = self._generate_radiance_file(enmap_ImageL2, path_outdir)
        fp_loc = self._generate_loc_file(enmap_ImageL2, path_outdir)
        fp_obs = self._generate_obs_file(enmap_ImageL2, path_outdir)
        fp_wvl = self._generate_wavelength_file(enmap_ImageL2, path_outdir)
        fp_surf = self._generate_surface_file(fp_wvl, path_outdir)
        fp_lut = self._generate_lut_file(path_outdir,
                                         enmap_ImageL2.meta.geom_sun_zenith)  # TODO: set LUT to None in case of 6S

        return fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut

    @staticmethod
    def _generate_radiance_file(enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        # ISOFIT expects radiance in uW/cm²/sr/nm, EnPT provides mW/m²/sr/nm
        # 1000 uW/10000 cm²/sr/nm corresponds to mW/m²/sr/nm
        mask_nodata = ~enmap_ImageL2.data.mask_nodata[:]
        radiance = enmap_ImageL2.data[:] / 10.0
        radiance[mask_nodata] = -9999

        fp_out = pjoin(path_outdir, f"{enmap_ImageL2.meta.scene_basename}_rdn")
        gA = GeoArray(radiance, enmap_ImageL2.data.gt, enmap_ImageL2.data.prj, nodata=-9999)
        gA.meta.band_meta = enmap_ImageL2.data.meta.band_meta
        gA.save(fp_out)

        return fp_out

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
        # FIXME: shape of elev does not match shape of lons/lats if no DEM is provided
        # elev = np.zeros_like(lons)  # use mean elevation

        loc_data = np.dstack([lons, lats, elev])
        loc_data[~enmap_ImageL2.data.mask_nodata[:]] = -9999

        fp_out = pjoin(path_outdir, f"{enmap_ImageL2.meta.scene_basename}_loc")  # no file extension supported
        GeoArray(loc_data,
                 enmap_ImageL2.data.gt, enmap_ImageL2.data.prj,
                 bandnames=[
                     'Longitude (WGS-84)',
                     'Latitude (WGS-84)',
                     'Elevation (m)'  # used as first guess in case this is defined as state_vector component in the config
                 ],
                 nodata = -9999
                 ).save(fp_out)

        return fp_out

    def _generate_obs_file(self, enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        path_length = np.full(enmap_ImageL2.data.shape[:2], fill_value=650000)  # from ~650km EnMAP flight height
        vaa = enmap_ImageL2.meta.geom_view_azimuth_array
        vza = enmap_ImageL2.meta.geom_view_zenith_array
        saa = enmap_ImageL2.meta.geom_sun_azimuth_array
        sza = enmap_ImageL2.meta.geom_sun_zenith_array
        phase = self._compute_solar_phase(vaa, vza, saa, sza)
        slope = np.full(enmap_ImageL2.data.shape[:2], fill_value=90)
        aspect = np.zeros(enmap_ImageL2.data.shape[:2])
        cos_i = self._compute_cos_i(saa, sza, slope=90, aspect=0)
        utc = enmap_ImageL2.meta.aqtime_utc_array  # TODO pixel-wise values
        earth_sun_dist = np.full(enmap_ImageL2.data.shape[:2], fill_value=enmap_ImageL2.meta.earthSunDist)  # TODO pixel-wise values

        obs_data = np.dstack([path_length, vaa, vza, saa, sza, phase, slope, aspect, cos_i, utc, earth_sun_dist])
        obs_data[~enmap_ImageL2.data.mask_nodata[:]] = -9999

        fp_out = pjoin(path_outdir, f"{enmap_ImageL2.meta.scene_basename}_obs")  # no file extension supported
        GeoArray(obs_data,
                 enmap_ImageL2.data.gt, enmap_ImageL2.data.prj,
                 bandnames=[
                     'Path length (m)',
                     'To-sensor azimuth (0 to 360 degrees cw from N)',
                     'To-sensor zenith (0 to 90 degrees from zenith)',
                     'To-sun azimuth (0 to 360 degrees cw from N)',
                     'To-sun zenith (0 to 90 degrees from zenith)',
                     'Solar phase',
                     'Slope',
                     'Aspect',
                     'Cosine(i)',
                     'UTC Time',
                     'Earth-sun distance (AU)'
                 ],
                 nodata=-9999
                 ).save(fp_out)

        return fp_out

    @staticmethod
    def _generate_wavelength_file(enmap_ImageL2: EnMAPL2Product_MapGeo, path_outdir: str):
        fp_out = pjoin(path_outdir, 'enmap_wavelength_fwhm.txt')
        wvl = enmap_ImageL2.meta.wvl_center
        fwhm = enmap_ImageL2.meta.fwhm
        df = DataFrame(np.hstack([wvl.reshape(-1, 1), fwhm.reshape(-1, 1)]))
        df.to_csv(fp_out, header=False, sep='\t')

        return fp_out

    @staticmethod
    def _generate_surface_file(path_wavelength_file: str, path_outdir: str):
        fp_out = pjoin(path_outdir, 'surface_enmap.mat')
        # fp_surfjson = pjoin(path_enptlib, 'options', 'isofit_surface_default.json')
        fp_surfjson = pjoin(path_enptlib, 'options', 'isofit_surface_20240103_REE.json')

        with (ZipFile(pjoin(path_enptlib, 'resources', 'isofit', 'isofit_surface_spectra.zip'), "r") as zf,
              TemporaryDirectory() as td):
            zf.extractall(td)
            fp_surfjson_tmp = pjoin(td, os.path.basename(fp_surfjson))
            shutil.copyfile(fp_surfjson, fp_surfjson_tmp)

            surface_model(
                config_path=fp_surfjson_tmp,
                wavelength_path=path_wavelength_file,
                output_path=fp_out
            )
            assert os.path.isfile(fp_out)

        return fp_out

    @staticmethod
    def _generate_lut_file(path_outdir: str, sza_scene: float):
        # TODO: By re-using either the unpacked lut.zip or the LUT_ISOFIT.nc,
        #       the processing time can be reduced by ~20-60 sec.
        fp_out = pjoin(path_outdir, 'EnMAP_LUT_MOD5_ISOFIT_formatted_1nm.nc')

        with (ZipFile(pjoin(path_enptlib, 'resources', 'isofit', 'lut.zip'), 'r') as zf,
              TemporaryDirectory() as td):
            zf.extractall(td)

            LUT_Transformer(
                path_lut=os.path.join(td, 'EnMAP_LUT_MOD5_formatted_1nm'),
                sza_scene=sza_scene
            ).read_binary_modtran_lut(
                path_out_nc=fp_out
            )

            assert os.path.isfile(fp_out)

        return fp_out

    @staticmethod
    def _compute_solar_phase(vaa, vza, saa, sza):
        vaa_r, vza_r, saa_r, sza_r = (np.deg2rad(i) for i in (vaa, vza, saa, sza))
        vx = np.sin(vza_r) * np.cos(vaa_r)
        vy = np.sin(vza_r) * np.sin(vaa_r)
        vz = np.cos(vza_r)
        sx = np.sin(sza_r) * np.cos(saa_r)
        sy = np.sin(sza_r) * np.sin(saa_r)
        sz = np.cos(sza_r)
        phase_r = np.arccos(
            (vx * sx + vy * sy + vz * sz) / np.sqrt((vx ** 2 + vy ** 2 + vz ** 2) * (sx ** 2 + sy ** 2 + sz ** 2)))
        phase = np.rad2deg(phase_r)

        return phase

    @staticmethod
    def _compute_cos_i(saa, sza, slope: float, aspect: float):
        saa_r, sza_r, slope_r = (np.deg2rad(i) for i in (saa, sza, slope))
        saa_asp_r = saa_r if aspect == 0 else np.deg2rad(saa - aspect)
        cos_i = np.cos(sza_r) * np.cos(slope_r) + np.sin(sza_r) * np.sin(slope_r) * np.cos(saa_asp_r)

        return cos_i

    def _apply_oe(self,
                  input_radiance: str,
                  input_loc: str,
                  input_obs: str,
                  working_directory: str,
                  surface_path: str,
                  sensor: str = 'enmap',
                  copy_input_files: bool = False,
                  modtran_path: str = None,
                  wavelength_path: str = None,
                  surface_category: str = 'multicomponent_surface',
                  aerosol_climatology_path: str = None,
                  rdn_factors_path: str = None,
                  atmosphere_type: str = 'ATM_MIDLAT_SUMMER',
                  channelized_uncertainty_path: str = None,
                  model_discrepancy_path: str = None,
                  lut_config_file: str = None,
                  multiple_restarts: bool = False,
                  logging_level: str = None,
                  log_file: str = None,
                  n_cores: int = 1,
                  presolve: bool = False,
                  empirical_line: bool = False,
                  analytical_line: bool = False,
                  ray_temp_dir: str = "/tmp/ray",
                  emulator_base: str = None,
                  segmentation_size: int = 40,
                  num_neighbors: Tuple[int] = (),
                  atm_sigma: Tuple[float] = (2.0,),
                  pressure_elevation: bool = False,
                  prebuilt_lut: str = None,
                  no_min_lut_spacing: bool = False,
                  inversion_windows: List[float] = None,
                  config_only: bool = False,
                  interpolate_bad_rdn=False,
                  interpolate_inplace=False,
                  ):
        logging_level = logging_level or self.log_level
        params = {k: v for k, v in locals().items() if not k.startswith('__')}

        try:
            apply_oe(**params)

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
            import ray
            ray.shutdown()  # FIXME: This should be done by ISOFIT itself (calling ray stop --force is not sufficient)

    def apply_oe_on_sensor_geometry(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        with TemporaryDirectory() as td:
            self._apply_oe()

    def apply_oe_on_map_geometry(self, enmap_ImageL2: EnMAPL2Product_MapGeo):
        with TemporaryDirectory() as td:
            path_indir = pjoin(td, 'input')
            fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut = self.generate_input_files(enmap_ImageL2, path_indir)

            os.makedirs(pjoin(td, 'workdir'), exist_ok=True)
            os.makedirs(pjoin(td, 'input'), exist_ok=True)
            os.makedirs(pjoin(td, 'output'), exist_ok=True)

            self._apply_oe(
                input_radiance=fp_rad,
                input_loc=fp_loc,
                input_obs=fp_obs,
                working_directory=pjoin(td, 'workdir'),
                surface_path=fp_surf,
                wavelength_path=fp_wvl,
                log_file=pjoin(td, 'output', 'isofit.log'),
                presolve=True,
                emulator_base=pjoin(Path.home(), '.isofit', 'srtmnet', 'sRTMnet_v120.h5'),
                n_cores=30,  # FIXME hardcoded
                prebuilt_lut=fp_lut
            )

            # read the AC results back into memory
            boa_rfl = GeoArray(glob(pjoin(td, 'output', 'estimated_reflectance.bsq'))[0])
            # state = GeoArray(glob(pjoin(td, 'output', 'estimated_state.bsq'))[0])
            # uncert = GeoArray(glob(pjoin(td, 'output', 'posterior_uncertainty.bsq'))[0])
            boa_rfl.to_mem()

            return boa_rfl

    def _run(self,
             path_toarad: str,
             path_loc: str,
             path_obs: str,
             path_outdir: str,
             path_workdir: str,
             path_enmap_wavelengths: str,
             path_emulator_basedir: str,
             path_surface_file: str,
             path_lut: str = None,
             aot: float = None,
             cwv: float = None,
             segmentation: bool = False,
             segmentation_size: int = 40,
             n_cores: int = cpu_count()
             ) -> dict:
        enmap_timestamp = os.path.basename(path_toarad).split('____')[1].split('_')[1]
        path_isocfg_default = pjoin(path_enptlib, 'options', 'isofit_config_default.json')
        path_isocfg = pjoin(path_workdir, 'config', 'isofit_config.json')
        path_data = os.path.abspath(pjoin(Path.home(), '.isofit', 'data'))
        path_examples = os.path.abspath(pjoin(Path.home(), '.isofit', 'examples'))
        path_logfile = pjoin(path_outdir, f'{enmap_timestamp}_isofit.log')

        if os.path.isdir(path_workdir):
            if path_lut and path_lut.startswith(path_workdir):
                raise ValueError(path_lut, "The given prebuilt LUT must not be within the given "
                                           "working directory as it is deleted before running ISOFIT.")
            shutil.rmtree(path_workdir)

        for d in [
            path_workdir,
            pjoin(path_workdir, 'config'),
            pjoin(path_workdir, 'lut_full'),
            path_outdir,
            os.path.dirname(path_toarad)  # input directory
        ]:
            os.makedirs(d, exist_ok=True)

        with open(path_isocfg_default) as json_file:
            isocfg_default = json.load(json_file)

        updatedict = dict(
            ISOFIT_base=os.path.dirname(isofit.__path__[0]),
            forward_model=dict(
                instrument=dict(
                    # SNR=500,  # use noise file instead of static SNR
                    parametric_noise_file=pjoin(path_enptlib, 'resources', 'EnMAP_Sensor', 'parametric_noise.txt'),
                    wavelength_file=path_enmap_wavelengths,
                ),
                radiative_transfer=dict(
                    radiative_transfer_engines=dict(
                        vswir=dict(
                            aerosol_model_file=pjoin(path_data, 'aerosol_model.txt'),
                            aerosol_template_file=pjoin(path_data, 'aerosol_template.json'),
                            earth_sun_distance_file=pjoin(path_data, 'earth_sun_distance.txt'),
                            emulator_aux_file=pjoin(path_emulator_basedir, 'sRTMnet_v120_aux.npz'),
                            emulator_file=pjoin(path_emulator_basedir, 'sRTMnet_v120.h5'),
                            engine_base_dir=pjoin(Path.home(), '.isofit', 'sixs'),
                            interpolator_base_path=pjoin(path_workdir, 'lut_full', 'sRTMnet_v120_vi'),
                            irradiance_file=pjoin(path_examples, '20151026_SantaMonica/data/prism_optimized_irr.dat'),
                            lut_path=path_lut or pjoin(path_workdir, 'lut_full', 'lut.nc'),  # if existing -> use this one, otherwise simulate to lut.nc
                            sim_path=pjoin(path_workdir, 'lut_full'),
                            template_file=pjoin(path_workdir, 'config', f'{enmap_timestamp}_modtran_tpl.json')
                        )
                    ),
                    statevector=dict(
                        AERFRAC_2=dict(
                            init=aot,
                            prior_mean=aot,
                        ) if aot else isocfg_default['forward_model']['radiative_transfer']['statevector']['AERFRAC_2'],
                        H2OSTR=dict(
                            init=cwv,
                            prior_mean=cwv
                        ) if cwv else isocfg_default['forward_model']['radiative_transfer']['statevector']['H2OSTR'],
                    ),
                ),
                surface=dict(
                    surface_file=path_surface_file
                ),
            ),
            implementation=dict(
                # debug_mode=False,
                debug_mode=True,  # TODO deactivate if done
                n_cores=n_cores,
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
        paths = Pathnames(
            input_radiance=path_toarad,
            input_loc=path_loc,
            input_obs=path_obs,
            working_directory=os.path.abspath(pjoin(path_workdir, '..')),
            surface_path=path_surface_file,
            aerosol_climatology_path=None,
            sensor='enmap',
            copy_input_files=False,
            channelized_uncertainty_path=None,
            model_discrepancy_path=None,
            modtran_path=None,
            rdn_factors_path=None,
            ray_temp_dir='/tmp/ray',  # FIXME not Windows-compatible
            interpolate_inplace=False
        )

        with EnvContextManager(
                MKL_NUM_THREADS='1',
                OMP_NUM_THREADS='1',
                ISOFIT_DEBUG='0'
        ):
            try:
                # Superpixel segmentation
                if segmentation:
                    if not os.path.exists(paths.lbl_working_path) or \
                       not os.path.exists(paths.radiance_working_path
                    ):
                        # logging.info("Segmenting...")
                        segment(
                            spectra=(paths.radiance_working_path, paths.lbl_working_path),
                            nodata_value=-9999,  # as set in self._generate_radiance_file()
                            npca=5,
                            segsize=segmentation_size,
                            nchunk=CHUNKSIZE,
                            n_cores=n_cores,
                            loglevel=self.log_level,
                            logfile=path_logfile,
                        )

                    # Extract input data per segment
                    for inp, outp in [
                        (paths.radiance_working_path, paths.rdn_subs_path),
                        (paths.obs_working_path, paths.obs_subs_path),
                        (paths.loc_working_path, paths.loc_subs_path),
                    ]:
                        if not os.path.exists(outp):
                            if not os.path.exists(os.path.dirname(outp)):
                                os.makedirs(os.path.dirname(outp))

                            # logging.info("Extracting " + outp)
                            extractions(
                                inputfile=inp,
                                labels=paths.lbl_working_path,
                                output=outp,
                                chunksize=CHUNKSIZE,
                                flag=-9999,
                                n_cores=n_cores,
                                loglevel=self.log_level,
                                logfile=path_logfile,
                            )

                    # enable segmentation and update input/output files accordingly
                    update_nested_dict(
                        isocfg,
                        dict(
                            forward_model=dict(
                                instrument=dict(
                                    integrations=segmentation_size
                                )
                            ),
                            input=dict(
                                measured_radiance_file=paths.rdn_subs_path,
                                loc_file=paths.loc_subs_path,
                                obs_file=paths.obs_subs_path,
                            ),
                            output=dict(
                                estimated_reflectance_file=paths.rfl_subs_path,
                                estimated_state_file=paths.state_subs_path,
                                posterior_uncertainty_file=paths.uncert_subs_path,
                                atmospheric_coefficients_file=paths.atm_coeff_path
                            )
                        )
                    )

                with open(path_isocfg, 'w') as json_file:
                    json.dump(isocfg, json_file, skipkeys=False, indent=4)

                self._build_modtran_template_file(
                    path_emulator_basedir,
                    path_obs,
                    path_loc,
                    path_workdir,
                    enmap_timestamp
                )

                Isofit(
                    config_file=path_isocfg,
                    level=self.log_level,
                    logfile=path_logfile
                ).run(row_column=None)

                if segmentation:
                    # logging.info("Analytical line inference")
                    analytical_line(
                        rdn_file=paths.radiance_working_path,
                        loc_file=paths.loc_working_path,
                        obs_file=paths.obs_working_path,
                        isofit_dir=paths.working_directory,
                        isofit_config=path_isocfg,  # FIXME equalize with paths.isofit_full_config_path
                        segmentation_file=paths.lbl_working_path,
                        n_atm_neighbors=[int(round(3950 / 9 - 35 / 36 * segmentation_size))],
                        n_cores=n_cores,
                        smoothing_sigma=[2],
                        output_rfl_file=paths.rfl_working_path,
                        output_unc_file=paths.uncert_working_path,
                        # atm_file=None,
                        loglevel=self.log_level,
                        logfile=path_logfile
                    )

                    return dict(
                        estimated_reflectance_file=paths.rfl_working_path,
                        estimated_state_file=paths.state_working_path,
                        posterior_uncertainty_file=paths.uncert_working_path,
                    )
                else:
                    return isocfg['output']

            finally:
                # if os.environ.get('ISOFIT_DEBUG') != '1':
                print('Stopping ray.')
                import ray
                ray.shutdown()  # FIXME: This should be done by ISOFIT itself (calling ray stop --force is not sufficient)

    def run_on_map_geometry(self,
                            enmap_ImageL2: EnMAPL2Product_MapGeo,
                            segmentation: bool = False,
                            n_cores: int = cpu_count()
                            ) -> GeoArray:
        with TemporaryDirectory() as td:
            path_indir = pjoin(td, 'input')
            fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut = self.generate_input_files(enmap_ImageL2, path_indir)

            paths_output = \
                self._run(
                    path_toarad=fp_rad,
                    path_loc=fp_loc,
                    path_obs=fp_obs,
                    path_outdir=pjoin(td, 'output'),
                    path_workdir=pjoin(td, 'workdir'),
                    path_enmap_wavelengths=fp_wvl,
                    path_emulator_basedir=pjoin(Path.home(), '.isofit', 'srtmnet'),
                    path_surface_file=fp_surf,
                    path_lut=fp_lut,
                    aot=enmap_ImageL2.meta.aot,
                    cwv=enmap_ImageL2.meta.water_vapour,
                    segmentation=segmentation,
                    n_cores=n_cores
                )

            # read the AC results back into memory
            boa_rfl = GeoArray(paths_output['estimated_reflectance_file'])
            # state = GeoArray(paths_output['estimated_state_file'])
            # uncert = GeoArray(paths_output['posterior_uncertainty_file'])
            # atm_coef = GeoArray(paths_output['atmospheric_coefficients_file'])  # not always present
            boa_rfl.to_mem()

            return boa_rfl


class LUT_Transformer(object):
    def __init__(self, path_lut: str, sza_scene: float):
        """

        :param path_lut:    atmospheric look-up-table in binary format as created by Luis
        """
        self.p_lut_bin = path_lut
        self._offset = 0
        self.sza_scene = sza_scene

    def read_binary_modtran_lut(self, path_out_nc: str):
        """
        Read MODTRAN® LUT.

        :param path_out_nc: path to output LUT
        :return:         LUT of atmospheric functions, x and y axes grid points, LUT wavelengths
        """
        def read_int16(data: np.ndarray, count):
            val = np.array(data[self._offset:self._offset + count * 2].view('int16'))
            self._offset += count * 2
            return val

        def read_float32(data: np.ndarray, count):
            val = np.array(data[self._offset:self._offset + count * 4].view('f4'))
            self._offset += count * 4
            return val

        with open(self.p_lut_bin, 'rb') as fd:
            data = np.frombuffer(fd.read(), dtype=np.uint8)  # Read all data as bytes

        wvl, vza, sza, alt, aot, raa, cwv = [read_float32(data, count = read_int16(data, count=1)[0]) for _ in range(7)]
        npar1, npar2 = [read_int16(data, count=1)[0] for _ in range(2)]

        # LUT1 containing path radiance
        # NOTE: The rhoatm values do not vary much with changing CVW contents,
        #       therefore the LUT was created with a single CWV value.
        # lut1_axnames = ['VZA', 'SZA', 'ALT', 'AOT', 'RAA', 'CWV', 'WVL', 'PAR1']
        lut1 = read_float32(data, count=10584000).reshape(
            (len(vza), len(sza), len(alt), len(aot), len(raa), 1, len(wvl), npar1))
        # LUT2 containing downwelling direct and diffuse surface solar irradiance, spherical albedo and transmittance
        # NOTE: The values in LUT2 do not vary much with changing RAAs,
        #       therefore the LUT was created with a single RAA value.
        # lut2_axnames = ['VZA', 'SZA', 'ALT', 'AOT', 'RAA', 'CWV', 'WVL', 'PAR2']
        lut2 = read_float32(data, count=42336000).reshape(
            (len(vza), len(sza), len(alt), len(aot), 1, len(cwv), len(wvl), npar2))

        ##############################
        # TRANSFORM TO ISOFIT FORMAT #
        ##############################

        # Get data from LUT1 and process it
        fwhm = np.zeros_like(wvl)
        fwhm[:len(wvl) - 1] = np.diff(wvl)
        fwhm[len(wvl) - 1] = fwhm[len(wvl) - 2]
        cos_sza = np.cos(np.deg2rad(sza))[None, :, None, None, None, None, None]

        # NOTE: ISOFIT expects a LUT with all parameters simulated for multiple CWV values.
        #       Since the binary MODTRAN LUT1 is created with a single CWV value (see above),
        #       LUT1 has to be replicated 7 times for the CWV axis (axis=5). Same applies
        #       to LUT2 and RAA (axis=4).
        _lut1 = np.repeat(lut1, 7, axis=5)
        _lut2 = np.repeat(lut2, 7, axis=4)

        # extract all the parameters from the 2 LUTs
        rhoatm = _lut1[..., 0]  # rhoatm
        edir = _lut2[..., 0]  # transm_down_dir
        edif = _lut2[..., 1]  # transm_down_dif
        sab = _lut2[..., 2]  # spherical albedo
        # The binary MODTRAN LUT has values in _lut2[..., 3], which should be transm_up_dir + transm_up_dif.
        # However, this is already contained in the Edir/Edif variables, thus we need to set it to 1 here.
        var5 = 1  # transm_up_dir + transm_up_dif
        var6 = 0  # not needed ... already 0 and the transm_up_dir includes the transm_up_dif already

        # convert parameters to what ISOFIT expects
        _unit_scalefactor = 1000  # adapt MODTRAN LUT unit to ISOFIT
        rhoatm *= _unit_scalefactor  # path radiance
        sphalb = sab  # atmospheric spherical albedo
        # tg = tdir + tdif
        # transm_down_dir  = tg * cos(sza) * edir /  PI
        # transm_down_dif  = tg * edif /  PI
        _constant = (var5 + var6) / np.pi * _unit_scalefactor
        transm_down_dir = edir * cos_sza * _constant  # direct downward transmittance
        transm_down_dif = edif * _constant  # (diffuse downward transmittance)

        del lut1, lut2, _lut1, _lut2, edir, edif, sab, data

        # extrapolate data at 8 km elevation due to a bug in the MODTRAN-LUT (data at 8km = data at 0 km)
        for arr in [rhoatm, sphalb, transm_down_dir, transm_down_dif]:
            self.extrapolate_8km(arr, alt)

        transm_up_dir = np.zeros_like(rhoatm)  # direct upward transmittance
        transm_up_dif = np.zeros_like(rhoatm)  # diffuse upward transmittance
        thermal_upwelling = np.zeros_like(rhoatm)  # thermal up-welling
        thermal_downwelling = np.zeros_like(rhoatm)  # thermal down-welling

        # write LUT to NetCDF file
        with nc.Dataset(path_out_nc, 'w', format='NETCDF4') as nc_out:
            nc_out.setncattr("RT_mode", "rdn")

            var = nc_out.createVariable('coszen', 'f4', ())
            var.assignValue(np.cos(np.deg2rad(self.sza_scene)))

            # Add dimensions
            for t, v in zip(('H2OSTR', 'AERFRAC_2', 'surface_elevation_km', 'solar_zenith',
                             'observer_zenith', 'relative_azimuth', 'wl'),
                            (cwv, aot, alt, sza, vza, raa, wvl)):
                nc_out.createDimension(t, v.size)

            # Add data
            dims = ('observer_zenith', 'solar_zenith', 'surface_elevation_km',
                    'AERFRAC_2', 'relative_azimuth', 'H2OSTR', 'wl')
            for t, d, v in [
                # 1D data
                ('H2OSTR', ('H2OSTR',), cwv),
                ('AERFRAC_2', ('AERFRAC_2',), aot),
                ('surface_elevation_km', ('surface_elevation_km',), alt),
                ('solar_zenith', ('solar_zenith',), sza),
                ('observer_zenith', ('observer_zenith',), vza),
                ('relative_azimuth', ('relative_azimuth',), raa),
                ('wl', ('wl',), wvl),
                ('fwhm', ('wl',), fwhm),
                ('solar_irr', ('wl',), np.zeros_like(wvl)),  # not used by ISOFIT)
                # 7D data
                ('rhoatm', dims, rhoatm),
                ('sphalb', dims, sphalb),
                ('transm_down_dir', dims, transm_down_dir),
                ('transm_down_dif', dims, transm_down_dif),
                ('transm_up_dir', dims, transm_up_dir),
                ('transm_up_dif', dims, transm_up_dif),
                ('thermal_upwelling', dims, thermal_upwelling),
                ('thermal_downwelling', dims, thermal_downwelling),
            ]:
                var = nc_out.createVariable(t, "f4", d)
                var[:] = v

    @staticmethod
    def extrapolate_8km(var: np.ndarray, alt_grid: np.ndarray):
        """Replace data at 8km altitude by linearly extrapolating with 0.7km and 2.5km as grid points."""
        var[:, :, 3, :, :, :, :] = \
            interp1d(x=alt_grid[1:3],
                     y=var[:, :, 1:3, :, :, :, :],
                     axis=2,
                     fill_value='extrapolate',
                     bounds_error=False)(
                np.array([8.])
            )[:, :, 0, :, :, :, :]

        return var

        # Adjust boundaries
        for arr, dim in zip([vza, sza, alt, aot, raa, cwv], dim_arr):
            arr[0] += 0.0001
            arr[dim - 1] -= 0.0001

        # Extract LUTs
        luts = l0_lut, edir_lut, edif_lut, sab_lut = [
            np.squeeze(lut1[..., 0], axis=5),  # l0 LUT
            np.squeeze(lut2[..., 0], axis=4),  # edir LUT
            np.squeeze(lut2[..., 1], axis=4),  # edif LUT
            np.squeeze(lut2[..., 2], axis=4)   # sab LUT
        ]

        # Define axes
        axes_x = [
            [vza, sza, alt, aot, raa],  # axes x l0
            [vza, sza, alt, aot, cwv]   # axes x e s
        ]
        axes_y = [
            [np.arange(i) for i in luts[0].shape[:-1]],  # axes y l0
            [np.arange(i) for i in luts[0].shape[:-1]]   # axes y e s
        ]

        return luts, axes_x, axes_y, wvl, lut1, lut2, xnodes, 2 ** ndim, ndim, x_cell