# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz.de)
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
from tempfile import TemporaryDirectory, mkdtemp, gettempdir
from typing import Tuple, List
from fnmatch import fnmatch
import os
from os.path import isdir, isfile, join as pjoin, abspath as pabs
from pathlib import Path
from glob import glob
import json
from collections.abc import Mapping
from datetime import datetime
from zipfile import ZipFile
from multiprocessing import cpu_count
import logging
import weakref
from warnings import warn

import numpy as np
from pyproj.crs import CRS
from pandas import DataFrame

from ...utils import EnvContextManager
with EnvContextManager(ISOFIT_DEBUG='0',
                       MKL_NUM_THREADS='1',
                       OMP_NUM_THREADS='1'
                       # ISOFIT_NO_SET_THREADS='1'
                       ):
    import isofit
    from isofit.core.isofit import Isofit
    from isofit.utils import surface_model, analytical_line as run_analytical_line, extractions, segment
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

from ._isofit_downloads import download_isofit_resources
from ._isofit_lut_preparation import LUTTransformer
from ...model.images import EnMAPL2Product_MapGeo
from ...options.config import EnPTConfig, path_enptlib
from ...utils.logging import EnPT_Logger

__author__ = 'Daniel Scheffler'


class IsofitEnMAP(object):
    """Class to perform atmospheric correction of EnMAP data using ISOFIT."""

    def __init__(self,
                 config: EnPTConfig = None,
                 log_level: str = None,
                 tempdir: str = None,
                 ) -> None:
        """Create an instance of IsofitEnMAP.

        :param config:      instance of EnPTConfig
        :param log_level:   logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.cfg = config
        self.log_level = log_level or (config.log_level if config else 'INFO')
        self.logger = self._initialize_logging(logger=None)  # default logger without FileHandler (overridden later)
        self.cpus = config.CPUs if config else cpu_count()
        self._tmpdir = (
            tempdir if tempdir and os.path.exists(tempdir) else
            pjoin(self.cfg.working_dir, 'isofit') if config else
            mkdtemp(dir=tempdir)
        )

        # setup a finalizer that destroys remaining data (directories, etc.) in case of unexpected exit
        self._finalizer = weakref.finalize(self, self._cleanup, self._tmpdir,
                                           warn_message="Implicitly cleaning up {!r}".format(self))

        # leave at least 4 cores free to ensure optimal performance (in case there are more than 6 cores available)
        if cpu_count() > 6 and self.cpus > (cpu_count() - 4):
            self.logger.debug(f"Reduced number of CPU cores to be used for ISOFIT from {self.cpus} to "
                              f"{cpu_count() - 4} to achieve optimal performance without breaking execution.")
            self.cpus = cpu_count() - 4

        os.environ['SIXS_DIR'] = pjoin(Path.home(), '.isofit', 'sixs')
        # os.environ['EMULATOR_PATH'] = '/home/gfz-fe/scheffler/srtmnet/sRTMnet_v120.h5'  # duplicate of emulator_base

        # make sure EnPT data for ISOFIT are downloaded (not contained in EnPT package distribution)
        self.path_surf_spec_zip = pjoin(path_enptlib, 'resources', 'isofit', 'isofit_surface_spectra.zip')
        self.path_lut_zip = pjoin(path_enptlib, 'resources', 'isofit', 'lut.zip')
        self.logger.info("Downloading EnPT-internal resources for ISOFIT...")
        download_isofit_resources(pjoin(path_enptlib, 'resources', 'isofit'), self.logger)

        # make sure ISOFIT's extra-files are downloaded
        self.logger.info('Downloading ISOFIT extra-files...')
        download_data(path=None, tag="latest")
        download_examples(path=None, tag="latest")

    @classmethod
    def _cleanup(cls, tempDir, warn_message):
        """Clean up implicitly (not to be called directly)."""
        if tempDir and os.path.exists(tempDir):
            shutil.rmtree(tempDir)
            warn(warn_message, ResourceWarning)

    def _initialize_logging(self, logger: EnPT_Logger = None):
        """
        Initialize logging and forward the ISOFIT logs to EnPT.

        :param logger: EnPT_Logger instance to be used or None (default) to create a new one.
        """
        # get root logger (used by ISOFIT) and remove all StreamHandlers to avoid duplicated log lines
        # NOTE: ray workers have their own root logger which is NOT captured here (logs are missing)  # FIXME
        root_logger = logging.getLogger()
        root_logger.setLevel(self.log_level)
        for h in root_logger.handlers[:]:
            if isinstance(h, logging.StreamHandler):
                root_logger.removeHandler(h)

        logger = logger or EnPT_Logger('log__EnPT_ISOFIT', fmt_suffix=None, log_level=self.log_level)

        # attach EnPT handlers to the root logger so that Isofit logs go through the EnPT logger
        for handler in logger.handlers:
            root_logger.addHandler(handler)

        self.logger = logger

        return logger

    @staticmethod
    def _build_modtran_template_file(path_emulator_basedir: str,
                                     path_obs: str,
                                     path_loc: str,
                                     path_workdir: str,
                                     enmap_timestamp: str):
        """
        Build a MODTRAN template file based on the given input files and parameters.

        :param path_emulator_basedir:   Path to the ISOFIT emulator base directory.
        :param path_obs:                Path to the observation file.
        :param path_loc:                Path to the location file.
        :param path_workdir:            Path to the working directory.
        :param enmap_timestamp:         Timestamp of the EnMAP data.
        :return: None
        """
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

    def generate_input_files(self, enmap: EnMAPL2Product_MapGeo, path_outdir: str):
        """
        Generate the necessary input files for ISOFIT.

        :param enmap:           EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :param path_outdir:     Output directory path.
        :return: A tuple containing file paths to the generated radiance, location,
                 observation, wavelength, surface, and LUT files.
        """
        fp_rad = self._generate_radiance_file(enmap, path_outdir)
        fp_loc = self._generate_loc_file(enmap, path_outdir)
        fp_obs = self._generate_obs_file(enmap, path_outdir)
        fp_wvl = self._generate_wavelength_file(enmap, path_outdir)
        fp_surf = self._generate_surface_file(fp_wvl, path_outdir)
        fp_lut = self._generate_lut_file(path_outdir, enmap.meta.geom_sun_zenith)  # TODO: set LUT to None in case of 6S

        return fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut

    def _generate_radiance_file(self, enmap: EnMAPL2Product_MapGeo, path_outdir: str):
        """
        Generate the radiance file for ISOFIT.

        :param enmap:   EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :param path_outdir: Output directory path.
        :return: A file path to the generated radiance file.
        """
        self.logger.info("Generating radiance file...")
        # ISOFIT expects radiance in uW/cm²/sr/nm, EnPT provides mW/m²/sr/nm
        # 1000 uW/10000 cm²/sr/nm corresponds to mW/m²/sr/nm
        mask_nodata = ~enmap.data.mask_nodata[:]
        radiance = enmap.data[:] / 10.0
        radiance[mask_nodata] = -9999

        fp_out = pjoin(path_outdir, f"{enmap.meta.scene_basename}_rdn")
        rad = GeoArray(radiance, enmap.data.gt, enmap.data.prj, nodata=-9999)
        rad.meta.band_meta = enmap.data.meta.band_meta
        rad.save(fp_out)

        return fp_out

    def _generate_loc_file(self, enmap: EnMAPL2Product_MapGeo, path_outdir: str):
        """
        Generate the location file for ISOFIT containing longitude, latitude, and elevation data.

        This function computes a grid of coordinates based on the bounding box and grid size of the
        input EnMAP data. It transforms the coordinates to WGS-84 if necessary and combines them with
        elevation data extracted from the DEM. The resulting 3D data array is saved as a location file
        with specified metadata and nodata values.

        :param enmap: EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :param path_outdir: Output directory path where the location file will be saved.
        :return: A file path to the generated location file.
        """
        self.logger.info("Generating location file...")
        xmin, xmax, ymin, ymax = enmap.data.box.boundsMap
        xgsd, ygsd = (enmap.data.xgsd, enmap.data.ygsd)
        x_grid, ygrid = get_coord_grid((xmin, ymax), (xmax, ymin), (xgsd, -ygsd))
        if enmap.data.epsg == 4326:
            lons = x_grid
            lats = ygrid
        else:
            lons, lats = transform_coordArray(
                CRS(enmap.data.epsg).to_wkt(),
                CRS(4326).to_wkt(),
                x_grid, ygrid
            )

        elev = enmap.dem[:]  # FIXME nodata value 0
        assert elev.shape == lons.shape == lats.shape

        loc_data = np.dstack([lons, lats, elev])
        loc_data[~enmap.data.mask_nodata[:]] = -9999

        fp_out = pjoin(path_outdir, f"{enmap.meta.scene_basename}_loc")  # no file extension supported
        GeoArray(loc_data,
                 enmap.data.gt, enmap.data.prj,
                 bandnames=[
                     'Longitude (WGS-84)',
                     'Latitude (WGS-84)',
                     'Elevation (m)'  # used as first guess in case this is defined as state_vector component in config
                 ],
                 nodata=-9999
                 ).save(fp_out)

        return fp_out

    def _generate_obs_file(self, enmap: EnMAPL2Product_MapGeo, path_outdir: str):
        """
        Generate observation file for ISOFIT containing angular, slope, aspect, time, and ESD information.

        :param enmap: EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :param path_outdir: Output directory path where the observation file will be saved.
        :return: A file path to the generated observation file.
        """
        self.logger.info("Generating observation file...")
        path_length = np.full(enmap.data.shape[:2], fill_value=650000)  # from ~650km EnMAP flight height
        vaa = enmap.meta.geom_view_azimuth_array
        vza = enmap.meta.geom_view_zenith_array
        saa = enmap.meta.geom_sun_azimuth_array
        sza = enmap.meta.geom_sun_zenith_array
        phase = self._compute_solar_phase(vaa, vza, saa, sza)
        slope = np.full(enmap.data.shape[:2], fill_value=0)
        aspect = np.zeros(enmap.data.shape[:2])
        cos_i = self._compute_cos_i(saa, sza, slope=90, aspect=0)
        utc = enmap.meta.aqtime_utc_array  # TODO pixel-wise values
        # TODO pixel-wise values
        earth_sun_dist = np.full(enmap.data.shape[:2], fill_value=enmap.meta.earthSunDist)

        obs_data = np.dstack([path_length, vaa, vza, saa, sza, phase, slope, aspect, cos_i, utc, earth_sun_dist])
        obs_data[~enmap.data.mask_nodata[:]] = -9999

        fp_out = pjoin(path_outdir, f"{enmap.meta.scene_basename}_obs")  # no file extension supported
        GeoArray(obs_data,
                 enmap.data.gt, enmap.data.prj,
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

    def _generate_wavelength_file(self, enmap: EnMAPL2Product_MapGeo, path_outdir: str):
        """
        Generates a file containing central wavelength and full width at half maximum (FWHM) information for ISOFIT.

        :param enmap: EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :param path_outdir: Output directory path where the wavelength file will be saved.
        :return: A file path to the generated wavelength file.
        """
        self.logger.info("Generating wavelength file...")
        fp_out = pjoin(path_outdir, 'enmap_wavelength_fwhm.txt')
        wvl = enmap.meta.wvl_center
        fwhm = enmap.meta.fwhm
        df = DataFrame(np.hstack([wvl.reshape(-1, 1), fwhm.reshape(-1, 1)]))
        df.to_csv(fp_out, header=False, sep='\t')

        return fp_out

    def _generate_surface_file(self, path_wavelength_file: str, path_outdir: str):
        """
        Generate a file containing surface coverage type information for ISOFIT.

        :param path_wavelength_file: Path to the file containing central wavelength and bandwidths (FWHM) information.
        :param path_outdir: Output directory path where the surface file will be saved.
        :return: A file path to the generated surface file.
        """
        surf_preset = self.cfg.isofit_surface_category if self.cfg else 'multicomponent_surface'
        fp_out = pjoin(path_outdir, 'surface_enmap.mat')

        if surf_preset == 'multicomponent_surface':
            self.logger.info("Generating surface file for default set of surface coverage types...")
            fp_surfjson = pjoin(path_enptlib, 'options', 'isofit_surface_default.json')
        elif surf_preset == 'ree':
            self.logger.info("Generating surface file with specific optimization for rare earth elements (REE)...")
            fp_surfjson = pjoin(path_enptlib, 'options', 'isofit_surface_20240103_REE.json')
        else:  # 'custom'
            self.logger.info(f"Generating surface file based on user provided input "
                             f"({self.cfg.path_isofit_surface_config})...")
            fp_surfjson = self.cfg.path_isofit_surface_config

        path_surf_spec = (
                self.cfg.path_isofit_surface_priors if self.cfg and self.cfg.path_isofit_surface_priors else
                self.path_surf_spec_zip
        )
        with (ZipFile(path_surf_spec, "r") as zf,
              TemporaryDirectory(dir=self._tmpdir, prefix='surface__') as td):
            zf.extractall(td)
            fp_surfjson_tmp = pjoin(td, os.path.basename(fp_surfjson))
            shutil.copyfile(fp_surfjson, fp_surfjson_tmp)

            # generate surface_enmap.mat within /<self._tmpdir>/surface__tmp*/
            # with multiprocessing in KMeans disabled (avoids a UnicodeDecodeError on Windows and is even faster)
            with EnvContextManager(OMP_NUM_THREADS='1'):
                surface_model(
                    config_path=fp_surfjson_tmp,
                    wavelength_path=path_wavelength_file,
                    output_path=fp_out
                )
            assert os.path.isfile(fp_out)

        return fp_out

    def _generate_lut_file(self, path_outdir: str, sza_scene: float):
        """
        Generate the LUT file for ISOFIT.

        :param path_outdir: Output directory path where the LUT file will be saved.
        :param sza_scene: Solar zenith angle of the scene in degrees.
        :return: A file path to the generated LUT file.
        """
        self.logger.info("Generating LUT file...")
        # TODO: By re-using either the unpacked lut.zip or the LUT_ISOFIT.nc,
        #       the processing time can be reduced by ~20-60 sec.
        fp_out = pjoin(path_outdir, 'EnMAP_LUT_MOD5_ISOFIT_formatted_1nm.nc')

        with (ZipFile(self.path_lut_zip, 'r') as zf,
              TemporaryDirectory(dir=self._tmpdir, prefix='lut__') as td):
            zf.extractall(td)

            LUTTransformer(
                path_lut=os.path.join(td, 'EnMAP_LUT_MOD5_formatted_1nm'),
                sza_scene=sza_scene,
                logger=self.logger
            ).read_binary_modtran_lut(
                path_out_nc=fp_out
            )

            assert os.path.isfile(fp_out)

        return fp_out

    @staticmethod
    def _compute_solar_phase(vaa: np.ndarray, vza: np.ndarray, saa: np.ndarray, sza: np.ndarray):
        """
        Compute the solar phase angle given the following angles in degrees:

        :param vaa: View azimuth angle array in degrees.
        :param vza: View zenith angle array in degrees.
        :param saa: Solar azimuth angle array in degrees.
        :param sza: Solar zenith angle array in degrees.
        :return: The solar phase angle in degrees.
        """
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
    def _compute_cos_i(saa: np.ndarray, sza: np.ndarray, slope: float, aspect: float):
        """
        Compute the cosine of the illumination angle (i) given the following angles in degrees:

        :param saa: Solar azimuth angle array in degrees.
        :param sza: Solar zenith angle array in degrees.
        :param slope: Slope of the terrain in degrees.
        :param aspect: Aspect of the terrain in degrees.
        :return: The cosine of the illumination angle (i).
        """
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
                  n_cores: int = None,
                  presolve: bool = False,
                  empirical_line: bool = False,
                  analytical_line: bool = False,
                  ray_temp_dir=pjoin(gettempdir(), "ray"),
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
        n_cores = n_cores if n_cores is not None else self.cpus
        params = {k: v for k, v in locals().items() if not k.startswith('__') and k != 'self'}

        try:
            self.logger.info("Running Isofit.apply_oe().")
            apply_oe(**params)
            self.logger.info("Finished running Isofit.apply_oe().")

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
            self.logger.info('Stopping ray.')
            import ray
            ray.shutdown()  # FIXME: This should be done by ISOFIT itself (calling ray stop --force is not sufficient)

    # def apply_oe_on_sensor_geometry(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
    #     with TemporaryDirectory() as td:
    #         self._apply_oe()

    def apply_oe_on_map_geometry(self, enmap: EnMAPL2Product_MapGeo):
        """
        Apply the ISOFIT atmospheric correction on EnMAP L2 data in map geometry using Isofit.apply_oe().

        This method prepares necessary input files and executes the atmospheric correction
        using ISOFIT on the given EnMAP Level 2 product with map geometry.

        :param enmap: EnMAPL2Product_MapGeo instance containing the EnMAP Level 2 data.
        :return: A GeoArray object containing the bottom-of-atmosphere (BOA) reflectance.
        """
        self._initialize_logging(enmap.logger)  # use enmap.logger instead if self.logger which has no FileHandler

        path_indir = pjoin(self._tmpdir, 'input')
        fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut = self.generate_input_files(enmap, path_indir)

        os.makedirs(pjoin(self._tmpdir, 'workdir'), exist_ok=True)
        os.makedirs(pjoin(self._tmpdir, 'input'), exist_ok=True)
        os.makedirs(pjoin(self._tmpdir, 'output'), exist_ok=True)

        self._apply_oe(
            input_radiance=fp_rad,
            input_loc=fp_loc,
            input_obs=fp_obs,
            working_directory=pjoin(self._tmpdir, 'workdir'),
            surface_path=fp_surf,
            wavelength_path=fp_wvl,
            log_file=pjoin(self._tmpdir, 'output', 'isofit.log'),
            presolve=True,
            emulator_base=pjoin(Path.home(), '.isofit', 'srtmnet', 'sRTMnet_v120.h5'),
            n_cores=self.cpus,
            prebuilt_lut=fp_lut
        )

        # read the AC results back into memory
        boa_rfl = GeoArray(glob(pjoin(self._tmpdir, 'output', 'estimated_reflectance.bsq'))[0])
        # state = GeoArray(glob(pjoin(self._tmpdir, 'output', 'estimated_state.bsq'))[0])
        # uncert = GeoArray(glob(pjoin(self._tmpdir, 'output', 'posterior_uncertainty.bsq'))[0])
        boa_rfl.to_mem()

        return boa_rfl

    def _run(self,
             path_toarad: str,
             path_loc: str,
             path_obs: str,
             path_outdir: str,
             path_workdir: str,
             path_enmap_wavelengths: str,
             path_surface_file: str,
             path_emulator_basedir: str = None,
             path_lut: str = None,
             aot: float = None,
             cwv: float = None,
             segmentation: bool = False,
             segmentation_size: int = 40,
             n_cores: int = None
             ) -> dict:  # noqa
        """
        Run the atmospheric correction process using Isofit.run().

        This executes the atmospheric correction process using the provided
        paths and parameters. It reads input data, applies the correction, and returns
        the results as a dictionary.

        :param path_toarad:             Path to the TOA radiance file.
        :param path_loc:                Path to the location file.
        :param path_obs:                Path to the observation file.
        :param path_outdir:             Path to the output directory.
        :param path_workdir:            Path to the working directory.
        :param path_enmap_wavelengths:  Path to the EnMAP wavelengths file.
        :param path_surface_file:       Path to the surface file.
        :param path_emulator_basedir:   Path to the emulator base directory.
        :param path_lut:                Path to the lookup table.
        :param aot:                     Aerosol optical thickness value.
        :param cwv:                     Columnar water vapor value.
        :param segmentation:            Flag to enable segmentation.
        :param segmentation_size:       Size of the segmentation.
        :param n_cores:                 Number of cores to use.
        :return: Dictionary containing the results of the atmospheric correction.
        """
        enmap_timestamp = os.path.basename(path_toarad).split('____')[1].split('_')[1]
        path_isocfg_default = pjoin(path_enptlib, 'options', 'isofit_config_default_MOD5.json')
        path_isocfg = pjoin(path_workdir, 'config', 'isofit_config.json')
        path_data = pabs(pjoin(Path.home(), '.isofit', 'data'))
        path_examples = pabs(pjoin(Path.home(), '.isofit', 'examples'))
        path_outdir = path_outdir or pjoin(self._tmpdir, 'output')
        path_workdir = path_workdir or pjoin(self._tmpdir, 'workdir')
        path_logfile = pjoin(path_outdir, f'{enmap_timestamp}_isofit.log')
        n_cores = n_cores if n_cores is not None and n_cores < (cpu_count() - 4) else self.cpus

        use_6s = not path_lut or not isfile(path_lut)
        if use_6s:
            if not path_emulator_basedir or not isdir(path_emulator_basedir):
                raise ValueError(path_emulator_basedir,
                                 "'path_emulator_basedir' does not point to an existing directory.")

        if os.path.isdir(path_workdir):
            if path_lut and path_lut.startswith(path_workdir):
                raise ValueError(path_lut, "The given prebuilt LUT must not be within the given "
                                           "working directory as it is deleted before running ISOFIT.")
            shutil.rmtree(path_workdir)

        for d in [
            path_workdir,
            pjoin(path_workdir, 'config'),
            pjoin(path_workdir, 'lut_full'),  # needed for 6S simulations
            path_outdir
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
                    wavelength_file=pabs(path_enmap_wavelengths),
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
                            # use lut_path if existing, otherwise simulate to lut.nc
                            lut_path=path_lut or pjoin(path_workdir, 'lut_full', 'lut.nc'),
                            sim_path=pjoin(path_workdir, 'lut_full'),
                            template_file=pjoin(path_workdir, 'config', f'{enmap_timestamp}_modtran_tpl.json')
                        ) if use_6s else dict(
                            earth_sun_distance_file=pjoin(path_data, 'earth_sun_distance.txt'),
                            irradiance_file=pjoin(path_examples, '20151026_SantaMonica/data/prism_optimized_irr.dat'),
                            # use lut_path if existing, otherwise simulate to lut.nc
                            lut_path=path_lut
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
                    surface_file=pabs(path_surface_file)
                ),
            ),
            implementation=dict(
                debug_mode=False,
                # debug_mode=True,  # TODO deactivate if done
                n_cores=n_cores,
                ray_temp_dir=pjoin(gettempdir(), "ray"),
            ),
            input=dict(
                measured_radiance_file=pabs(path_toarad),
                loc_file=pabs(path_loc),
                obs_file=pabs(path_obs)
            ),
            output=dict(
                estimated_reflectance_file=pjoin(pabs(path_outdir), f'{enmap_timestamp}_estimated_reflectance.bsq'),
                estimated_state_file=pjoin(pabs(path_outdir), f'{enmap_timestamp}_estimated_state.bsq'),
                posterior_uncertainty_file=pjoin(pabs(path_outdir), f'{enmap_timestamp}_posterior_uncertainty.bsq'),
            )
        )

        def update_nested_dict(dic, u):
            for k, v in u.items():
                if isinstance(v, Mapping):
                    dic[k] = update_nested_dict(dic.get(k, {}), v)
                else:
                    dic[k] = v
            return dic

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
            ray_temp_dir=pjoin(gettempdir(), "ray"),
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
                       not os.path.exists(paths.radiance_working_path):
                        self.logger.info("Running forward segmentation...")
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

                            self.logger.info("Extracting " + outp)
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

                if use_6s:
                    self.logger.info("Building template file for 6S...")
                    self._build_modtran_template_file(
                        path_emulator_basedir,
                        path_obs,
                        path_loc,
                        path_workdir,
                        enmap_timestamp
                    )

                self.logger.info("Running ISOFIT...")
                Isofit(
                    config_file=path_isocfg,
                    level=self.log_level,
                    logfile=path_logfile
                ).run(row_column=None)

                if segmentation:
                    self.logger.info("Running inverse-segmentation through analytical line inference")
                    run_analytical_line(
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
                        atm_file=paths.state_working_path,
                        loglevel=self.log_level,
                        logfile=path_logfile
                    )

                    self.logger.info("ISOFIT finished.")
                    return dict(
                        estimated_reflectance_file=paths.rfl_working_path,
                        estimated_state_file=paths.state_working_path,
                        posterior_uncertainty_file=paths.uncert_working_path
                    )
                else:
                    self.logger.info("ISOFIT finished.")
                    return isocfg['output']

            finally:
                # if os.environ.get('ISOFIT_DEBUG') != '1':
                # FIXME: This should be done by ISOFIT itself (calling ray stop --force is not sufficient)
                try:
                    import ray
                    if ray.is_initialized():
                        self.logger.info('Stopping ray.')
                        ray.shutdown()
                except Exception as e:
                    self.logger.info(f"Ray shutdown failed: {e}")

    def run_on_map_geometry(self,
                            enmap: EnMAPL2Product_MapGeo,
                            segmentation: bool = False,
                            n_cores: int = None
                            ) -> (GeoArray, GeoArray, GeoArray):
        """
        Run ISOFIT atmospheric correction on map geometry.

        This function initializes logging, sets up input files, and runs the ISOFIT
        atmospheric correction process on the provided EnMAP L2 product in map geometry.

        :param enmap:           The EnMAP L2 product in map geometry to process.
        :param segmentation:    Flag to determine if segmentation should be applied.
        :param n_cores:         Number of cores to use during processing.
        :return: A tuple containing the estimated reflectance, atmospheric state, and uncertainty GeoArrays.
        """
        self._initialize_logging(enmap.logger)  # use enmap.logger instead if self.logger which has no FileHandler
        self.logger.info("Initializing ISOFIT run on map geometry...")

        path_indir = pjoin(self._tmpdir, 'input')
        fp_rad, fp_loc, fp_obs, fp_wvl, fp_surf, fp_lut = self.generate_input_files(enmap, path_indir)

        paths_output = \
            self._run(
                path_toarad=fp_rad,
                path_loc=fp_loc,
                path_obs=fp_obs,
                path_outdir=pjoin(self._tmpdir, 'output'),
                path_workdir=pjoin(self._tmpdir, 'workdir'),
                path_enmap_wavelengths=fp_wvl,
                path_emulator_basedir=pjoin(Path.home(), '.isofit', 'srtmnet'),
                path_surface_file=fp_surf,
                path_lut=fp_lut,
                aot=enmap.meta.aot,
                cwv=enmap.meta.water_vapour,
                segmentation=segmentation,
                n_cores=n_cores if n_cores is not None else self.cpus
            )

        # read the AC results back into memory
        boa_rfl = GeoArray(paths_output['estimated_reflectance_file'])
        state = GeoArray(paths_output['estimated_state_file'])
        uncert = GeoArray(paths_output['posterior_uncertainty_file'])
        # atm_coef = GeoArray(paths_output['atmospheric_coefficients_file'])  # not always present
        boa_rfl.to_mem()
        state.to_mem()
        uncert.to_mem()

        return boa_rfl, state, uncert
