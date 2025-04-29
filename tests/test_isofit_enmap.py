#!/usr/bin/env python
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

"""
test_isofit_enmap
-----------------

Tests for `atmospheric_correction._isofit_enmap` module.
"""

import unittest
import os
from os.path import join as pjoin
from pathlib import Path
import tempfile
import shutil
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from importlib.util import find_spec
from fnmatch import filter as fnfilter
from glob import glob

import numpy as np
import pytest
from geoarray import GeoArray

from enpt.io.reader import L1B_Reader
from enpt.processors.orthorectification import Orthorectifier
from enpt.options.config import EnPTConfig, config_for_testing, config_for_testing_dlr
from enpt.processors.atmospheric_correction._isofit_enmap import IsofitEnMAP, LUTTransformer

__author__ = 'Daniel Scheffler'

path_enptlib = os.path.dirname(find_spec("enpt").origin)
dir_isofit_data = os.path.join(path_enptlib, '..', 'tests', 'data', 'isofit')


class Test_ISOFIT_EnMAP(unittest.TestCase):
    """Tests for L1B_Reader class."""

    def setUp(self):
        self.config = EnPTConfig(**config_for_testing)

        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        if os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        if os.path.isdir(self.config.output_dir):
            shutil.rmtree(self.config.output_dir)

    @staticmethod
    def _get_enmap_l2a_obj():
        cfg_dict = dict(config_for_testing_dlr, **dict(target_projection_type='UTM'))
        cfg = EnPTConfig(**cfg_dict)

        with ZipFile(cfg.path_l1b_enmap_image, "r") as zf, \
                TemporaryDirectory(cfg.working_dir) as td:
            zf.extractall(td)
            L1_obj = L1B_Reader(config=cfg).read_inputdata(
                root_dir_main=td,
                compute_snr=False)

            L2_obj = Orthorectifier(config=cfg).run_transformation(L1_obj)

        return L2_obj

    @pytest.mark.skip(reason="too slow for running in CI")
    def test_apply_oe__direct_call(self):
        # from geoarray import GeoArray
        # gA = GeoArray('/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_rdn_sub.bil')
        # gA.show()

        from time import time
        t0 = time()

        IsofitEnMAP()._apply_oe(
            # input_radiance='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_rdn_sub.bil',
            # input_radiance='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/ENMAP01-____L1X-DT000000XXXX_20220712T000000Z_00x_VXXXXXX_XXXXXXTXXXXXXZ.bil',
            # input_loc='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_loc_sub.bil',
            # input_obs='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub.bil',

            input_radiance='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/ENMAP01-____L1X-DT000000XXXX_20220712T000000Z_00x_VXXXXXX_XXXXXXTXXXXXXZ__subX0-10Y0-10.bsq',
            input_loc='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_loc_sub__subX0-10Y0-10.bsq',
            input_obs='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub__subX0-10Y0-10.bsq',
            working_directory='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_out__subX0-10Y0-10/',
            surface_path='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/surface/surface_20221020_EnMAP.mat',
            wavelength_path='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/sensor_new/enmap_wavelengths.txt',
            log_file='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_out/isofit.log',
            presolve=True,
            emulator_base=pjoin(Path.home(), '.isofit', 'srtmnet', 'sRTMnet_v120.h5'),
            n_cores=30
        )

        t = time() - t0
        print(t / 60)
        a = 1

    @pytest.mark.skip(reason="too slow for running in CI")
    def test_apply_oe_on_map_geometry(self):
        IsofitEnMAP().apply_oe_on_map_geometry(self._get_enmap_l2a_obj())

    def test__run__backtransformed_l2_spectra_lut_mod5(self):
        """4x3 spectra within 10x10 image. Lines: different forward sim.: #1: UVEG LUT; #2: Luis LUT, #3: MODTRAN."""
        with (TemporaryDirectory() as td,
              ZipFile(pjoin(dir_isofit_data, 'isofit_testdata.zip'), "r") as zf):

            p_extr = pjoin(td, 'extracted')
            zf.extractall(p_extr)
            files = glob(pjoin(p_extr, 'backtransformed_l2_spectra_v9', '*'))

            IR = IsofitEnMAP()
            IR._run(
                path_toarad=fnfilter(files, '*ENMAP*__subX800-810Y370-380')[0],
                path_loc=fnfilter(files, '*ENMAP*loc')[0],
                path_obs=fnfilter(files, '*ENMAP*obs*v2')[0],
                path_outdir=pjoin(td, 'output/'),
                path_workdir=pjoin(td, 'workdir/'),
                path_enmap_wavelengths=pjoin(p_extr, 'enmap_wavelengths.txt'),
                # path_emulator_basedir=pjoin(Path.home(), '.isofit', 'srtmnet'),
                path_surface_file=pjoin(p_extr, 'surface_20221020_EnMAP.mat'),
                path_lut=IR._generate_lut_file(td, 45),
                segmentation=False,
            )
            a = 1

    @pytest.mark.skip(reason="too slow for running in CI")
    def test__run__backtransformed_l2_spectra_6s(self):
        """4x3 spectra within 10x10 image. Lines: different forward sim.: #1: UVEG LUT; #2: Luis LUT, #3: MODTRAN."""
        with (TemporaryDirectory() as td,
              ZipFile(pjoin(dir_isofit_data, 'isofit_testdata.zip'), "r") as zf):

            p_extr = pjoin(td, 'extracted')
            zf.extractall(p_extr)
            files = glob(pjoin(p_extr, 'backtransformed_l2_spectra_v9', '*'))

            IR = IsofitEnMAP()
            IR._run(
                path_toarad=fnfilter(files, '*ENMAP*__subX800-810Y370-380')[0],
                path_loc=fnfilter(files, '*ENMAP*loc')[0],
                path_obs=fnfilter(files, '*ENMAP*obs*v2')[0],
                path_outdir=pjoin(td, 'output/'),
                path_workdir=pjoin(td, 'workdir/'),
                path_enmap_wavelengths=pjoin(p_extr, 'enmap_wavelengths.txt'),
                path_emulator_basedir=pjoin(Path.home(), '.isofit', 'srtmnet'),
                path_surface_file=pjoin(p_extr, 'surface_20221020_EnMAP.mat'),
                path_lut=None,
                segmentation=False,
            )

    @pytest.mark.skip(reason="too slow for running in CI")
    def test_run_on_map_geometry(self):
        IsofitEnMAP().run_on_map_geometry(self._get_enmap_l2a_obj(), segmentation=True)

    @pytest.mark.skip(reason="too slow for running in CI")
    def test_run_on_map_geometry_no_segmentation(self):
        IsofitEnMAP().run_on_map_geometry(self._get_enmap_l2a_obj(), segmentation=False)

    def test_generate_input_files(self):
        with TemporaryDirectory() as td:
            IsofitEnMAP().generate_input_files(self._get_enmap_l2a_obj(), td)

    def test__compute_solar_phase(self):
        with (TemporaryDirectory() as td,
              ZipFile(pjoin(dir_isofit_data, 'isofit_testdata.zip'), "r") as zf):

            p_extr = pjoin(td, 'extracted')
            zf.extractall(p_extr)
            files = glob(pjoin(p_extr, 'c_translation_test', '*'))
            gA = GeoArray(fnfilter(files, '*ENMAP*__subX800-810Y370-380')[0])

            vaa = gA[:, :, 1]
            vza = gA[:, :, 2]
            saa = gA[:, :, 3]
            sza = gA[:, :, 4]
            phase_karl = gA[:, :, 5]

        phase_py = IsofitEnMAP()._compute_solar_phase(vaa, vza, saa, sza)
        assert np.allclose(phase_karl, phase_py)

    def test__compute_cos_i(self):
        with (TemporaryDirectory() as td,
              ZipFile(pjoin(dir_isofit_data, 'isofit_testdata.zip'), "r") as zf):

            p_extr = pjoin(td, 'extracted')
            zf.extractall(p_extr)
            files = glob(pjoin(p_extr, 'c_translation_test', '*'))
            gA = GeoArray(fnfilter(files, '*ENMAP*__subX800-810Y370-380')[0])

            saa = gA[:, :, 3]
            sza = gA[:, :, 4]
            cos_i_karl = gA[:, :, 8]

        cos_i_py = IsofitEnMAP()._compute_cos_i(saa, sza, slope=90, aspect=0)
        assert np.allclose(cos_i_karl, cos_i_py)


class Test_LUT_Transformer(unittest.TestCase):
    """Tests for L1B_Reader class."""

    def setUp(self):
        self.p_lut_bin = '/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/SCAPE_M/EnMAP_LUT_MOD5_formatted_1nm'

    def test_modtran_lut_to_netcdf(self):
        LUTTransformer(self.p_lut_bin, sza_scene=40).read_binary_modtran_lut(self.p_nc_enpt)


if __name__ == '__main__':
    pytest.main()
