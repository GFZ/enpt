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
---------------

Tests for `atmospheric_correction._isofit_enmap` module.
"""

import unittest
import os
import tempfile
import shutil
from tempfile import TemporaryDirectory

import numpy as np

from enpt.options.config import EnPTConfig, config_for_testing, config_for_testing_dlr
from enpt.processors.atmospheric_correction._isofit_enmap import IsofitEnMAP

__author__ = 'Daniel Scheffler'


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
        from zipfile import ZipFile
        from enpt.io.reader import L1B_Reader
        from enpt.processors.orthorectification import Orthorectifier
        from enpt.options.config import config_for_testing_dlr
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
            emulator_base='/home/gfz-fe/scheffler/sRTMnet_v100/sRTMnet_v100',  # FIXME why not /home/gfz-fe/scheffler/sRTMnet_v100/
            n_cores=30
        )

        t = time() - t0
        print(t / 60)
        a = 1

    def test_apply_oe_on_map_geometry(self):
        IsofitEnMAP().apply_oe_on_map_geometry(self._get_enmap_l2a_obj())

    def test__run(self):
        IsofitEnMAP()._run(
            path_toarad='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/ENMAP01-____L1X-DT000000XXXX_20220712T000000Z_00x_VXXXXXX_XXXXXXTXXXXXXZ__subX0-10Y0-10.bsq',
            path_loc='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_loc_sub__subX0-10Y0-10.bsq',
            path_obs='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub__subX0-10Y0-10.bsq',
            path_outdir='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/core_run/output/',
            path_workdir='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/core_run/workdir/',
            path_enmap_wavelengths='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/sensor_new/enmap_wavelengths.txt',
            path_emulator_basedir='/home/gfz-fe/scheffler/sRTMnet_v100/sRTMnet_v100',
            path_surface_file='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/surface/surface_20221020_EnMAP.mat'
        )

    def test_run_on_map_geometry(self):
        IsofitEnMAP().run_on_map_geometry(self._get_enmap_l2a_obj())

    def test_generate_input_files(self):
        with TemporaryDirectory() as td:
            IsofitEnMAP().generate_input_files(self._get_enmap_l2a_obj(), td)

    def test__compute_solar_phase(self):
        from geoarray import GeoArray
        gA = GeoArray('/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub__subX0-10Y0-10.bsq')
        vaa = gA[:, :, 1]
        vza = gA[:, :, 2]
        saa = gA[:, :, 3]
        sza = gA[:, :, 4]
        phase_karl = gA[:, :, 5]

        phase_py = IsofitEnMAP()._compute_solar_phase(vaa, vza, saa, sza)
        assert np.allclose(phase_karl, phase_py)

    def test__compute_cos_i(self):
        from geoarray import GeoArray
        gA = GeoArray('/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub__subX0-10Y0-10.bsq')
        saa = gA[:, :, 3]
        sza = gA[:, :, 4]
        cos_i_karl = gA[:, :, 8]

        cos_i_py = IsofitEnMAP()._compute_cos_i(saa, sza, slope=90, aspect=0)
        assert np.allclose(cos_i_karl, cos_i_py)


if __name__ == '__main__':
    import pytest
    pytest.main()
