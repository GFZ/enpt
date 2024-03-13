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

    def test_run_isofit(self):
        with TemporaryDirectory() as td:
            from geoarray import GeoArray
            gA = GeoArray('/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_rdn_sub.bil')
            gA.show()

            IsofitEnMAP()._apply_oe(
                # input_radiance='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_rdn_sub.bil',
                input_radiance='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/ENMAP01-____L1X-DT000000XXXX_20220712T000000Z_00x_VXXXXXX_XXXXXXTXXXXXXZ.bil',
                input_loc='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_loc_sub.bil',
                input_obs='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/data_in/emp20220712t184754_obs_sub.bil',
                working_directory=td,
                surface_path='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/surface/surface_20221020_EnMAP.mat',
                wavelength_path='/home/gfz-fe/scheffler/temp/EnPT/isofit_implementation/sensor_new/enmap_wavelengths.txt',
                log_file=os.path.join(td, 'isofit.log'),
                # lut_config_file=os.path.join(td, 'lut_config_file.nc'),
            )


if __name__ == '__main__':
    import pytest
    pytest.main()
