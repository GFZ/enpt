#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2023 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# Details can be found here: https://github.com/tqdm/tqdm/blob/main/LICENCE.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
test_spatial_optimization
-------------------------

Tests for `processors.spatial_optimization.spatial_optimization` module.
"""

from unittest import TestCase
from zipfile import ZipFile
import tempfile
import shutil
import numpy as np

from enpt.processors.spatial_optimization import Spatial_Optimizer
from enpt.options.config import config_for_testing_dlr, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.model.images.images_sensorgeo import EnMAPL1Product_SensorGeo

__author__ = 'Daniel Scheffler'


class Test_Spatial_Optimizer(TestCase):
    def setUp(self):
        self.config = EnPTConfig(**config_for_testing_dlr)

        # create a temporary directory
        # NOTE: This must exist during the whole runtime of Test_Spatial_Optimizer, otherwise
        #       Spatial_Optimizer.optimize_geolayer will fail to read some files.
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

        # get EnMAP L1 object in sensor geometry
        with ZipFile(self.config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(self.tmpdir)
            self.L1_obj = L1B_Reader(config=self.config).read_inputdata(
                root_dir_main=self.tmpdir,
                compute_snr=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_optimize_geolayer(self):
        SO = Spatial_Optimizer(config=self.config)
        L1_obj = SO.optimize_geolayer(self.L1_obj)

        self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
        self.assertNotEqual(np.mean(L1_obj.meta.vnir.lons), 0)
        self.assertNotEqual(np.std(L1_obj.meta.vnir.lons), 0)
        self.assertNotEqual(np.mean(L1_obj.meta.vnir.lats), 0)
        self.assertNotEqual(np.std(L1_obj.meta.vnir.lats), 0)


if __name__ == '__main__':
    import pytest
    pytest.main()
