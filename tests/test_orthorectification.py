#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tools - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
#
# This software was developed within the context of the EnMAP project supported
# by the DLR Space Administration with funds of the German Federal Ministry of
# Economic Affairs and Energy (on the basis of a decision by the German Bundestag:
# 50 EE 1529) and contributions from DLR, GFZ and OHB System AG.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
test_orthorectification
-----------------------

Tests for `processors.orthorectification.orthorectification` module.
"""

import os
from unittest import TestCase
from zipfile import ZipFile
import tempfile
import shutil

from enpt.processors.orthorectification import Orthorectifier
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL2Product_MapGeo


class Test_Orthorectifier(TestCase):
    def setUp(self):
        self.config = EnPTConfig(**config_for_testing)  # FIXME still the Alpine dataset

        # create a temporary directory
        # NOTE: This must exist during the whole runtime of Test_Orthorectifier, otherwise
        #       Orthorectifier.run_transformation will fail to read some files.
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

        # get lons / lats
        with ZipFile(self.config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(self.tmpdir)
            self.L1_obj = L1B_Reader(config=self.config).read_inputdata(
                root_dir_main=os.path.join(self.tmpdir,
                                           os.path.splitext(os.path.basename(self.config.path_l1b_enmap_image))[0]),
                compute_snr=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_run_transformation(self):
        # FIXME td does not exist here anymore
        OR = Orthorectifier(config=self.config)
        L2_obj = OR.run_transformation(self.L1_obj)

        self.assertIsInstance(L2_obj, EnMAPL2Product_MapGeo)
        self.assertTrue(L2_obj.data.is_map_geo)
