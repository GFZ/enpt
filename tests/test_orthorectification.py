#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_orthorectification
-----------------------

Tests for `processors.orthorectification.orthorectification` module.
"""

import os
from unittest import TestCase
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from enpt.processors.orthorectification import Orthorectifier
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL2Product_MapGeo


class Test_Orthorectifier(TestCase):
    def setUp(self):
        config = EnPTConfig(**config_for_testing)

        # get lons / lats
        with TemporaryDirectory() as td, ZipFile(config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(td)
            self.L1_obj = L1B_Reader(config=config).read_inputdata(
                root_dir_main=os.path.join(td, os.path.splitext(os.path.basename(config.path_l1b_enmap_image))[0]),
                lon_lat_smpl=(1000, 100),
                compute_snr=False)

    def test_run_transformation(self):
        OR = Orthorectifier(config=config_for_testing)
        L2_obj = OR.run_transformation(self.L1_obj)

        self.assertIsInstance(L2_obj, EnMAPL2Product_MapGeo)
        self.assertTrue(L2_obj.data.is_map_geo)
