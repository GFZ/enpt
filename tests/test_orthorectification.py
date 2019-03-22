#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
