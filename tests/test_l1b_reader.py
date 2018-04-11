#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_l1b_reader
---------------

Tests for `l1b_reader` module.
"""

import unittest
from glob import glob
import os
from os import path
import tempfile
import zipfile
from datetime import datetime
import shutil

from enpt.options.config import EnPTConfig
from . import config_for_testing


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.pathList_testimages = glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip"))
        self.config = EnPTConfig(**config_for_testing)
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        shutil.rmtree(self.config.output_dir)

    def test_read_and_save_inputdata(self):
        from enpt.io.reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_SensorGeo

        print("Test reading EnMAP Level-1B products")

        rd = L1B_Reader(config=self.config)

        for l1b_file in self.pathList_testimages:
            print("Tmp dir: %s" % self.tmpdir)
            with zipfile.ZipFile(l1b_file, "r") as zf:
                zf.extractall(self.tmpdir)

                root_dir = os.path.join(self.tmpdir, os.path.basename(l1b_file).split(".zip")[0])

                ###############
                # without snr #
                ###############

                # read and write L1 data
                L1_obj = rd.read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10), compute_snr=False)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNone(L1_obj.swir.detector_meta.snr)
                root_dir_written_L1_data = L1_obj.save(path.join(self.tmpdir, "no_snr"))

                # read self written L1 data
                L1_obj = rd.read_inputdata(root_dir_written_L1_data, observation_time=datetime(2015, 12, 7, 10),
                                           compute_snr=False)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                ############
                # with snr #
                ############

                # read and write L1 data
                L1_obj = rd.read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
                root_dir_written_L1_data = L1_obj.save(path.join(self.config.output_dir, "with_snr"))

                # read self written L1 data
                L1_obj = rd.read_inputdata(root_dir_written_L1_data, observation_time=datetime(2015, 12, 7, 10))
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)


if __name__ == "__main__":
    unittest.main()
