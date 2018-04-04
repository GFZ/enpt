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


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.pathList_testimages = glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip"))
        self.config = EnPTConfig()

        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_read_and_save_inputdata(self):
        from enpt.io.reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_SensorGeo

        print("Test reading EnMAP Level-1B products")

        rd = L1B_Reader(config=self.config)

        for l1b_file in self.pathList_testimages:
            print("Tmp dir: %s" % self.tmpdir)
            with zipfile.ZipFile(self.config.path_l1b_snr_model, "r") as zf:
                zf.extractall(self.tmpdir)

                with zipfile.ZipFile(l1b_file, "r") as zf:
                    zf.extractall(self.tmpdir)
                    # without snr
                    root_dir = os.path.join(self.tmpdir, os.path.basename(l1b_file).split(".zip")[0])
                    L1_obj = rd.read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))
                    self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                    root_dir = L1_obj.save(path.join(self.tmpdir, "_no_snr"))
                    L1_obj = rd.read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))
                    self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                    # with snr
                    L1_obj = rd.read_inputdata(
                        root_dir, observation_time=datetime(2015, 12, 7, 10),
                        snr_vnir=os.path.join(self.tmpdir, "SNR_D1.hdr"),
                        snr_swir=os.path.join(self.tmpdir, "SNR_D2.hdr"))
                    self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                    root_dir = L1_obj.save(path.join(self.tmpdir), "with_snr")
                    L1_obj = rd.read_inputdata(
                        root_dir, observation_time=datetime(2015, 12, 7, 10),
                        snr_vnir=os.path.join(self.tmpdir, "SNR_D1.hdr"),
                        snr_swir=os.path.join(self.tmpdir, "SNR_D2.hdr"))
                    self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)


if __name__ == "__main__":
    unittest.main()
