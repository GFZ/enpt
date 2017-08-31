#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_l1b_reader
----------------------------------

Tests for `l1b_reader` module.
"""

import unittest
from glob import glob
import os
import tempfile
import zipfile
from datetime import datetime


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.pathList_testimages = glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip"))
        self.l1b_snr_file = glob(os.path.join(os.path.dirname(__file__),
                                              "data", "EnMAP_Sensor", "EnMAP_Level_1B_SNR.zip"))[0]
        self.user_config = dict()

    def tearDown(self):
        pass

    def test_read_inputdata(self):
        from enpt.io.reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_SensorGeo

        print("Test reading EnMAP Level-1B products")
        for l1b_file in self.pathList_testimages:
            with tempfile.TemporaryDirectory() as tmpdir:
                print("Tmp dir: %s" % tmpdir)
                with zipfile.ZipFile(self.l1b_snr_file, "r") as zf:
                    zf.extractall(tmpdir)

                    with zipfile.ZipFile(l1b_file, "r") as zf:
                        zf.extractall(tmpdir)

                        root_dir = os.path.join(tmpdir, os.path.basename(l1b_file).split(".zip")[0])
                        L1_obj_no_snr = L1B_Reader(**self.user_config)\
                            .read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))

                        L1_obj_with_snr = L1B_Reader(**self.user_config) \
                            .read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10),
                                            snr_vnir=os.path.join(tmpdir, "SNR_D1.hdr"),
                                            snr_swir=os.path.join(tmpdir, "SNR_D2.hdr"))

            self.assertIsInstance(L1_obj_no_snr, EnMAPL1Product_SensorGeo)
            self.assertIsInstance(L1_obj_with_snr, EnMAPL1Product_SensorGeo)


if __name__ == "__main__":
    unittest.main()
