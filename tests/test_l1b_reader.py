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
        pass

    def tearDown(self):
        pass

    def test_read_inputdata(self):
        from enpt.io.l1b_reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_ImGeo

        user_config = dict()

        print("Test reading EnMAP Level-1B products")
        for l1b_file in glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip")):
            with tempfile.TemporaryDirectory() as tmpdir:
                print("Tmp dir: %s" % tmpdir)
                with zipfile.ZipFile(l1b_file, "r") as zf:
                    zf.extractall(tmpdir)

                    root_dir = os.path.join(tmpdir, os.listdir(tmpdir)[0])
                    L1B_obj = L1B_Reader(**user_config)\
                        .read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))
                    L1B_obj.vnir.arr.show()

            self.assertIsInstance(L1B_obj, EnMAPL1Product_ImGeo)


if __name__ == "__main__":
    unittest.main()
