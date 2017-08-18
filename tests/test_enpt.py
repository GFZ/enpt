#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_enpt
----------------------------------

Tests for `enpt` module.
"""

import unittest
from glob import glob
from os import path
import tempfile
import zipfile
from datetime import datetime


class TestEnPT(unittest.TestCase):
    """EnPT integration tests."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_IO_EnMAPL1BProduct(self):
        from enpt.io import EnMAPL1BProduct
        print("Test reading EnMAP Level-1B products")
        for l1b_file in glob(path.join(path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip")):
            print("File: %s" % l1b_file)
            with tempfile.TemporaryDirectory() as tmpdir:
                print("Tmp dir: %s" % tmpdir)
                with zipfile.ZipFile(l1b_file, "r") as zf:
                    zf.extractall(tmpdir)
                    l1b_header_fn = glob(path.join(tmpdir, "*", "*_header.xml"))[0]
                    EnMAPL1BProduct(header_fn=l1b_header_fn,
                                    observation_time=datetime(2015, 12, 7, 10))


if __name__ == "__main__":
    unittest.main()
