#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2022 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
test_l1b_reader
---------------

Tests for `l1b_reader` module.
"""

import unittest
import os
from os import path
import tempfile
import zipfile
import shutil

from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL1Product_SensorGeo
from enpt.options.config import EnPTConfig, config_for_testing, config_for_testing_dlr

__author__ = 'Daniel Scheffler'


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class."""

    def setUp(self):
        self.config = EnPTConfig(**config_for_testing)

        # don't drop bands - otherwise we can't run write-read-tests as the writer does not include the full bandlist
        self.config.drop_bad_bands = False

        self.pathList_testimages = [self.config.path_l1b_enmap_image,
                                    self.config.path_l1b_enmap_image_gapfill]
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)
        os.makedirs(self.config.output_dir, exist_ok=True)

        # unzip both test images in dummy format
        for l1b_file in self.pathList_testimages:
            with zipfile.ZipFile(l1b_file, "r") as zf:
                zf.extractall(self.tmpdir)

        self.testproducts = [os.path.join(self.tmpdir, os.path.basename(self.pathList_testimages[i]).split(".zip")[0])
                             for i in range(len(self.pathList_testimages))]

        self.RD = L1B_Reader(config=self.config)

    def tearDown(self):
        if os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        if os.path.isdir(self.config.output_dir):
            shutil.rmtree(self.config.output_dir)

    def test_read_and_save_single_image_no_snr(self):
        """Test to read test image 1, save it and read the saved result again - without SNR."""
        with tempfile.TemporaryDirectory(dir=self.config.output_dir) as tempdir:
            # read
            L1_obj = self.RD.read_inputdata(self.testproducts[0], compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNone(L1_obj.swir.detector_meta.snr)

            # save
            root_dir_written_L1_data = L1_obj.save(tempdir)

            # read saved result
            L1_obj = self.RD.read_inputdata(root_dir_written_L1_data, compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

    def test_read_and_save_single_image_with_snr(self):
        """Test to read test image 1, save it and read the saved result again - with SNR."""
        with tempfile.TemporaryDirectory(dir=self.config.output_dir) as tempdir:
            # read
            L1_obj = self.RD.read_inputdata(self.testproducts[0], compute_snr=True)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNotNone(L1_obj.swir.detector_meta.snr)

            # save
            root_dir_written_L1_data = L1_obj.save(tempdir)

            # read saved result
            L1_obj = self.RD.read_inputdata(root_dir_written_L1_data, compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

    def _test_append_n_lines(self, *reader_args, **reader_kwargs):
        with tempfile.TemporaryDirectory(dir=self.config.output_dir) as tempdir:
            # read images and test append method
            L1_obj = self.RD.read_inputdata(*reader_args, **reader_kwargs)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            if reader_kwargs['compute_snr']:
                self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
            else:
                self.assertIsNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNone(L1_obj.swir.detector_meta.snr)

            # save
            root_dir_written_L1_data = L1_obj.save(path.join(tempdir))

            # read saved result
            L1_obj = self.RD.read_inputdata(root_dir_written_L1_data, compute_snr=reader_kwargs['compute_snr'])
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

    def _test_append_n_lines_allimagecombinations_withwithoutSNR(self, n_lines):
        # append second test image to first (with and without SNR)
        self._test_append_n_lines(self.testproducts[0], self.testproducts[1], n_line_ext=n_lines, compute_snr=False)
        self._test_append_n_lines(self.testproducts[0], self.testproducts[1], n_line_ext=n_lines, compute_snr=True)

        # append first test image to second (with and without SNR)
        self._test_append_n_lines(self.testproducts[1], self.testproducts[0], n_line_ext=n_lines, compute_snr=False)
        self._test_append_n_lines(self.testproducts[1], self.testproducts[0], n_line_ext=n_lines, compute_snr=True)

    def test_append_all_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=None)

    def test_append_10_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=10)

    def test_append_50_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=50)

    def test_append_80_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=80)

    def test_append_100_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=100)

    def test_append_150_lines(self):
        self._test_append_n_lines_allimagecombinations_withwithoutSNR(n_lines=150)


class Test_L1B_Reader_DLR(unittest.TestCase):
    """Tests for L1B_Reader class."""

    def setUp(self):
        self.config = EnPTConfig(**config_for_testing_dlr)
        self.pathList_testimages = [self.config.path_l1b_enmap_image,
                                    self.config.path_l1b_enmap_image_gapfill]
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        shutil.rmtree(self.config.output_dir)

    def test_read_inputdata(self):
        with zipfile.ZipFile(self.pathList_testimages[0], "r") as zf:
            zf.extractall(self.tmpdir)

        RD = L1B_Reader(config=self.config)

        L1_obj = RD.read_inputdata(self.tmpdir, compute_snr=False)
        L1_obj.save(path.join(self.config.output_dir, "no_snr"))

    def test_read_inputdata_dont_drop_bad_bands(self):
        with zipfile.ZipFile(self.pathList_testimages[0], "r") as zf:
            zf.extractall(self.tmpdir)

        cfg = self.config
        cfg.drop_bad_bands = False
        RD = L1B_Reader(config=cfg)

        L1_obj = RD.read_inputdata(self.tmpdir, compute_snr=False)
        self.assertEqual(L1_obj.swir.detector_meta.nwvl, 130)


if __name__ == '__main__':
    import pytest
    pytest.main()
