#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018–2025 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz.de)
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
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""
test_l1b_reader
---------------

Tests for `l1b_reader` module.
"""

import unittest
import os
from os import path
from pathlib import Path
import tempfile
import zipfile
import shutil
from glob import glob

import numpy as np

from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL1Product_SensorGeo
from enpt.options.config import EnPTConfig, path_enptlib

__author__ = 'Daniel Scheffler'


class Test_L1B_Reader_DLR(unittest.TestCase):
    """Tests for L1B_Reader class."""

    @classmethod
    def setUpClass(cls):
        path_l1b_testimages = (Path(path_enptlib) / ".." / "tests" / "data" / "EnMAP_Level_1B").resolve()
        cls.config = EnPTConfig(
            path_l1b_enmap_image=str(
                path_l1b_testimages / "ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z"
                                      "__rows700-799.zip"
            ),
            path_l1b_enmap_image_gapfill=str(
                path_l1b_testimages / "ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z"
                                      "__rows800-899.zip"
            ),
            output_dir=str((Path(path_enptlib) / ".." / "tests" / "data" / "test_outputs" / 'test_reader').resolve())
        )
        cls.config.drop_bad_bands = False  # otherwise the read/write/read tests will fail
        cls.pathList_testimages = [cls.config.path_l1b_enmap_image,
                                   cls.config.path_l1b_enmap_image_gapfill]
        cls.tmpdir = tempfile.mkdtemp(dir=cls.config.working_dir)
        os.makedirs(cls.config.output_dir, exist_ok=True)

        # unzip both test images
        for l1b_file in cls.pathList_testimages:
            with zipfile.ZipFile(l1b_file, "r") as zf:
                zf.extractall(Path(cls.tmpdir) / Path(l1b_file).stem)

        cls.testproducts = glob(os.path.join(cls.tmpdir, '*'))
        cls.RD = L1B_Reader(config=cls.config)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmpdir)
        shutil.rmtree(cls.config.output_dir)

    def test_read_inputdata_dont_drop_bad_bands(self):
        L1_obj = self.RD.read_inputdata(self.testproducts[0], compute_snr=False)
        assert L1_obj.swir.detector_meta.nwvl == 130

    def _test_read_and_save_single_image(self, compute_snr: bool):
        with tempfile.TemporaryDirectory(dir=self.config.output_dir) as tempdir:
            # read
            L1_obj = self.RD.read_inputdata(self.testproducts[0], compute_snr=compute_snr)
            assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)
            if compute_snr:
                assert isinstance(L1_obj.vnir.detector_meta.snr, np.ndarray)
                assert isinstance(L1_obj.swir.detector_meta.snr, np.ndarray)
                assert L1_obj.vnir.detector_meta.snr.shape == L1_obj.vnir.data.shape
                assert L1_obj.swir.detector_meta.snr.shape == L1_obj.swir.data.shape
            else:
                assert L1_obj.vnir.detector_meta.snr is None
                assert L1_obj.swir.detector_meta.snr is None

            # save
            L1_obj.save(tempdir)
            root_dir_written_L1_data = path.join(tempdir, L1_obj.meta.scene_basename)

            # read saved result
            L1_obj = self.RD.read_inputdata(root_dir_written_L1_data, compute_snr=compute_snr)
            assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)

    def test_read_and_save_single_image_no_snr(self):
        """Test to read test image 1, save it and read the saved result again - without SNR."""
        self._test_read_and_save_single_image(compute_snr=False)

    def test_read_and_save_single_image_with_snr(self):
        """Test to read test image 1, save it and read the saved result again - with SNR."""
        self._test_read_and_save_single_image(compute_snr=True)

    def _test_append_n_lines(self, *reader_args, **reader_kwargs):
        with tempfile.TemporaryDirectory(dir=self.config.output_dir) as tempdir:
            # read images and test append method
            L1_obj = self.RD.read_inputdata(*reader_args, **reader_kwargs)
            assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)
            if reader_kwargs['compute_snr']:
                assert L1_obj.vnir.detector_meta.snr is not None
                assert L1_obj.swir.detector_meta.snr is not None
            else:
                assert L1_obj.vnir.detector_meta.snr is None
                assert L1_obj.swir.detector_meta.snr is None

            # save
            L1_obj.save(tempdir)
            root_dir_written_L1_data = path.join(tempdir, L1_obj.meta.scene_basename)

            # read saved result
            L1_obj = self.RD.read_inputdata(root_dir_written_L1_data,
                                            compute_snr=reader_kwargs['compute_snr'])
            assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)

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


if __name__ == '__main__':
    import pytest
    pytest.main()
