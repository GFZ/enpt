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
test_orthorectification
-----------------------

Tests for `processors.orthorectification.orthorectification` module.
"""

import os
from unittest import TestCase
from zipfile import ZipFile
import tempfile
import shutil
from copy import deepcopy

from pyproj import CRS
import numpy as np
from geoarray import GeoArray

from enpt.processors.orthorectification import Orthorectifier, VNIR_SWIR_Stacker
from enpt.options.config import config_for_testing, config_for_testing_dlr, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL2Product_MapGeo

__author__ = 'Daniel Scheffler'


class Test_Orthorectifier(TestCase):
    def setUp(self):
        self.config = EnPTConfig(**config_for_testing)

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
        OR = Orthorectifier(config=self.config)
        L2_obj = OR.run_transformation(self.L1_obj)

        self.assertIsInstance(L2_obj, EnMAPL2Product_MapGeo)
        self.assertTrue(L2_obj.data.is_map_geo)
        self.assertGreater(L2_obj.data.shape[0], self.L1_obj.vnir.data.shape[0])
        self.assertNotEqual(L2_obj.data.shape[1], self.L1_obj.vnir.data.shape[1])
        self.assertEqual(L2_obj.data.ndim, self.L1_obj.vnir.data.ndim)
        self.assertTrue(np.isclose(np.mean(self.L1_obj.vnir.data[:, :, 0]),
                                   np.mean(L2_obj.data[:, :, 0][L2_obj.data[:, :, 0] != L2_obj.data.nodata]),
                                   rtol=0.01
                                   ))


class Test_Orthorectifier_DLR(TestCase):
    def setUp(self):
        self.config = EnPTConfig(**config_for_testing_dlr)

        # create a temporary directory
        # NOTE: This must exist during the whole runtime of Test_Orthorectifier, otherwise
        #       Orthorectifier.run_transformation will fail to read some files.
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

        # get lons / lats
        with ZipFile(self.config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(self.tmpdir)
            self.L1_obj = L1B_Reader(config=self.config).read_inputdata(
                root_dir_main=self.tmpdir,
                compute_snr=False)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_run_transformation(self):
        OR = Orthorectifier(config=self.config)
        L2_obj = OR.run_transformation(self.L1_obj)

        self.assertIsInstance(L2_obj, EnMAPL2Product_MapGeo)
        self.assertTrue(L2_obj.data.is_map_geo)
        self.assertGreater(L2_obj.data.shape[0], self.L1_obj.vnir.data.shape[0])
        self.assertNotEqual(L2_obj.data.shape[1], self.L1_obj.vnir.data.shape[1])
        self.assertEqual(L2_obj.data.ndim, self.L1_obj.vnir.data.ndim)
        self.assertTrue(np.isclose(np.mean(self.L1_obj.vnir.data[:, :, 0]),
                                   np.mean(L2_obj.data[:, :, 0][L2_obj.data[:, :, 0] != L2_obj.data.nodata]),
                                   rtol=0.01
                                   ))


class Test_VNIR_SWIR_Stacker(TestCase):
    def setUp(self):
        self.vnir_gA = GeoArray(np.random.randint(0, 255, (10, 10, 10)),
                                geotransform=(331185.0, 30.0, -0.0, 5840115.0, -0.0, -30.0),
                                projection=CRS(32633).to_wkt())
        self.swir_gA = GeoArray(np.random.randint(0, 255, (10, 10, 20)),
                                geotransform=(331185.0, 30.0, -0.0, 5840115.0, -0.0, -30.0),
                                projection=CRS(32633).to_wkt())
        self.vnir_wvls = np.arange(900, 1000, 10)
        self.swir_wvls = np.arange(935, 1135, 10)

        self.VSSt = VNIR_SWIR_Stacker(vnir=self.vnir_gA, swir=self.swir_gA,
                                      vnir_wvls=self.vnir_wvls, swir_wvls=self.swir_wvls)

    def test_validate_input(self):
        # unequal geotransform
        swir_gA = deepcopy(self.swir_gA)
        swir_gA.gt = (331185.0, 10.0, -0.0, 5840115.0, -0.0, -10.0)
        with self.assertRaises(ValueError):
            VNIR_SWIR_Stacker(vnir=self.vnir_gA, swir=swir_gA,
                              vnir_wvls=self.vnir_wvls, swir_wvls=self.swir_wvls)

        # unequal projection
        swir_gA = deepcopy(self.swir_gA)
        swir_gA.prj = CRS(32632).to_wkt()
        with self.assertRaises(ValueError):
            VNIR_SWIR_Stacker(vnir=self.vnir_gA, swir=swir_gA,
                              vnir_wvls=self.vnir_wvls, swir_wvls=self.swir_wvls)

        # wrong length of provided wavelength
        with self.assertRaises(ValueError):
            VNIR_SWIR_Stacker(vnir=self.vnir_gA, swir=self.swir_gA,
                              vnir_wvls=np.array(list(self.vnir_wvls) + [1]), swir_wvls=self.swir_wvls)
        with self.assertRaises(ValueError):
            VNIR_SWIR_Stacker(vnir=self.vnir_gA, swir=self.swir_gA,
                              vnir_wvls=self.vnir_wvls, swir_wvls=np.array(list(self.swir_wvls) + [1]))

    def validate_output(self, gA_stacked: GeoArray):
        self.assertIsInstance(gA_stacked, GeoArray)
        self.assertEqual(gA_stacked.gt, self.vnir_gA.gt)
        self.assertEqual(gA_stacked.prj, self.vnir_gA.prj)
        self.assertEqual(gA_stacked.shape[:2], self.vnir_gA.shape[:2])
        self.assertTrue('wavelength' in gA_stacked.meta.band_meta and
                        gA_stacked.meta.band_meta['wavelength'])
        self.assertEqual(gA_stacked.bands, len(gA_stacked.meta.band_meta['wavelength']))

    def test_get_stack_order_by_wvl(self):
        gA_stacked = self.VSSt.compute_stack(algorithm='order_by_wvl')
        self.validate_output(gA_stacked)

    def test_get_stack_average(self):
        gA_stacked = self.VSSt.compute_stack(algorithm='average')
        self.validate_output(gA_stacked)

    def test_get_stack_vnir_only(self):
        gA_stacked = self.VSSt.compute_stack(algorithm='vnir_only')
        self.validate_output(gA_stacked)

    def test_get_stack_swir_only(self):
        gA_stacked = self.VSSt.compute_stack(algorithm='swir_only')
        self.validate_output(gA_stacked)

    def test_compute_stack(self):
        # wrong input algorithm
        with self.assertRaises(ValueError):
            self.VSSt.compute_stack(algorithm='mean')


if __name__ == '__main__':
    import pytest
    pytest.main()
