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
test_images_sensorgeo
---------------------

Tests for `model.images.images_sensorgeo` module.
"""

from unittest import TestCase
import pytest
from zipfile import ZipFile
import tempfile
import shutil
import numpy as np
from geoarray import GeoArray

from enpt.options.config import config_for_testing_dlr, EnPTConfig
from enpt.io.reader import L1B_Reader


__author__ = 'Daniel Scheffler'


class Test_EnMAPL1Product_SensorGeo(TestCase):
    config = None
    tmpdir = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.config = EnPTConfig(**config_for_testing_dlr)

        cls.tmpdir = tempfile.mkdtemp(dir=cls.config.working_dir)

        # get an EnMAPL1Product_SensorGeo instance
        with ZipFile(cls.config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(cls.tmpdir)
            cls.L1_obj = L1B_Reader(config=cls.config).read_inputdata(
                root_dir_main=cls.tmpdir,
                compute_snr=False)

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.tmpdir)  # noqa

    def test_set_SWIRattr_with_transformedVNIRattr(self):
        assert self.L1_obj.swir.mask_snow is None

        self.L1_obj.set_SWIRattr_with_transformedVNIRattr('mask_snow')

        assert isinstance(self.L1_obj.swir.mask_snow, GeoArray)

    def test_set_SWIRattr_with_transformedVNIRattr__attrIsNone(self):
        del self.L1_obj.vnir.mask_haze

        with pytest.raises(RuntimeError, match='.vnir.mask_haze has not yet been set.'):
            self.L1_obj.set_SWIRattr_with_transformedVNIRattr('mask_haze')

    def test_transform_vnir_to_swir_raster_with_keystone(self):
        vnirdata_swirgeo = self.L1_obj.transform_vnir_to_swir_raster(self.L1_obj.vnir.data[:],
                                                                     respect_keystone=True)

        assert isinstance(vnirdata_swirgeo, np.ndarray)
        assert vnirdata_swirgeo.shape == self.L1_obj.vnir.data.shape

        # at least the last 10 lines must be zero due to the VNIR/SWIR shift
        assert np.mean(vnirdata_swirgeo[-10:, :, :]) == 0

    def test_transform_swir_to_vnir_raster_no_keystone(self):
        swirdata_vnirgeo = self.L1_obj.transform_swir_to_vnir_raster(self.L1_obj.swir.data[:],
                                                                     respect_keystone=False)

        assert isinstance(swirdata_vnirgeo, np.ndarray)
        assert swirdata_vnirgeo.shape == self.L1_obj.swir.data.shape

        # at least the first 10 lines must be zero due to the VNIR/SWIR shift
        assert np.mean(swirdata_vnirgeo[:10, :, :]) == 0

    def test_transform_vnir_to_swir_raster__no_geolayer(self):
        vnir_lons = self.L1_obj.meta.vnir.lons

        try:
            self.L1_obj.meta.vnir.lons = None

            with pytest.raises(RuntimeError, match='The VNIR/SWIR geolayers must be computed first '
                                                   'to transform arrays from VNIR to SWIR sensor geometry.'):
                self.L1_obj.transform_vnir_to_swir_raster(self.L1_obj.vnir.data[:])

        finally:
            self.L1_obj.meta.vnir.lons = vnir_lons
            assert isinstance(self.L1_obj.meta.vnir.lons, np.ndarray)

    def test_transform_swir_to_vnir_raster__no_geolayer(self):
        vnir_lons = self.L1_obj.meta.vnir.lons

        try:
            self.L1_obj.meta.vnir.lons = None

            with pytest.raises(RuntimeError, match='The VNIR/SWIR geolayers must be computed first '
                                                   'to transform arrays from VNIR to SWIR sensor geometry.'):
                self.L1_obj.transform_swir_to_vnir_raster(self.L1_obj.swir.data[:])

        finally:
            self.L1_obj.meta.vnir.lons = vnir_lons
            assert isinstance(self.L1_obj.meta.vnir.lons, np.ndarray)


if __name__ == '__main__':
    pytest.main()
