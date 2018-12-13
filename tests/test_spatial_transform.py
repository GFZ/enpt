#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_spatial_transform
----------------------

Tests for `processors.spatial_transform.spatial_transform` module.
"""

import os
from unittest import TestCase
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from geoarray import GeoArray


from enpt.processors.spatial_transform.spatial_transform import Geometry_Transformer
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader


class Test_Geometry_Transformer(TestCase):
    def setUp(self):
        config = EnPTConfig(**config_for_testing)

        self.gA2transform = GeoArray(config.path_dem)

        # get lons / lats
        with TemporaryDirectory() as td, ZipFile(config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(td)
            L1_obj = L1B_Reader(config=config).read_inputdata(
                root_dir_main=os.path.join(td, os.path.splitext(os.path.basename(config.path_l1b_enmap_image))[0]),
                compute_snr=False)

        R, C = L1_obj.vnir.data.shape[:2]
        self.lons_vnir = L1_obj.meta.vnir.interpolate_corners(*L1_obj.meta.vnir.lon_UL_UR_LL_LR, nx=C, ny=R)
        self.lats_vnir = L1_obj.meta.vnir.interpolate_corners(*L1_obj.meta.vnir.lat_UL_UR_LL_LR, nx=C, ny=R)

    def test_to_map_geometry(self):
        GT = Geometry_Transformer(self.gA2transform, lons=self.lons_vnir, lats=self.lats_vnir)
        GT.to_map_geometry(tgt_prj=32632)
