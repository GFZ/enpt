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
import numpy as np

from geoarray import GeoArray
from py_tools_ds.geo.coord_grid import is_point_on_grid

from enpt.processors.spatial_transform.spatial_transform import Geometry_Transformer
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.options.config import enmap_coordinate_grid


class Test_Geometry_Transformer(TestCase):
    def setUp(self):
        config = EnPTConfig(**config_for_testing)

        # get lons / lats
        with TemporaryDirectory() as td, ZipFile(config.path_l1b_enmap_image, "r") as zf:
            zf.extractall(td)
            L1_obj = L1B_Reader(config=config).read_inputdata(
                root_dir_main=os.path.join(td, os.path.splitext(os.path.basename(config.path_l1b_enmap_image))[0]),
                compute_snr=False)

        R, C = L1_obj.vnir.data.shape[:2]
        self.lons_vnir = L1_obj.meta.vnir.interpolate_corners(*L1_obj.meta.vnir.lon_UL_UR_LL_LR, nx=C, ny=R)
        self.lats_vnir = L1_obj.meta.vnir.interpolate_corners(*L1_obj.meta.vnir.lat_UL_UR_LL_LR, nx=C, ny=R)

        self.gA2transform_sensorgeo = L1_obj.vnir.data[:, :, 50]  # a single VNIR band in sensor geometry
        self.gA2transform_mapgeo = GeoArray(config.path_dem)  # a DEM in map geometry given by the user

    def test_to_map_geometry(self):
        # transforming map to map geometry must raise a RuntimeError
        with self.assertRaises(RuntimeError):
            GT = Geometry_Transformer(self.gA2transform_mapgeo, lons=self.lons_vnir, lats=self.lats_vnir)
            GT.to_map_geometry(tgt_prj=32632)

        # test transformation to UTM zone 32
        GT = Geometry_Transformer(self.gA2transform_sensorgeo, lons=self.lons_vnir, lats=self.lats_vnir)
        data_mapgeo, gt, prj = GT.to_map_geometry(tgt_prj=32632)
        self.assertEqual((gt[1], -gt[5]), (np.ptp(enmap_coordinate_grid['x']),
                                           np.ptp(enmap_coordinate_grid['x'])))  # 30m output
        self.assertTrue(is_point_on_grid((gt[0], gt[3]),
                                         xgrid=enmap_coordinate_grid['x'], ygrid=enmap_coordinate_grid['y']))

        # test transformation to LonLat
        GT = Geometry_Transformer(self.gA2transform_sensorgeo, lons=self.lons_vnir, lats=self.lats_vnir)
        GT.to_map_geometry(tgt_prj=4326)

    def test_to_sensor_geometry(self):
        # transforming sensor to sensor geometry must raise a RuntimeError
        with self.assertRaises(RuntimeError):
            GT = Geometry_Transformer(self.gA2transform_sensorgeo, lons=self.lons_vnir, lats=self.lats_vnir)
            GT.to_sensor_geometry()

        GT = Geometry_Transformer(self.gA2transform_mapgeo, lons=self.lons_vnir, lats=self.lats_vnir)
        data_sensorgeo = GT.to_sensor_geometry()

        self.assertEqual(data_sensorgeo.shape, self.gA2transform_sensorgeo.shape)
