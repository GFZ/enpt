#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_spatial_transform
----------------------

Tests for `processors.spatial_transform.spatial_transform` module.
"""

import os
from typing import Tuple  # noqa: F401
from unittest import TestCase, SkipTest
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import pickle
import numpy as np

from geoarray import GeoArray
from py_tools_ds.geo.coord_grid import is_point_on_grid

from enpt.processors.spatial_transform import Geometry_Transformer, RPC_Geolayer_Generator
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader
from enpt.options.config import enmap_coordinate_grid


path_testdata = os.path.abspath(os.path.join(__file__, '..', 'data'))


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
        GT = Geometry_Transformer(lons=self.lons_vnir, lats=self.lats_vnir)

        # transforming map to map geometry must raise a RuntimeError
        with self.assertRaises(RuntimeError):
            GT.to_map_geometry(self.gA2transform_mapgeo, tgt_prj=32632)

        # test transformation to UTM zone 32
        data_mapgeo, gt, prj = GT.to_map_geometry(self.gA2transform_sensorgeo, tgt_prj=32632)
        self.assertEqual((gt[1], -gt[5]), (np.ptp(enmap_coordinate_grid['x']),
                                           np.ptp(enmap_coordinate_grid['x'])))  # 30m output
        self.assertTrue(is_point_on_grid((gt[0], gt[3]),
                                         xgrid=enmap_coordinate_grid['x'], ygrid=enmap_coordinate_grid['y']))

        # test transformation to LonLat
        GT.to_map_geometry(self.gA2transform_sensorgeo, tgt_prj=4326)

    def test_to_sensor_geometry(self):
        GT = Geometry_Transformer(lons=self.lons_vnir, lats=self.lats_vnir)

        # transforming sensor to sensor geometry must raise a RuntimeError
        with self.assertRaises(RuntimeError):
            GT.to_sensor_geometry(self.gA2transform_sensorgeo)

        data_sensorgeo = GT.to_sensor_geometry(self.gA2transform_mapgeo)

        self.assertEqual(data_sensorgeo.shape, self.gA2transform_sensorgeo.shape)


class Test_RPC_Geolayer_Generator(TestCase):
    def setUp(self):
        with open(os.path.join(path_testdata, 'rpc_coeffs_B200.pkl'), 'rb') as pklF:
            self.rpc_coeffs = pickle.load(pklF)

        # bounding polygon DLR test data
        self.lats = np.array([47.7872236, 47.7232358, 47.5195676, 47.4557831])
        self.lons = np.array([10.7966311, 11.1693436, 10.7111131, 11.0815993])
        corner_coords = tuple(zip(self.lons, self.lats))  # type: Tuple[Tuple[float, float]]

        # spatial coverage of datatake DLR test data
        # self.lats = np.array([47.7870358956, 47.723060779, 46.9808418244, 46.9174014681])
        # self.lons = np.array([10.7968099213, 11.1693752478, 10.5262233116, 10.8932492494])

        self.heights = np.array([764, 681, 842, 1539])  # values from ASTER DEM
        # TODO validate dem before passing to RPC_Geolayer_Generator
        self.dem = os.path.join(path_testdata, 'DLR_L2A_DEM_UTM32.bsq')
        self.dims_sensorgeo = (1024, 1000)

        self.RPCGG = RPC_Geolayer_Generator(self.rpc_coeffs,
                                            dem=self.dem,
                                            enmapIm_cornerCoords=corner_coords,
                                            enmapIm_dims_sensorgeo=self.dims_sensorgeo)

    def test_normalize_coordinates(self):
        lon_norm, lat_norm, height_norm = \
            self.RPCGG._normalize_map_coordinates(lon=self.lons, lat=self.lats, height=self.heights)
        self.assertEquals(lon_norm.shape, self.lons.shape)
        self.assertEquals(lat_norm.shape, self.lats.shape)
        self.assertEquals(height_norm.shape, self.heights.shape)

    def test_compute_normalized_image_coordinates(self):
        row_norm, col_norm = self.RPCGG._compute_normalized_image_coordinates(
            lon_norm=np.array([-0.61827327,  1.02025641, -0.99423002,  0.63451233]),
            lat_norm=np.array([0.97834101,  0.59053179, -0.64383482, -1.0304119]),
            height_norm=np.array([-0.85741862, -0.79074519, -0.72407176, -0.65739833]),
        )

        rows, cols = self.RPCGG._denormalize_image_coordinates(row_norm, col_norm)

    def test_transform_LonLatHeight_to_RowCol(self):
        rows, cols = self.RPCGG.transform_LonLatHeight_to_RowCol(lon=self.lons, lat=self.lats, height=self.heights)

    def test_compute_geolayer(self):
        lons_interp, lats_interp = self.RPCGG.compute_geolayer()
        self.assertEquals(lons_interp.shape, lats_interp.shape)
        self.assertEquals(lons_interp.shape, self.dims_sensorgeo)
        self.assertFalse(np.isnan(lons_interp).any())
        self.assertFalse(np.isnan(lats_interp).any())


