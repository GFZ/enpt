#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_spatial_transform
----------------------

Tests for `processors.spatial_transform.spatial_transform` module.
"""

import os
from typing import Tuple  # noqa: F401
from unittest import TestCase
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import numpy as np

from geoarray import GeoArray
from py_tools_ds.geo.coord_grid import is_point_on_grid

from enpt.processors.spatial_transform import Geometry_Transformer, RPC_Geolayer_Generator
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
        # self.rpc_coeffs = dict(
        #     row_off=512.0, col_off=500.0, lat_off=47.6257991503, lon_off=10.9372682697,
        #     height_off=1785.9975032126, row_scale=517.1199951172, col_scale=504.9999952316,
        #     lat_scale=0.1649981434, lon_scale=0.2274676529, height_scale=1499.847892108,
        #     row_num_coeffs=np.array(
        #         [3.57047913e-03, -2.74808704e-01, -1.16541133e+00, -2.86627308e-04, -2.68319570e-01,
        #          3.48077310e-04, 1.64050321e-03, -7.46389036e-02, -3.24470589e-01, 1.05686138e-04,
        #          5.15176091e-03, 9.04990036e-03, 4.59846719e-02, 7.66279483e-03, 5.89784062e-02,
        #          8.07669634e-02, 3.24461847e-02, 6.60198890e-04, 9.61177339e-03, 7.03402953e-06]),
        #     row_den_coeffs=np.array(
        #         [1.00000000e+00, 1.75171420e-01, 2.74865227e-01, -1.45692178e-03, -2.09325421e-02,
        #          -2.43456641e-03, -8.18513259e-03, -5.04700287e-02, -7.07800273e-02, -2.78124516e-02,
        #          -1.31225254e-05, 3.02270705e-04, -1.85584425e-03, -3.16658113e-06, 3.20667033e-04,
        #          -1.76663639e-03, 8.62457306e-06, 1.58896500e-04, 1.06977978e-04, 3.45796280e-05]),
        #     col_num_coeffs=np.array([
        #         -1.95581836e-02, 1.14502576e+00, -2.76461594e-01, 1.70538121e-04, 2.64971210e-01,
        #         1.36009430e-03, -2.16623832e-04, 2.01329933e-01, -7.28906110e-02, 5.49497388e-04,
        #         -8.08011145e-03, -5.83776374e-02, -7.47574926e-02, -3.19433457e-02, -1.06122613e-02,
        #         2.17612444e-02, 7.50293073e-03, -2.28926437e-03, 1.96028015e-03, -4.49157359e-06]),
        #     col_den_coeffs=np.array([
        #         1.00000000e+00, 1.75171420e-01, 2.74865227e-01, -1.45692178e-03, -2.09325421e-02,
        #         -2.43456641e-03, -8.18513259e-03, -5.04700287e-02, -7.07800273e-02, -2.78124516e-02,
        #         -1.31225254e-05, 3.02270705e-04, -1.85584425e-03, -3.16658113e-06, 3.20667033e-04,
        #         -1.76663639e-03, 8.62457306e-06, 1.58896500e-04, 1.06977978e-04, 3.45796280e-05])
        # )
        import pickle
        with open('/home/gfz-fe/scheffler/temp/enpt_testing/dlr_l1b_test_data/rpc_coeffs_B200.pkl', 'rb') as dillF:
            self.rpc_coeffs = pickle.load(dillF)

        # bounding polygon DLR test data
        self.lats = np.array([47.7872236, 47.7232358, 47.5195676, 47.4557831])
        self.lons = np.array([10.7966311, 11.1693436, 10.7111131, 11.0815993])
        corner_coords = tuple(zip(self.lons, self.lats))  # type: Tuple[Tuple[float, float]]

        # spatial coverage of datatake DLR test data
        # self.lats = np.array([47.7870358956, 47.723060779, 46.9808418244, 46.9174014681])
        # self.lons = np.array([10.7968099213, 11.1693752478, 10.5262233116, 10.8932492494])

        self.heights = np.array([764, 681, 842, 1539])  # values from ASTER DEM
        # TODO validate dem before passing to RPC_Geolayer_Generator
        self.dem = '/home/gfz-fe/scheffler/temp/enpt_testing/dlr_l1b_test_data/DLR_L2A_DEM_UTM32.bsq'
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


