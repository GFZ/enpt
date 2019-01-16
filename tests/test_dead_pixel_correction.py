# -*- coding: utf-8 -*-

from unittest import TestCase

import numpy as np
from geoarray import GeoArray

from enpt.processors.dead_pixel_correction import \
    Dead_Pixel_Corrector, \
    interp_nodata_along_axis_2d, \
    interp_nodata_along_axis


class Test_Dead_Pixel_Corrector(TestCase):

    def setUp(self):
        """Set up the needed test data"""
        self.DPC = Dead_Pixel_Corrector(algorithm='spectral', interp='linear', filter_halfwidth=2)

        # create test data
        self.im = np.random.randint(1, 10, (50, 1000, 88), np.int16)  # VNIR test size

        # create 2D dead pixel map
        self.deadpixelmap_2D = np.zeros((self.im.shape[2], self.im.shape[1]), np.bool)
        for band, column in \
            [[0, 0],  # first band, first column
             [0, 2],  # first band, any column
             [1, 2],  # second band, same column
             [50, 4],  # 2 adjacent bands
             [51, 4],  # 2 adjacent bands
             [60, 20],  # single dead column
             [86, 50],  # second last band, same column
             [87, 50],  # last band, same column
             [87, 2]]:  # single dead column, last band
            self.im[:, column, band] = 0
            self.deadpixelmap_2D[band, column] = 1

        # create 3D dead pixel map
        self.deadpixelmap_3D = self.im == 0

    def test_correct_using_2D_deadpixelmap(self):
        corrected = self.DPC.correct(self.im, self.deadpixelmap_2D)

        # output assertions
        self.assertIsInstance(corrected, (GeoArray, np.ndarray))
        # self.assertEqual(np.mean(self.im[:, 0, 0]), 0)  # first band, first column (currently not corrected)
        self.assertNotEqual(np.mean(self.im[:, 2, 0]), 0)  # first band, any column
        self.assertNotEqual(np.mean(self.im[:, 2, 1]), 0)  # second band, same column
        self.assertNotEqual(np.mean(self.im[:, 4, 50]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 4, 10]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 20, 60]), 0)  # single dead column
        self.assertNotEqual(np.mean(self.im[:, 50, 86]), 0)  # second last band, same column
        self.assertNotEqual(np.mean(self.im[:, 50, 87]), 0)  # last band, same column
        self.assertNotEqual(np.mean(self.im[:, 2, 87]), 0)  # single dead column, last band

    def test_correct_using_3D_deadpixelmap(self):
        corrected = self.DPC.correct(self.im, self.deadpixelmap_3D)

        # output assertions
        self.assertIsInstance(corrected, (GeoArray, np.ndarray))
        self.assertEqual(self.im.shape, corrected.shape)
        self.assertNotEqual(np.mean(self.im[:, 2, 0]), 0)  # first band, any column
        self.assertNotEqual(np.mean(self.im[:, 2, 1]), 0)  # first band, same column
        self.assertNotEqual(np.mean(self.im[:, 4, 50]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 4, 10]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 20, 60]), 0)  # single dead column
        self.assertNotEqual(np.mean(self.im[:, 50, 86]), 0)  # second last band, same column
        self.assertNotEqual(np.mean(self.im[:, 50, 87]), 0)  # last band, same column
        self.assertNotEqual(np.mean(self.im[:, 2, 87]), 0)  # single dead column, last band


class Test_interp_nodata_along_axis_2d(TestCase):
    @staticmethod
    def get_data2d():
        return np.array([[0, 0, 2],
                         [3, np.nan, 5],
                         [np.nan, 10, 8]])

    def test_axis_0(self):
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, method='linear')
        arr_exp = np.array([[0, 0, 2], [3, 5, 5], [6, 10, 8]])
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        mask_nodata = ~np.isfinite(self.get_data2d())
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, nodata=mask_nodata, method='linear')
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, method='linear', fill_value=-1)
        arr_exp = np.array([[0, 0, 2], [3, 5, 5], [-1, 10, 8]])
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

    def test_axis_1(self):
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, method='linear')
        arr_exp = np.array([[0, 0, 2], [3, 4, 5], [12, 10, 8]])
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        mask_nodata = ~np.isfinite(self.get_data2d())
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, nodata=mask_nodata, method='linear')
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, method='linear', fill_value=-1)
        arr_exp = np.array([[0, 0, 2], [3, 4, 5], [-1, 10, 8]])
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

    def test_bad_args(self):
        with self.assertRaises(ValueError):
            interp_nodata_along_axis_2d(self.get_data2d(), axis=3)
        with self.assertRaises(ValueError):
            interp_nodata_along_axis_2d(np.dstack([self.get_data2d(), self.get_data2d()]))


class Test_interp_nodata_along_axis(TestCase):
    @staticmethod
    def get_data3d():
        data3d = np.zeros((3, 3, 3))
        data3d[:, :, 0] = [[0, 0, 2],
                           [3, np.nan, 5],
                           [np.nan, 10, 8]]
        data3d[:, :, 1] = [[10, 10, 12],
                           [13, np.nan, 15],
                           [16, 110, np.nan]]
        data3d[:, :, 2] = [[20, 20, 22],
                           [23, np.nan, 25],
                           [np.nan, 210, 20]]

        return data3d

    def test_3d_axis_0(self):
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=0, method='linear')
        arr_exp = np.zeros((3, 3, 3))
        arr_exp[:, :, 0] = [[0, 0, 2], [3, 5, 5], [6, 10, 8]]
        arr_exp[:, :, 1] = [[10, 10, 12], [13, 60, 15], [16, 110, 18]]
        arr_exp[:, :, 2] = [[20, 20, 22], [23, 115, 25], [26, 210, 20]]
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        mask_nodata = ~np.isfinite(self.get_data3d())
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=0, nodata=mask_nodata, method='linear')
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

    def test_3d_axis_1(self):
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=1, method='linear')
        arr_exp = np.zeros((3, 3, 3))
        arr_exp[:, :, 0] = [[0, 0, 2], [3, 4, 5], [12, 10, 8]]
        arr_exp[:, :, 1] = [[10, 10, 12], [13, 14, 15], [16, 110, 204]]
        arr_exp[:, :, 2] = [[20, 20, 22], [23, 24, 25], [400, 210, 20]]
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

        mask_nodata = ~np.isfinite(self.get_data3d())
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=1, nodata=mask_nodata, method='linear')
        self.assertTrue(np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int)

    def test_3d_axis_2(self):
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=2, method='linear')
        arr_exp = np.zeros((3, 3, 3))
        arr_exp[:, :, 0] = [[0, 0, 2], [3, np.nan, 5], [np.nan, 10, 8]]
        arr_exp[:, :, 1] = [[10, 10, 12], [13, np.nan, 15], [16, 110, 14]]
        arr_exp[:, :, 2] = [[20, 20, 22], [23, np.nan, 25], [np.nan, 210, 20]]
        np.testing.assert_array_equal(data_int, arr_exp, 'Computed %s.' % data_int)

        mask_nodata = ~np.isfinite(self.get_data3d())
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=2, nodata=mask_nodata, method='linear')
        np.testing.assert_array_equal(data_int, arr_exp, 'Computed %s.' % data_int)

    def test_2d(self):
        data_int = interp_nodata_along_axis(Test_interp_nodata_along_axis_2d.get_data2d())
        self.assertTrue(np.array_equal(data_int, np.array([[0, 0, 2],
                                                           [3, 5, 5],
                                                           [6, 10, 8]])), 'Computed %s.' % data_int)

    def test_bad_args(self):
        with self.assertRaises(ValueError):
            interp_nodata_along_axis(np.array([1, 2, 3]))
