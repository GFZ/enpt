#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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

from unittest import TestCase
import pytest

import numpy as np
from geoarray import GeoArray

from enpt.processors.dead_pixel_correction import \
    Dead_Pixel_Corrector, \
    interp_nodata_along_axis_2d, \
    interp_nodata_along_axis, \
    interp_nodata_spatially_2d

__author__ = 'Daniel Scheffler'


class Test_Dead_Pixel_Corrector(TestCase):

    def setUp(self):
        """Set up the needed test data"""
        # create test data
        self.im = np.random.randint(1, 10, (50, 1000, 88)).astype(float)  # VNIR test size

        # create 2D dead pixel map
        self.deadpixelmap_2D = np.zeros((self.im.shape[2], self.im.shape[1]), bool)
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
            self.im[:, column, band] = np.nan
            self.deadpixelmap_2D[band, column] = 1

        # create 3D dead pixel map
        self.deadpixelmap_3D = np.isnan(self.im)

    def validate_output_spectral_interp(self, output):
        assert isinstance(output, (GeoArray, np.ndarray))
        assert self.im.shape == output.shape
        assert np.mean(output[:, 0, 0]) != 0  # first band, first column
        assert np.mean(output[:, 2, 0]) != 0  # first band, any column
        assert np.mean(output[:, 2, 1]) != 0  # second band, same column
        assert np.mean(output[:, 4, 50]) != 0  # 2 adjacent bands
        assert np.mean(output[:, 4, 10]) != 0  # 2 adjacent bands
        assert np.mean(output[:, 20, 60]) != 0  # single dead column
        assert np.mean(output[:, 50, 86]) != 0  # second last band, same column
        assert np.mean(output[:, 50, 87]) != 0  # last band, same column
        assert np.mean(output[:, 2, 87]) != 0  # single dead column, last band

    def test_correct_using_2D_deadpixelmap_spectral(self):
        DPC = Dead_Pixel_Corrector(algorithm='spectral', interp_spectral='linear')
        corrected = DPC.correct(self.im, self.deadpixelmap_2D)

        # output assertions
        self.validate_output_spectral_interp(corrected)

    def test_correct_using_3D_deadpixelmap_spectral(self):
        DPC = Dead_Pixel_Corrector(algorithm='spectral', interp_spectral='linear')
        corrected = DPC.correct(self.im, self.deadpixelmap_3D)

        # output assertions
        self.validate_output_spectral_interp(corrected)

    def test_correct_using_2D_deadpixelmap_spatial(self):
        DPC = Dead_Pixel_Corrector(algorithm='spatial', interp_spatial='linear')
        corrected = DPC.correct(self.im, self.deadpixelmap_2D)

        # output assertions
        self.validate_output_spectral_interp(corrected)

    def test_correct_using_3D_deadpixelmap_spatial(self):
        DPC = Dead_Pixel_Corrector(algorithm='spatial', interp_spatial='linear')
        corrected = DPC.correct(self.im, self.deadpixelmap_3D)

        # output assertions
        self.validate_output_spectral_interp(corrected)


class Test_interp_nodata_along_axis_2d(TestCase):
    @staticmethod
    def get_data2d():
        return np.array([[0, 0, 2],
                         [3, np.nan, 5],
                         [np.nan, 10, 8]])

    def test_axis_0(self):
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, method='linear')
        arr_exp = np.array([[0, 0, 2], [3, 5, 5], [6, 10, 8]])
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        mask_nodata = ~np.isfinite(self.get_data2d())
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, nodata=mask_nodata, method='linear')
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=0, method='linear', fill_value=-1)
        arr_exp = np.array([[0, 0, 2], [3, 5, 5], [-1, 10, 8]])
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

    def test_axis_1(self):
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, method='linear')
        arr_exp = np.array([[0, 0, 2], [3, 4, 5], [12, 10, 8]])
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        mask_nodata = ~np.isfinite(self.get_data2d())
        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, nodata=mask_nodata, method='linear')
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        data_int = interp_nodata_along_axis_2d(self.get_data2d(), axis=1, method='linear', fill_value=-1)
        arr_exp = np.array([[0, 0, 2], [3, 4, 5], [-1, 10, 8]])
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

    def test_bad_args(self):
        with pytest.raises(ValueError):
            interp_nodata_along_axis_2d(self.get_data2d(), axis=3)
        with pytest.raises(ValueError):
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
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        mask_nodata = ~np.isfinite(self.get_data3d())
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=0, nodata=mask_nodata, method='linear')
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

    def test_3d_axis_1(self):
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=1, method='linear')
        arr_exp = np.zeros((3, 3, 3))
        arr_exp[:, :, 0] = [[0, 0, 2], [3, 4, 5], [12, 10, 8]]
        arr_exp[:, :, 1] = [[10, 10, 12], [13, 14, 15], [16, 110, 204]]
        arr_exp[:, :, 2] = [[20, 20, 22], [23, 24, 25], [400, 210, 20]]
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

        mask_nodata = ~np.isfinite(self.get_data3d())
        data_int = interp_nodata_along_axis(self.get_data3d(), axis=1, nodata=mask_nodata, method='linear')
        assert np.array_equal(data_int, arr_exp), 'Computed %s.' % data_int

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
        assert np.array_equal(data_int, np.array([[0, 0, 2],
                                                  [3, 5, 5],
                                                  [6, 10, 8]])), \
            'Computed %s.' % data_int

    def test_bad_args(self):
        with pytest.raises(ValueError):
            interp_nodata_along_axis(np.array([1, 2, 3]))


class Test_interp_nodata_spatially_2d(TestCase):
    @staticmethod
    def get_data2d():
        return np.array([[0, 0, 2, 12],
                         [3, np.nan, 5, np.nan],
                         [np.nan, 20, 8, 3]])

    def test_interpolation_scipy(self):
        data_int_scipy = interp_nodata_spatially_2d(self.get_data2d(), nodata=np.nan, method='linear',
                                                    fill_value=np.nan, implementation='scipy')
        arr_exp_scipy = np.array([[0, 0, 2, 12], [3, 10, 5, 7.5], [np.nan, 20, 8, 3]])
        np.testing.assert_array_equal(data_int_scipy, arr_exp_scipy, 'Computed %s.' % data_int_scipy)

        mask_nodata = ~np.isfinite(self.get_data2d())
        data_int_scipy = interp_nodata_spatially_2d(self.get_data2d(), nodata=mask_nodata, method='linear',
                                                    fill_value=np.nan, implementation='scipy')
        np.testing.assert_array_equal(data_int_scipy, arr_exp_scipy, 'Computed %s.' % data_int_scipy)

    def test_interpolation_pandas(self):
        data_int_pandas = interp_nodata_spatially_2d(self.get_data2d(), nodata=np.nan, method='linear',
                                                     fill_value=np.nan, implementation='pandas')
        arr_exp_pandas = np.array([[0, 0, 2, 12], [3, 10, 5, 7.5], [3, 20, 8, 3]])
        np.testing.assert_array_equal(data_int_pandas, arr_exp_pandas, 'Computed %s.' % data_int_pandas)

    def test_bad_args(self):
        with pytest.raises(ValueError):
            interp_nodata_spatially_2d(np.array([1, 2, 3]))
        with pytest.raises(ValueError):
            interp_nodata_spatially_2d(self.get_data2d(), nodata=np.array([1, 2, 3]))
        with pytest.raises(ValueError):
            interp_nodata_spatially_2d(self.get_data2d(), implementation='invalid')


if __name__ == '__main__':
    pytest.main()
