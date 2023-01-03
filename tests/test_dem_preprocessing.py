#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2023 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
test_dem_preprocessing
----------------------

Tests for `processors.dem_preprocessing` module.
"""

from unittest import TestCase

import numpy as np
from pyproj import CRS
from geoarray import GeoArray

from enpt.processors.dem_preprocessing import DEM_Processor
from enpt.options.config import config_for_testing
from enpt.model.metadata import EnMAP_Metadata_L1B_Detector_SensorGeo

__author__ = 'Daniel Scheffler'


class Test_DEM_Processor(TestCase):
    def setUp(self):
        self.path_demfile = config_for_testing['path_dem']

        # get lons/lats
        lat_UL_UR_LL_LR = [47.545359963328366, 47.48153190433143, 47.505282507365813, 47.441546248160961]
        lon_UL_UR_LL_LR = [10.701359191637021, 11.072698329235017, 10.686064194247395, 11.056744608586392]
        self.ll_cornerCoords = tuple(zip(lon_UL_UR_LL_LR, lat_UL_UR_LL_LR))
        self.lats = EnMAP_Metadata_L1B_Detector_SensorGeo.interpolate_corners(*lat_UL_UR_LL_LR, 1000, 100)
        self.lons = EnMAP_Metadata_L1B_Detector_SensorGeo.interpolate_corners(*lon_UL_UR_LL_LR, 1000, 100)

        self.DP_mapgeo = DEM_Processor(self.path_demfile, enmapIm_cornerCoords=self.ll_cornerCoords)

    def test_init_nomapinfo(self):
        dem = GeoArray(np.array([1, 2]))

        # no map info, no projection
        with self.assertRaises(ValueError):
            DEM_Processor(dem, enmapIm_cornerCoords=self.ll_cornerCoords)

        # no projection
        dem.gt = (10.6, 0.00036, -0.0, 47.5, -0.0, -0.00036)  # can be anything
        with self.assertRaises(ValueError):
            DEM_Processor(dem, enmapIm_cornerCoords=self.ll_cornerCoords)

    def test_init_noWGS84(self):
        # NAD83 projection ({'proj': 'longlat', 'ellps': 'GRS80', 'towgs84': '0,0,0,0,0,0,0'})
        with self.assertRaises(ValueError):
            dem = GeoArray(np.array([1, 2]),
                           geotransform=(10.6, 0.00036, -0.0, 47.5, -0.0, -0.00036),  # can be anything
                           projection=CRS(4269).to_wkt())  # NAD83
            DEM_Processor(dem, enmapIm_cornerCoords=self.ll_cornerCoords)

    def test_init_demTooSmall(self):
        # covers only 100x100 px in the upper left (<5%)
        dem = GeoArray(np.random.randint(0, 500, (100, 100)),
                       geotransform=(626938.928052, 30.0, 0, 5267203.56579, 0, -30.0),
                       projection=CRS(32632).to_wkt())

        with self.assertRaises(ValueError):
            DEM_Processor(dem, enmapIm_cornerCoords=self.ll_cornerCoords)

    def test_get_flat_dem(self):
        DP = DEM_Processor.get_flat_dem_from_average_elevation(corner_coords_lonlat=self.ll_cornerCoords,
                                                               average_elevation=50)
        self.assertIsInstance(DP.dem, GeoArray)
        self.assertEqual(np.mean(DP.dem), 50)

    def test_fill_gaps(self):
        pass

    def test_compute_slopes(self):
        pass

    def test_compute_aspect(self):
        pass

    def test_to_sensor_geometry(self):
        dem_sensor_geo = self.DP_mapgeo.to_sensor_geometry(lons=self.lons, lats=self.lats)

        self.assertEqual(dem_sensor_geo.shape, (100, 1000))


if __name__ == '__main__':
    import pytest
    pytest.main()
