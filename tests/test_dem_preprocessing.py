#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_dem_preprocessing
----------------------

Tests for `processors.dem_preprocessing` module.
"""

from unittest import TestCase

from py_tools_ds.geo.projection import prj_equal
from geoarray import GeoArray

from enpt.processors.dem_preprocessing import DEM_Processor
from enpt.options.config import config_for_testing
from enpt.model.metadata import EnMAP_Metadata_L1B_Detector_SensorGeo


class Test_DEM_Processor(TestCase):
    def setUp(self):
        self.path_demfile = config_for_testing['path_dem']
        self.DP_mapgeo = DEM_Processor(self.path_demfile)

        # get lons/lats
        lat_UL_UR_LL_LR = [47.545359963328366, 47.48153190433143, 47.505282507365813, 47.441546248160961]
        lon_UL_UR_LL_LR = [10.701359191637021, 11.072698329235017, 10.686064194247395, 11.056744608586392]
        self.lats = EnMAP_Metadata_L1B_Detector_SensorGeo.interpolate_corners(*lat_UL_UR_LL_LR, 1000, 100)
        self.lons = EnMAP_Metadata_L1B_Detector_SensorGeo.interpolate_corners(*lon_UL_UR_LL_LR, 1000, 100)

    def test_fill_gaps(self):
        pass

    def test_compute_slopes(self):
        pass

    def test_compute_aspect(self):
        pass

    def test_to_sensor_geometry(self):
        dem_sensor_geo = self.DP_mapgeo.to_sensor_geometry(lons=self.lons, lats=self.lats)

        self.assertEquals(dem_sensor_geo.shape, (100, 1000))

    def test_to_map_geometry(self):
        dem_sensor_geo = self.DP_mapgeo.to_sensor_geometry(lons=self.lons, lats=self.lats)

        DP_sensorgeo = DEM_Processor(GeoArray(dem_sensor_geo))
        dem_map_geo_gA = DP_sensorgeo.to_map_geometry(lons=self.lons, lats=self.lats, tgt_prj=32632)  # UTM32

        self.assertTrue(prj_equal(dem_map_geo_gA.prj, 32632))
