# -*- coding: utf-8 -*-
"""EnPT pre-processing module for digital elevation models."""

from typing import Union  # noqa: F401
from multiprocessing import cpu_count
import numpy as np

from geoarray import GeoArray

from ..spatial_transform import Geometry_Transformer


class DEM_Processor(object):
    def __init__(self, dem_path_geoarray: Union[str, GeoArray], CPUs: int = None):
        self.dem = GeoArray(dem_path_geoarray)
        self.CPUs = CPUs or cpu_count()

        self._validate_input()

    def _validate_input(self):
        # TODO - check geographic datum
        #      - check overlap
        pass

    def fill_gaps(self):
        pass

    def compute_slopes(self):
        # compute on map geometry (as provided)
        pass

    def compute_aspect(self):
        # compute on map geometry (as provided)
        pass

    def to_map_geometry(self,
                        lons: np.ndarray,
                        lats: np.ndarray,
                        tgt_prj: Union[str, int] = None):
        GT = Geometry_Transformer(self.dem, lons=lons, lats=lats, nprocs=self.CPUs)
        data_mapgeo, gt, prj = GT.to_map_geometry(tgt_prj=tgt_prj)

        return GeoArray(data_mapgeo, geotransform=gt, projection=prj)

    def to_sensor_geometry(self,
                           lons: np.ndarray,
                           lats: np.ndarray):
        GT = Geometry_Transformer(self.dem, lons=lons, lats=lats, nprocs=self.CPUs)
        data_sensorgeo = GT.to_sensor_geometry()

        return GeoArray(data_sensorgeo)
