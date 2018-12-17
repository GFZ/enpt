# -*- coding: utf-8 -*-
"""EnPT pre-processing module for digital elevation models."""

from typing import Union, Tuple  # noqa: F401
from multiprocessing import cpu_count
import numpy as np

from geoarray import GeoArray
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict
from py_tools_ds.geo.coord_trafo import reproject_shapelyGeometry, transform_any_prj
from py_tools_ds.geo.vector.topology import get_footprint_polygon, get_overlap_polygon

from ..spatial_transform import Geometry_Transformer


class DEM_Processor(object):
    def __init__(self, dem_path_geoarray: Union[str, GeoArray],
                 enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                 CPUs: int = None):
        self.dem = GeoArray(dem_path_geoarray)
        self.enmapIm_cornerCoords = enmapIm_cornerCoords
        self.CPUs = CPUs or cpu_count()

        self._validate_input()

    def _validate_input(self):
        # check geocoding of DEM
        if not self.dem.is_map_geo:
            raise ValueError((self.dem.gt, self.dem.prj),
                             'The provided digital elevation model has no valid geo-coding or projection.')

        # check if provided as WGS-84
        proj4dict = proj4_to_dict(get_proj4info(proj=self.dem.prj))
        if 'datum' not in proj4dict or proj4dict['datum'] != 'WGS84':
            raise ValueError(proj4dict,
                             "The digital elevation model must be provided with 'WGS84' as geographic datum.")

        # check overlap
        dem_ll_mapPoly = reproject_shapelyGeometry(self.dem.footprint_poly, prj_src=self.dem.epsg, prj_tgt=4326)
        enmapIm_ll_mapPoly = get_footprint_polygon(self.enmapIm_cornerCoords, fix_invalid=True)
        overlap_dict = get_overlap_polygon(dem_ll_mapPoly, enmapIm_ll_mapPoly)
        overlap_perc = overlap_dict['overlap percentage']

        if overlap_perc < 100:
            # compute minimal extent in user provided projection
            cornersXY = np.array([transform_any_prj(4326, self.dem.epsg, x, y) for x, y in self.enmapIm_cornerCoords])
            xmin, xmax = cornersXY[:, 0].min(), cornersXY[:, 0].max()
            ymin, ymax = cornersXY[:, 1].min(), cornersXY[:, 1].max()

            raise ValueError('The provided digital elevation model covers %.1f percent of the EnMAP image. It must '
                             'cover it completely.The minimal needed extent in the provided projection is: \n'
                             'xmin: %s, xmax: %s, ymin: %s, ymax: %s.' % (overlap_perc, xmin, xmax, ymin, ymax))

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
