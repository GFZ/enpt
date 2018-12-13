# -*- coding: utf-8 -*-
""""""

import numpy as np
from geoarray import GeoArray

from py_tools_ds.geo.raster.reproject import SensorMapGeometryTransformer
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict
from py_tools_ds.geo.coord_grid import find_nearest

from ...options.config import enmap_coordinate_grid


class Geometry_Transformer(SensorMapGeometryTransformer):
    # use Sentinel-2 grid (30m grid with origin at 0/0)
    # EnMAP geolayer contains pixel center coordinate

    def __init__(self, path_or_geoarray, lons, lats, **opts):
        self._data2transform = GeoArray(path_or_geoarray)

        # TODO check if lons/lats represent pixel coordinate centers (needed as UL coordinates)

        super(Geometry_Transformer, self).__init__(self._data2transform[:], lons, lats, **opts)

    @property
    def is_sensor_geo(self):
        return self._data2transform.gt is None or self._data2transform.gt == [0, 1, 0, 0, 0, -1] \
               or not self._data2transform.prj

    @property
    def is_map_geo(self):
        return not self.is_sensor_geo

    def to_sensor_geometry(self, src_prj=None, src_extent=None):
        if self.is_sensor_geo:
            raise RuntimeError('The dataset to be transformed into sensor geometry already represents sensor geometry.')

        return super(Geometry_Transformer, self).to_sensor_geometry(
            src_prj=src_prj or self._data2transform.prj,
            src_extent=src_extent or list(np.array(self._data2transform.box.boundsMap)[[0, 2, 1, 3]]))

    def to_map_geometry(self, tgt_prj, tgt_extent=None, tgt_res=None):
        if self.is_map_geo:
            raise RuntimeError('The dataset to be transformed into map geometry already represents map geometry.')

        proj4dict = proj4_to_dict(get_proj4info(proj=tgt_prj))

        if 'units' in proj4dict and proj4dict['units'] == 'm':
            if not tgt_res:
                tgt_res = (np.ptp(enmap_coordinate_grid['x']), np.ptp(enmap_coordinate_grid['x']))

            if not tgt_extent:
                # use the extent computed by compute_output_shape and move it to the EnMAP coordinate grid
                x_size, y_size, out_gt, out_prj, out_extent = \
                    super(Geometry_Transformer, self).compute_output_shape(tgt_prj=tgt_prj)

                xmin, ymin, xmax, ymax = out_extent
                tgt_xgrid = enmap_coordinate_grid['x']
                tgt_ygrid = enmap_coordinate_grid['y']
                tgt_xmin = find_nearest(tgt_xgrid, xmin, roundAlg='off', extrapolate=True)
                tgt_xmax = find_nearest(tgt_xgrid, xmax, roundAlg='on', extrapolate=True)
                tgt_ymin = find_nearest(tgt_ygrid, ymin, roundAlg='off', extrapolate=True)
                tgt_ymax = find_nearest(tgt_ygrid, ymax, roundAlg='on', extrapolate=True)
                tgt_extent = list((tgt_xmin, tgt_ymin, tgt_xmax, tgt_ymax))

        out_data, out_gt, out_prj = \
            super(Geometry_Transformer, self).to_map_geometry(tgt_prj, tgt_extent=tgt_extent, tgt_res=tgt_res)

        return out_data, out_gt, out_prj
