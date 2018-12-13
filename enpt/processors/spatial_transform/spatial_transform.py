# -*- coding: utf-8 -*-
"""EnPT module 'spatial transform', containing everything related to spatial transformations."""

import numpy as np
from typing import Union, Tuple  # noqa: F401
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

    def to_sensor_geometry(self,
                           src_prj: Union[str, int] = None,
                           src_extent: Tuple[float, float, float, float] = None):
        if not self._data2transform.is_map_geo:
            raise RuntimeError('The dataset to be transformed into sensor geometry already represents sensor geometry.')

        return super(Geometry_Transformer, self).to_sensor_geometry(
            src_prj=src_prj or self._data2transform.prj,
            src_extent=src_extent or list(np.array(self._data2transform.box.boundsMap)[[0, 2, 1, 3]]))

    def to_map_geometry(self,
                        tgt_prj:  Union[str, int],
                        tgt_extent: Tuple[float, float, float, float] = None,
                        tgt_res: Tuple[float, float] = None):
        if self._data2transform.is_map_geo:
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
