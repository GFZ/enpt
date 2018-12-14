# -*- coding: utf-8 -*-
"""EnPT module 'spatial transform', containing everything related to spatial transformations."""

import numpy as np
from typing import Union, Tuple  # noqa: F401
from geoarray import GeoArray

from py_tools_ds.geo.raster.reproject import SensorMapGeometryTransformer, AreaDefinition
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict
from py_tools_ds.geo.coord_grid import find_nearest

from ...options.config import enmap_coordinate_grid


class Geometry_Transformer(SensorMapGeometryTransformer):
    # use Sentinel-2 grid (30m grid with origin at 0/0)
    # EnMAP geolayer contains pixel center coordinate

    def to_sensor_geometry(self,
                           path_or_geoarray_mapgeo: Union[str, GeoArray],
                           src_prj: Union[str, int] = None,
                           src_extent: Tuple[float, float, float, float] = None):
        data_mapgeo = GeoArray(path_or_geoarray_mapgeo)

        if not data_mapgeo.is_map_geo:
            raise RuntimeError('The dataset to be transformed into sensor geometry already represents sensor geometry.')

        return super(Geometry_Transformer, self).to_sensor_geometry(
            data_mapgeo[:],
            src_prj=src_prj or data_mapgeo.prj,
            src_extent=src_extent or list(np.array(data_mapgeo.box.boundsMap)[[0, 2, 1, 3]]))

    def to_map_geometry(self,
                        path_or_geoarray_sensorgeo: Union[str, GeoArray],
                        tgt_prj:  Union[str, int] = None,
                        tgt_extent: Tuple[float, float, float, float] = None,
                        tgt_res: Tuple[float, float] = None,
                        area_definition: AreaDefinition = None):
        data_sensorgeo = GeoArray(path_or_geoarray_sensorgeo)

        if data_sensorgeo.is_map_geo:
            raise RuntimeError('The dataset to be transformed into map geometry already represents map geometry.')

        if area_definition:
            self.area_definition = area_definition
        else:
            if not tgt_prj:
                raise ValueError(tgt_prj, 'Target projection must be given if area_definition is not given.')

            proj4dict = proj4_to_dict(get_proj4info(proj=tgt_prj))

            if 'units' in proj4dict and proj4dict['units'] == 'm':
                if not tgt_res:
                    tgt_res = (np.ptp(enmap_coordinate_grid['x']), np.ptp(enmap_coordinate_grid['x']))

                if not tgt_extent:
                    # use the extent computed by compute_output_shape and move it to the EnMAP coordinate grid
                    area_definition = self.compute_areadefinition_sensor2map(
                        data_sensorgeo[:], tgt_prj, tgt_res=tgt_res)

                    tgt_extent = move_extent_to_EnMAP_grid(tuple(area_definition.area_extent))

        out_data, out_gt, out_prj = \
            super(Geometry_Transformer, self).to_map_geometry(data_sensorgeo[:], tgt_prj=tgt_prj,
                                                              tgt_extent=tgt_extent, tgt_res=tgt_res,
                                                              area_definition=self.area_definition)

        return out_data, out_gt, out_prj


def move_extent_to_EnMAP_grid(extent_utm: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    """Move a UTM coordinate extent to the EnMAP coordinate grid (30m x 30m, origin at 0/0).

    :param extent_utm:  xmin, ymin, xmax, ymax coordinates
    """
    xmin, ymin, xmax, ymax = extent_utm
    tgt_xgrid = enmap_coordinate_grid['x']
    tgt_ygrid = enmap_coordinate_grid['y']
    tgt_xmin = find_nearest(tgt_xgrid, xmin, roundAlg='off', extrapolate=True)
    tgt_xmax = find_nearest(tgt_xgrid, xmax, roundAlg='on', extrapolate=True)
    tgt_ymin = find_nearest(tgt_ygrid, ymin, roundAlg='off', extrapolate=True)
    tgt_ymax = find_nearest(tgt_ygrid, ymax, roundAlg='on', extrapolate=True)

    return tgt_xmin, tgt_ymin, tgt_xmax, tgt_ymax
