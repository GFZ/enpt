# -*- coding: utf-8 -*-
"""
EnPT module 'orthorectification' for transforming an EnMAP image from sensor to map geometry
based on a pixel- and band-wise coordinate-layer (geolayer).
"""

from typing import Tuple  # noqa: F401

import numpy as np
from geoarray import GeoArray
from py_tools_ds.geo.coord_trafo import transform_any_prj
from py_tools_ds.geo.projection import EPSG2WKT, prj_equal

from ...options.config import EnPTConfig
from ...model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo
from ..spatial_transform import Geometry_Transformer, move_extent_to_EnMAP_grid


class Orthorectifier(object):
    def __init__(self, config: EnPTConfig = None):
        """Create an instance of Orthorectifier."""
        self.cfg = config

    @staticmethod
    def validate_input(enmap_ImageL1):
        # check type
        if not isinstance(enmap_ImageL1, EnMAPL1Product_SensorGeo):
            raise TypeError(enmap_ImageL1, "The Orthorectifier expects an instance of 'EnMAPL1Product_SensorGeo'."
                                           "Received a '%s' instance." % type(enmap_ImageL1))

        # check geolayer shapes
        datashape = enmap_ImageL1.vnir.data.shape
        msg = 'Expected a %s ' + 'geolayer shape of %s or %s. ' % (str(datashape[:2]), str(datashape)) + 'Received %s.'
        for detector in [enmap_ImageL1.vnir, enmap_ImageL1.swir]:
            for XY in [detector.detector_meta.lons, detector.detector_meta.lats]:
                if XY.shape not in [datashape, datashape[:2]]:
                    raise RuntimeError(msg % (detector.detector_name, str(XY.shape)))

    def run_transformation(self, enmap_ImageL1: EnMAPL1Product_SensorGeo) -> EnMAPL2Product_MapGeo:
        self.validate_input(enmap_ImageL1)

        enmap_ImageL1.logger.info('Starting orthorectification...')

        # get a new instance of EnMAPL2Product_MapGeo
        L2_obj = EnMAPL2Product_MapGeo(config=self.cfg)

        # geometric transformations #
        #############################

        # get target UTM zone and common extent   # TODO add optionally user defined UTM zone?
        tgt_epsg = self._get_tgt_UTMepsg(enmap_ImageL1)
        tgt_extent = self._get_common_extent(enmap_ImageL1, tgt_epsg, enmap_grid=True)

        # transform VNIR and SWIR to map geometry
        GT_vnir = Geometry_Transformer(lons=enmap_ImageL1.vnir.detector_meta.lons,
                                       lats=enmap_ImageL1.vnir.detector_meta.lats,
                                       nprocs=self.cfg.CPUs)
        vnir_mapgeo, vnir_gt, vnir_prj = GT_vnir.to_map_geometry(enmap_ImageL1.vnir.data[:],
                                                                 tgt_prj=tgt_epsg,
                                                                 tgt_extent=tgt_extent)

        GT_swir = Geometry_Transformer(lons=enmap_ImageL1.swir.detector_meta.lons,
                                       lats=enmap_ImageL1.swir.detector_meta.lats,
                                       nprocs=self.cfg.CPUs)
        swir_mapgeo, swir_gt, swir_prj = GT_swir.to_map_geometry(enmap_ImageL1.swir.data[:],
                                                                 tgt_prj=tgt_epsg,
                                                                 tgt_extent=tgt_extent)

        # combine VNIR and SWIR
        L2_obj.data = self._get_VNIR_SWIR_stack(vnir_mapgeo, swir_mapgeo, vnir_gt, swir_gt, vnir_prj, swir_prj)

        # TODO transform mask_clouds, mask_clouds_confidence, ac_errors

        # metadata adjustments #
        ########################

        # TODO

        return L2_obj

    @staticmethod
    def _get_tgt_UTMepsg(enmap_ImageL1: EnMAPL1Product_SensorGeo) -> int:
        # get center lon/lat
        lon = np.mean([np.min(enmap_ImageL1.vnir.detector_meta.lon_UL_UR_LL_LR),
                       np.max(enmap_ImageL1.vnir.detector_meta.lon_UL_UR_LL_LR)])
        lat = np.mean([np.min(enmap_ImageL1.vnir.detector_meta.lat_UL_UR_LL_LR),
                       np.max(enmap_ImageL1.vnir.detector_meta.lat_UL_UR_LL_LR)])

        zoneNr = int(1 + (lon + 180.0) / 6.0)
        isNorth = lat >= 0

        return int('326' + str(zoneNr)) if isNorth else int('327' + str(zoneNr))

    @staticmethod
    def _get_common_extent(enmap_ImageL1: EnMAPL1Product_SensorGeo,
                           tgt_epsg: int,
                           enmap_grid: bool = True) -> Tuple[float, float, float, float]:
        """Get common target extent for VNIR and SWIR.

        :param enmap_ImageL1:
        :param tgt_epsg:
        :param enmap_grid:
        :return:
        """
        vnir_coords_utm = np.array([transform_any_prj(EPSG2WKT(4326), EPSG2WKT(tgt_epsg), x, y)
                                    for x, y in zip(enmap_ImageL1.vnir.detector_meta.lon_UL_UR_LL_LR,
                                                    enmap_ImageL1.vnir.detector_meta.lat_UL_UR_LL_LR,)])
        swir_coords_utm = np.array([transform_any_prj(EPSG2WKT(4326), EPSG2WKT(tgt_epsg), x, y)
                                    for x, y in zip(enmap_ImageL1.swir.detector_meta.lon_UL_UR_LL_LR,
                                                    enmap_ImageL1.swir.detector_meta.lat_UL_UR_LL_LR, )])
        all_coords_utm = np.dstack([vnir_coords_utm, swir_coords_utm])

        common_extent = (all_coords_utm[:, 0, :].min(),  # xmin
                         all_coords_utm[:, 1, :].min(),  # ymin
                         all_coords_utm[:, 0, :].max(),  # xmax
                         all_coords_utm[:, 1, :].max())  # ymax

        if enmap_grid:
            common_extent = move_extent_to_EnMAP_grid(common_extent)

        return common_extent

    @staticmethod
    def _get_VNIR_SWIR_stack(vnir_data, swir_data, vnir_gt, swir_gt, vnir_prj, swir_prj):
        """Stack VNIR and SWIR bands with respect to their spectral overlap."""
        if vnir_gt != swir_gt:
            raise ValueError((vnir_gt, swir_gt), 'VNIR and SWIR geoinformation should be equal.')
        if not prj_equal(vnir_prj, swir_prj):
            raise ValueError((vnir_prj, swir_prj), 'VNIR and SWIR projection should be equal.')

        # TODO implement VNIR / SWIR fusion respecting their spectral overlap
        return GeoArray(np.dstack([vnir_data, swir_data]), geotransform=vnir_gt, projection=vnir_prj)


