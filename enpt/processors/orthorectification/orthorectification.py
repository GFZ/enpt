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
from ..spatial_transform import \
    Geometry_Transformer, \
    Geometry_Transformer_3D, \
    move_extent_to_EnMAP_grid, \
    get_UTMEPSG_from_LonLat_cornersXY


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
        for detector in [enmap_ImageL1.vnir, enmap_ImageL1.swir]:
            for XY in [detector.detector_meta.lons, detector.detector_meta.lats]:
                datashape = detector.data.shape
                if XY.shape not in [datashape, datashape[:2]]:
                    raise RuntimeError('Expected a %s geolayer shape of %s or %s. Received %s.'
                                       % (detector.detector_name, str(datashape), str(datashape[:2]), str(XY.shape)))

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
        GeoTransformer = \
            Geometry_Transformer if enmap_ImageL1.vnir.detector_meta.lons.ndim == 2 else Geometry_Transformer_3D

        GT_vnir = GeoTransformer(lons=enmap_ImageL1.vnir.detector_meta.lons,
                                 lats=enmap_ImageL1.vnir.detector_meta.lats,
                                 nprocs=self.cfg.CPUs)
        vnir_mapgeo, vnir_gt, vnir_prj = GT_vnir.to_map_geometry(enmap_ImageL1.vnir.data[:],
                                                                 tgt_prj=tgt_epsg,
                                                                 tgt_extent=tgt_extent)

        GT_swir = GeoTransformer(lons=enmap_ImageL1.swir.detector_meta.lons,
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

        # TODO ncols, nrows
        # TODO add CWLs and FWHM, nodata value to L2_obj.data.metadata

        return L2_obj

    @staticmethod
    def _get_tgt_UTMepsg(enmap_ImageL1: EnMAPL1Product_SensorGeo) -> int:
        return get_UTMEPSG_from_LonLat_cornersXY(lons=enmap_ImageL1.vnir.detector_meta.lon_UL_UR_LL_LR,
                                                 lats=enmap_ImageL1.vnir.detector_meta.lat_UL_UR_LL_LR)

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
