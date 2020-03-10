# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
EnPT module 'orthorectification' for transforming an EnMAP image from sensor to map geometry
based on a pixel- and band-wise coordinate-layer (geolayer).
"""

from typing import Tuple, Union  # noqa: F401

import numpy as np
from mvgavg import mvgavg
from geoarray import GeoArray
from py_tools_ds.geo.coord_trafo import transform_any_prj
from py_tools_ds.geo.projection import EPSG2WKT, prj_equal

from ...options.config import EnPTConfig
from ...model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo
from ...model.metadata import EnMAP_Metadata_L2A_MapGeo
from ..spatial_transform import \
    Geometry_Transformer, \
    Geometry_Transformer_3D, \
    move_extent_to_EnMAP_grid, \
    get_UTMEPSG_from_LonLat_cornersXY

__author__ = 'Daniel Scheffler'


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
        L2_obj = EnMAPL2Product_MapGeo(config=self.cfg, logger=enmap_ImageL1.logger)

        # geometric transformations #
        #############################

        lons_vnir, lats_vnir = enmap_ImageL1.vnir.detector_meta.lons, enmap_ImageL1.vnir.detector_meta.lats
        lons_swir, lats_swir = enmap_ImageL1.swir.detector_meta.lons, enmap_ImageL1.swir.detector_meta.lats

        # get target UTM zone and common extent  # TODO add optionally user defined UTM zone?
        tgt_epsg = self._get_tgt_UTMepsg(enmap_ImageL1)
        tgt_extent = self._get_common_extent(enmap_ImageL1, tgt_epsg, enmap_grid=True)
        kw_init = dict(resamp_alg=self.cfg.ortho_resampAlg,
                       nprocs=self.cfg.CPUs,
                       radius_of_influence=30 if not self.cfg.ortho_resampAlg == 'bilinear' else 45)
        kw_trafo = dict(tgt_prj=tgt_epsg, tgt_extent=tgt_extent)

        # transform VNIR and SWIR to map geometry
        GeoTransformer = Geometry_Transformer if lons_vnir.ndim == 2 else Geometry_Transformer_3D

        enmap_ImageL1.logger.info('Orthorectifying VNIR data...')
        GT_vnir = GeoTransformer(lons=lons_vnir, lats=lats_vnir, **kw_init)
        vnir_mapgeo, vnir_gt, vnir_prj = GT_vnir.to_map_geometry(enmap_ImageL1.vnir.data[:], **kw_trafo)

        enmap_ImageL1.logger.info('Orthorectifying SWIR data...')
        GT_swir = GeoTransformer(lons=lons_swir, lats=lats_swir, **kw_init)
        swir_mapgeo, swir_gt, swir_prj = GT_swir.to_map_geometry(enmap_ImageL1.swir.data[:], **kw_trafo)

        # combine VNIR and SWIR
        enmap_ImageL1.logger.info('Merging VNIR and SWIR data...')
        L2_obj.data = self._get_VNIR_SWIR_stack(vnir_mapgeo, swir_mapgeo, vnir_gt, swir_gt, vnir_prj, swir_prj,
                                                enmap_ImageL1.meta.vnir.wvl_center, enmap_ImageL1.meta.swir.wvl_center)

        # TODO allow to set geolayer band to be used for warping of 2D arrays
        GT_2D = Geometry_Transformer(lons=lons_vnir if lons_vnir.ndim == 2 else lons_vnir[:, :, 0],
                                     lats=lats_vnir if lats_vnir.ndim == 2 else lats_vnir[:, :, 0],
                                     ** kw_init)

        # FIXME cloud mask applies to BOTH detectors
        enmap_ImageL1.logger.info('Orthorectifying cloud mask...')
        L2_obj.mask_clouds = GeoArray(*GT_2D.to_map_geometry(enmap_ImageL1.vnir.mask_clouds, **kw_trafo))

        # TODO transform mask_clouds_confidence, ac_errors, pixel masks

        # metadata adjustments #
        ########################

        enmap_ImageL1.logger.info('Generating L2A metadata...')
        L2_obj.meta = EnMAP_Metadata_L2A_MapGeo(config=self.cfg,
                                                meta_l1b=enmap_ImageL1.meta,
                                                dims_mapgeo=L2_obj.data.shape,
                                                logger=L2_obj.logger)
        L2_obj.meta.add_band_statistics(L2_obj.data)

        L2_obj.data.meta.band_meta['wavelength'] = list(L2_obj.meta.wvl_center)
        L2_obj.data.meta.band_meta['bandwidths'] = list(L2_obj.meta.fwhm)
        L2_obj.data.meta.global_meta['wavelength_units'] = 'nanometers'

        # Get the paths according information delivered in the metadata
        L2_obj.paths = L2_obj.get_paths(self.cfg.output_dir)

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

        :para enmap_ImageL1:
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

    def _get_VNIR_SWIR_stack(self,
                             vnir_data: Union[np.ndarray, GeoArray],
                             swir_data: Union[np.ndarray, GeoArray],
                             vnir_gt: tuple,
                             swir_gt: tuple,
                             vnir_prj: str,
                             swir_prj: str,
                             wvls_vnir: np.ndarray,
                             wvls_swir: np.ndarray
                             ) -> GeoArray:
        """Stack VNIR and SWIR bands with respect to their spectral overlap.

        :param vnir_data:
        :param swir_data:
        :param vnir_gt:
        :param swir_gt:
        :param vnir_prj:
        :param swir_prj:
        :param wvls_vnir:
        :param wvls_swir:
        :return:
        """
        if vnir_gt != swir_gt:
            raise ValueError((vnir_gt, swir_gt), 'VNIR and SWIR geoinformation should be equal.')
        if not prj_equal(vnir_prj, swir_prj):
            raise ValueError((vnir_prj, swir_prj), 'VNIR and SWIR projection should be equal.')

        if self.cfg.vswir_overlap_algorithm in ['order_by_wvl', 'average']:
            # get band index order
            wvls_vswir = np.hstack([wvls_vnir, wvls_swir])
            wvls_vswir_sorted = np.array(sorted(wvls_vswir))
            bandidx_order = np.array([np.argmin(np.abs(wvls_vswir - cwl))
                                      for cwl in wvls_vswir_sorted])

            # stack bands ordered by wavelengths
            data_stacked = np.dstack([vnir_data, swir_data])[:, :, bandidx_order]

            if self.cfg.vswir_overlap_algorithm == 'average':
                # get wavelenghts and indices of overlapping bands
                wvls_overlap_vnir = wvls_vnir[wvls_vnir > wvls_swir.min()]
                wvls_overlap_swir = wvls_swir[wvls_swir < wvls_vnir.max()]
                wvls_overlap_all = np.array(sorted(np.hstack([wvls_overlap_vnir,
                                                              wvls_overlap_swir])))
                bandidxs_overlap = np.array([np.argmin(np.abs(wvls_vswir_sorted - cwl))
                                             for cwl in wvls_overlap_all])

                # apply a spectral moving average to the overlapping VNIR/SWIR bands
                filterwidth = 3  # number of bands to be included in the averaging - must be an uneven number
                bandidxs2average = np.array([bandidxs_overlap.min() - int((filterwidth - 1) / 2)] +
                                            list(bandidxs_overlap) +
                                            [bandidxs_overlap.max() + int((filterwidth - 1) / 2)])
                data2average = data_stacked[:, :, bandidxs2average]
                data_stacked[:, :, bandidxs_overlap] = mvgavg(data2average,
                                                              n=filterwidth,
                                                              axis=2).astype(data_stacked.dtype)

        elif self.cfg.vswir_overlap_algorithm == 'vnir_only':
            # get band index order
            wvls_swir_cut = wvls_swir[wvls_swir > wvls_vnir.max()]
            # wvls_vswir_sorted = np.hstack([wvls_vnir, wvls_swir_cut])
            idx_swir_firstband = np.argmin(np.abs(wvls_swir - wvls_swir_cut.min()))

            # stack bands while removing overlapping SWIR bands
            data_stacked = np.dstack([vnir_data, swir_data[:, :, idx_swir_firstband:]])

        elif self.cfg.vswir_overlap_algorithm == 'swir_only':
            # get band index order
            wvls_vnir_cut = wvls_vnir[wvls_vnir < wvls_swir.min()]
            # wvls_vswir_sorted = np.hstack([wvls_vnir_cut, wvls_swir])
            idx_vnir_lastband = np.argmin(np.abs(wvls_vnir - wvls_vnir_cut.max()))

            # stack bands while removing overlapping VNIR bands
            data_stacked = np.dstack([vnir_data[:, :, :idx_vnir_lastband], swir_data])

        else:
            raise ValueError(self.cfg.vswir_overlap_algorithm)

        return GeoArray(data_stacked, geotransform=vnir_gt, projection=vnir_prj)
