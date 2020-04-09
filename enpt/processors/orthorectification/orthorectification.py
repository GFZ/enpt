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
from types import SimpleNamespace

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
    def __init__(self, config: EnPTConfig):
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

        # FIXME So far, the fill value is set to 0. Is this meaningful?
        enmap_ImageL1.logger.info('Orthorectifying VNIR data...')
        GT_vnir = GeoTransformer(lons=lons_vnir, lats=lats_vnir, fill_value=0, **kw_init)
        vnir_mapgeo_gA = GeoArray(*GT_vnir.to_map_geometry(enmap_ImageL1.vnir.data[:], **kw_trafo),
                                  nodata=0)

        enmap_ImageL1.logger.info('Orthorectifying SWIR data...')
        GT_swir = GeoTransformer(lons=lons_swir, lats=lats_swir, fill_value=0, **kw_init)
        swir_mapgeo_gA = GeoArray(*GT_swir.to_map_geometry(enmap_ImageL1.swir.data[:], **kw_trafo),
                                  nodata=0)

        # combine VNIR and SWIR
        enmap_ImageL1.logger.info('Merging VNIR and SWIR data...')
        L2_obj.data = VNIR_SWIR_Stacker(vnir=vnir_mapgeo_gA,
                                        swir=swir_mapgeo_gA,
                                        vnir_wvls=enmap_ImageL1.meta.vnir.wvl_center,
                                        swir_wvls=enmap_ImageL1.meta.swir.wvl_center
                                        ).compute_stack(algorithm=self.cfg.vswir_overlap_algorithm)

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
                                                wvls_l2a=L2_obj.data.meta.band_meta['wavelength'],
                                                dims_mapgeo=L2_obj.data.shape,
                                                logger=L2_obj.logger)
        L2_obj.meta.add_band_statistics(L2_obj.data)

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


class VNIR_SWIR_Stacker(object):
    def __init__(self,
                 vnir: GeoArray,
                 swir: GeoArray,
                 vnir_wvls: Union[list, np.ndarray],
                 swir_wvls: Union[list, np.ndarray])\
            -> None:
        """Get an instance of VNIR_SWIR_Stacker.

        :param vnir:
        :param swir:
        :param vnir_wvls:
        :param swir_wvls:
        """
        self.vnir = vnir
        self.swir = swir
        self.wvls = SimpleNamespace(vnir=vnir_wvls, swir=swir_wvls)

        self.wvls.vswir = np.hstack([self.wvls.vnir, self.wvls.swir])
        self.wvls.vswir_sorted = np.array(sorted(self.wvls.vswir))

        self._validate_input()

    def _validate_input(self):
        if self.vnir.gt != self.swir.gt:
            raise ValueError((self.vnir.gt, self.swir.gt), 'VNIR and SWIR geoinformation should be equal.')
        if not prj_equal(self.vnir.prj, self.swir.prj):
            raise ValueError((self.vnir.prj, self.swir.prj), 'VNIR and SWIR projection should be equal.')
        if self.vnir.bands != len(self.wvls.vnir):
            raise ValueError("The number of VNIR bands must be equal to the number of elements in 'vnir_wvls': "
                             "%d != %d" % (self.vnir.bands, len(self.wvls.vnir)))
        if self.swir.bands != len(self.wvls.swir):
            raise ValueError("The number of SWIR bands must be equal to the number of elements in 'swir_wvls': "
                             "%d != %d" % (self.swir.bands, len(self.wvls.swir)))

    def _get_stack_order_by_wvl(self) -> Tuple[np.ndarray, np.ndarray]:
        """Stack bands ordered by wavelengths."""
        bandidx_order = np.array([np.argmin(np.abs(self.wvls.vswir - cwl))
                                  for cwl in self.wvls.vswir_sorted])

        return np.dstack([self.vnir[:], self.swir[:]])[:, :, bandidx_order], self.wvls.vswir_sorted

    def _get_stack_average(self, filterwidth: int = 3) -> Tuple[np.ndarray, np.ndarray]:
        """Stack bands and use averaging to compute the spectral information in the VNIR/SWIR overlap.

        :param filterwidth:     number of bands to be included in the averaging - must be an uneven number
        """
        # FIXME this has to respect nodata values - especially for pixels where one detector has no data.
        data_stacked = self._get_stack_order_by_wvl()[0]

        # get wavelenghts and indices of overlapping bands
        wvls_overlap_vnir = self.wvls.vnir[self.wvls.vnir > self.wvls.swir.min()]
        wvls_overlap_swir = self.wvls.swir[self.wvls.swir < self.wvls.vnir.max()]
        wvls_overlap_all = np.array(sorted(np.hstack([wvls_overlap_vnir,
                                                      wvls_overlap_swir])))
        bandidxs_overlap = np.array([np.argmin(np.abs(self.wvls.vswir_sorted - cwl))
                                     for cwl in wvls_overlap_all])

        # apply a spectral moving average to the overlapping VNIR/SWIR band
        bandidxs2average = np.array([bandidxs_overlap.min() - int((filterwidth - 1) / 2)] +
                                    list(bandidxs_overlap) +
                                    [bandidxs_overlap.max() + int((filterwidth - 1) / 2)])
        data2average = data_stacked[:, :, bandidxs2average]
        data_stacked[:, :, bandidxs_overlap] = mvgavg(data2average,
                                                      n=filterwidth,
                                                      axis=2).astype(data_stacked.dtype)

        return data_stacked, self.wvls.vswir_sorted

    def _get_stack_vnir_only(self) -> Tuple[np.ndarray, np.ndarray]:
        """Stack bands while removing overlapping SWIR bands."""
        wvls_swir_cut = self.wvls.swir[self.wvls.swir > self.wvls.vnir.max()]
        wvls_vswir_sorted = np.hstack([self.wvls.vnir, wvls_swir_cut])
        idx_swir_firstband = np.argmin(np.abs(self.wvls.swir - wvls_swir_cut.min()))

        return np.dstack([self.vnir[:], self.swir[:, :, idx_swir_firstband:]]), wvls_vswir_sorted

    def _get_stack_swir_only(self) -> Tuple[np.ndarray, np.ndarray]:
        """Stack bands while removing overlapping VNIR bands."""
        wvls_vnir_cut = self.wvls.vnir[self.wvls.vnir < self.wvls.swir.min()]
        wvls_vswir_sorted = np.hstack([wvls_vnir_cut, self.wvls.swir])
        idx_vnir_lastband = np.argmin(np.abs(self.wvls.vnir - wvls_vnir_cut.max()))

        return np.dstack([self.vnir[:, :, :idx_vnir_lastband + 1], self.swir[:]]), wvls_vswir_sorted

    def compute_stack(self, algorithm: str) -> GeoArray:
        """Stack VNIR and SWIR bands with respect to their spectral overlap.

        :param algorithm:   'order_by_wvl': keep spectral bands unchanged, order bands by wavelength
                            'average':      average the spectral information within the overlap
                            'vnir_only':    only use the VNIR bands (cut overlapping SWIR bands)
                            'swir_only':    only use the SWIR bands (cut overlapping VNIR bands)
        :return:    the stacked data cube as GeoArray instance
        """
        # TODO: This should also set an output nodata value.
        if algorithm == 'order_by_wvl':
            data_stacked, wvls = self._get_stack_order_by_wvl()
        elif algorithm == 'average':
            data_stacked, wvls = self._get_stack_average()
        elif algorithm == 'vnir_only':
            data_stacked, wvls = self._get_stack_vnir_only()
        elif algorithm == 'swir_only':
            data_stacked, wvls = self._get_stack_swir_only()
        else:
            raise ValueError(algorithm)

        gA_stacked = GeoArray(data_stacked, geotransform=self.vnir.gt, projection=self.vnir.prj)
        gA_stacked.meta.band_meta['wavelength'] = list(wvls)

        return gA_stacked
