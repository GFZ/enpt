# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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

"""EnPT module 'orthorectification' for transforming an EnMAP image from sensor to map geometry
based on a pixel- and band-wise coordinate-layer (geolayer).
"""


from typing import Tuple, Union  # noqa: F401
from types import SimpleNamespace

import numpy as np
from mvgavg import mvgavg
from geoarray import GeoArray
from py_tools_ds.geo.coord_trafo import transform_any_prj
from py_tools_ds.geo.projection import prj_equal

from ...options.config import EnPTConfig
from ...model.images import EnMAPL1Product_SensorGeo, EnMAPL2Product_MapGeo
from ...model.metadata import EnMAP_Metadata_L2A_MapGeo
from ..spatial_transform import \
    Geometry_Transformer, \
    Geometry_Transformer_3D, \
    move_extent_to_coord_grid

__author__ = 'Daniel Scheffler'


class Orthorectifier(object):
    def __init__(self, config: EnPTConfig):
        """Create an instance of Orthorectifier."""
        self.cfg = config

    @staticmethod
    def validate_input(enmap_ImageL1: EnMAPL1Product_SensorGeo):
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

        # get target EPSG code and common extent
        tgt_epsg = enmap_ImageL1.meta.vnir.epsg_ortho
        tgt_extent = self._get_common_extent(enmap_ImageL1, tgt_epsg, enmap_grid=True)
        kw_init = dict(resamp_alg=self.cfg.ortho_resampAlg,
                       nprocs=self.cfg.CPUs,
                       radius_of_influence=30 if not self.cfg.ortho_resampAlg == 'bilinear' else 45
                       )
        if self.cfg.ortho_resampAlg == 'bilinear':
            # increase that if the resampling result contains gaps (default is 32 but this is quite slow)
            kw_init['neighbours'] = 8

        kw_trafo = dict(tgt_prj=tgt_epsg, tgt_extent=tgt_extent,
                        tgt_coordgrid=((self.cfg.target_coord_grid['x'], self.cfg.target_coord_grid['y'])
                                       if self.cfg.target_coord_grid else None)
                        )

        # transform VNIR and SWIR to map geometry
        GeoTransformer = Geometry_Transformer if lons_vnir.ndim == 2 else Geometry_Transformer_3D

        # FIXME So far, the fill value is set to 0. Is this meaningful?
        enmap_ImageL1.logger.info("Orthorectifying VNIR data using '%s' resampling algorithm..."
                                  % self.cfg.ortho_resampAlg)
        GT_vnir = GeoTransformer(lons=lons_vnir, lats=lats_vnir, fill_value=0, **kw_init)
        vnir_mapgeo_gA = GeoArray(*GT_vnir.to_map_geometry(enmap_ImageL1.vnir.data[:], **kw_trafo),
                                  nodata=0)

        enmap_ImageL1.logger.info("Orthorectifying SWIR data using '%s' resampling algorithm..."
                                  % self.cfg.ortho_resampAlg)
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

        # transform masks #
        ###################

        # TODO allow to set geolayer band to be used for warping of 2D arrays
        GT_2D = Geometry_Transformer(lons=lons_vnir if lons_vnir.ndim == 2 else lons_vnir[:, :, 0],
                                     lats=lats_vnir if lats_vnir.ndim == 2 else lats_vnir[:, :, 0],
                                     ** kw_init)  # FIXME bilinear resampling for masks with discrete values?

        for attrName in ['mask_landwater', 'mask_clouds', 'mask_cloudshadow', 'mask_haze', 'mask_snow', 'mask_cirrus']:
            attr = getattr(enmap_ImageL1.vnir, attrName)

            if attr is not None:
                enmap_ImageL1.logger.info("Orthorectifying '%s' attribute..." % attrName)
                attr_ortho = GeoArray(*GT_2D.to_map_geometry(attr, **kw_trafo), nodata=attr.nodata)
                setattr(L2_obj, attrName, attr_ortho)

        # TODO transform dead pixel map, quality test flags?

        # set all pixels to nodata that don't have values in all bands #
        ################################################################

        enmap_ImageL1.logger.info("Setting all pixels to nodata that have values in the VNIR or the SWIR only...")
        mask_nodata_common = np.all(np.dstack([vnir_mapgeo_gA.mask_nodata[:],
                                               swir_mapgeo_gA.mask_nodata[:]]), axis=2)
        L2_obj.data[~mask_nodata_common] = L2_obj.data.nodata

        for attr_gA in [L2_obj.mask_landwater, L2_obj.mask_clouds, L2_obj.mask_cloudshadow, L2_obj.mask_haze,
                        L2_obj.mask_snow, L2_obj.mask_cirrus]:
            if attr_gA is not None:
                attr_gA[~mask_nodata_common] = attr_gA.nodata

        # metadata adjustments #
        ########################

        enmap_ImageL1.logger.info('Generating L2A metadata...')
        L2_obj.meta = EnMAP_Metadata_L2A_MapGeo(config=self.cfg,
                                                meta_l1b=enmap_ImageL1.meta,
                                                wvls_l2a=L2_obj.data.meta.band_meta['wavelength'],
                                                dims_mapgeo=L2_obj.data.shape,
                                                grid_res_l2a=(L2_obj.data.gt[1], abs(L2_obj.data.gt[5])),
                                                logger=L2_obj.logger)
        L2_obj.meta.add_band_statistics(L2_obj.data)

        L2_obj.data.meta.band_meta['fwhm'] = list(L2_obj.meta.fwhm)
        L2_obj.data.meta.global_meta['wavelength_units'] = 'nanometers'

        # Get the paths according information delivered in the metadata
        L2_obj.paths = L2_obj.get_paths(self.cfg.output_dir)

        return L2_obj

    def _get_common_extent(self,
                           enmap_ImageL1: EnMAPL1Product_SensorGeo,
                           tgt_epsg: int,
                           enmap_grid: bool = True) -> Tuple[float, float, float, float]:
        """Get common target extent for VNIR and SWIR.

        :para enmap_ImageL1:
        :param tgt_epsg:
        :param enmap_grid:
        :return:
        """
        # get geolayers - 2D for dummy data format else 3D
        V_lons, V_lats = enmap_ImageL1.meta.vnir.lons, enmap_ImageL1.meta.vnir.lats
        S_lons, S_lats = enmap_ImageL1.meta.swir.lons, enmap_ImageL1.meta.swir.lats

        # get Lon/Lat corner coordinates of geolayers
        V_UL_UR_LL_LR_ll = [(V_lons[y, x], V_lats[y, x]) for y, x in [(0, 0), (0, -1), (-1, 0), (-1, -1)]]
        S_UL_UR_LL_LR_ll = [(S_lons[y, x], S_lats[y, x]) for y, x in [(0, 0), (0, -1), (-1, 0), (-1, -1)]]

        # transform them to UTM
        V_UL_UR_LL_LR_prj = [transform_any_prj(4326, tgt_epsg, x, y) for x, y in V_UL_UR_LL_LR_ll]
        S_UL_UR_LL_LR_prj = [transform_any_prj(4326, tgt_epsg, x, y) for x, y in S_UL_UR_LL_LR_ll]

        # separate X and Y
        V_X_prj, V_Y_prj = zip(*V_UL_UR_LL_LR_prj)
        S_X_prj, S_Y_prj = zip(*S_UL_UR_LL_LR_prj)

        # in case of 3D geolayers, the corner coordinates have multiple values for multiple bands
        # -> use the innermost coordinates to avoid pixels with VNIR-only/SWIR-only values due to keystone
        #    (these pixels would be set to nodata later anyways, so we don't need to increase the extent for them)
        if V_lons.ndim == 3:
            V_X_prj = (V_X_prj[0].max(), V_X_prj[1].min(), V_X_prj[2].max(), V_X_prj[3].min())
            V_Y_prj = (V_Y_prj[0].min(), V_Y_prj[1].min(), V_Y_prj[2].max(), V_Y_prj[3].max())
            S_X_prj = (S_X_prj[0].max(), S_X_prj[1].min(), S_X_prj[2].max(), S_X_prj[3].min())
            S_Y_prj = (S_Y_prj[0].min(), S_Y_prj[1].min(), S_Y_prj[2].max(), S_Y_prj[3].max())

        # use inner coordinates of VNIR and SWIR as common extent
        xmin_prj = max([min(V_X_prj), min(S_X_prj)])
        ymin_prj = max([min(V_Y_prj), min(S_Y_prj)])
        xmax_prj = min([max(V_X_prj), max(S_X_prj)])
        ymax_prj = min([max(V_Y_prj), max(S_Y_prj)])
        common_extent_prj = (xmin_prj, ymin_prj, xmax_prj, ymax_prj)

        # move the extent to the EnMAP coordinate grid
        if enmap_grid and self.cfg.target_coord_grid:
            common_extent_prj = move_extent_to_coord_grid(common_extent_prj,
                                                          self.cfg.target_coord_grid['x'],
                                                          self.cfg.target_coord_grid['y'],)

        enmap_ImageL1.logger.info('Computed common target extent of orthorectified image (xmin, ymin, xmax, ymax in '
                                  'EPSG %s): %s' % (tgt_epsg, str(common_extent_prj)))

        return common_extent_prj


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
        bandidxs2average = np.array([np.min(bandidxs_overlap) - int((filterwidth - 1) / 2)] +
                                    list(bandidxs_overlap) +
                                    [np.max(bandidxs_overlap) + int((filterwidth - 1) / 2)])
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

        gA_stacked = GeoArray(data_stacked,
                              geotransform=self.vnir.gt, projection=self.vnir.prj, nodata=self.vnir.nodata)
        gA_stacked.meta.band_meta['wavelength'] = list(wvls)

        return gA_stacked
