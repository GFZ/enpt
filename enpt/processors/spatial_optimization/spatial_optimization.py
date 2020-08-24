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

"""EnPT spatial optimization module.

Adapts the EnMAP image geometry to a given Sentinel-2 L2A dataset.
Fits the VNIR detector data to the reference image. Corrects for keystone.
"""

__author__ = 'Daniel Scheffler'

import numpy as np

from arosics import COREG_LOCAL
from geoarray import GeoArray
from py_tools_ds.geo.coord_trafo import reproject_shapelyGeometry

from ...options.config import EnPTConfig
from ...model.images.images_sensorgeo import EnMAPL1Product_SensorGeo
from ..spatial_transform import Geometry_Transformer


class Spatial_Optimizer(object):
    def __init__(self, config: EnPTConfig):
        """Create an instance of Spatial_Optimizer"""
        self.cfg = config

    def _get_enmap_band_for_matching(self,
                                     enmap_ImageL1: EnMAPL1Product_SensorGeo)\
            -> GeoArray:
        """

        :param enmap_ImageL1:
        :return:
        """
        enmap_ImageL1.logger.warning('Statically using band 40 for co-registration.')
        bandidx = 39
        enmap_band_sensorgeo = enmap_ImageL1.vnir.data[:, :, bandidx]

        # transform from sensor to map geometry to make it usable for tie point detection
        enmap_ImageL1.logger.info('Temporarily transforming EnMAP band %d to map geometry for co-registration.'
                                  % (bandidx + 1))
        GT = Geometry_Transformer(lons=enmap_ImageL1.meta.vnir.lons[:, :, bandidx],
                                  lats=enmap_ImageL1.meta.vnir.lats[:, :, bandidx],
                                  fill_value=0,
                                  resamp_alg='gauss',
                                  radius_of_influence=30,
                                  nprocs=self.cfg.CPUs)

        ref_gA = GeoArray(self.cfg.path_reference_image)

        enmap_band_mapgeo = \
            GeoArray(*GT.to_map_geometry(enmap_band_sensorgeo,
                                         tgt_prj=ref_gA.prj,  # TODO correct?
                                         tgt_coordgrid=ref_gA.xygrid_specs),
                     nodata=0)

        return enmap_band_mapgeo

    def _get_enmap_mask_for_matching(self,
                                     enmap_ImageL1: EnMAPL1Product_SensorGeo)\
            -> GeoArray:
        """

        :param enmap_ImageL1:
        :return:
        """
        # use the water mask
        enmap_mask_sensorgeo = enmap_ImageL1.vnir.mask_landwater[:] == 2  # 2 is water here

        # transform from sensor to map geometry to make it usable for tie point detection
        enmap_ImageL1.logger.info('Temporarily transforming EnMAP water mask to map geometry for co-registration.')
        GT = Geometry_Transformer(lons=enmap_ImageL1.meta.vnir.lons[:, :, 39],  # FIXME hardcoded
                                  lats=enmap_ImageL1.meta.vnir.lats[:, :, 39],
                                  fill_value=0,
                                  resamp_alg='nearest',
                                  nprocs=self.cfg.CPUs)

        ref_gA = GeoArray(self.cfg.path_reference_image)

        enmap_mask_mapgeo = \
            GeoArray(*GT.to_map_geometry(enmap_mask_sensorgeo,
                                         tgt_prj=ref_gA.prj,  # TODO correct?
                                         tgt_coordgrid=ref_gA.xygrid_specs),
                     nodata=0)

        return enmap_mask_mapgeo

    def _compute_tie_points(self,
                            enmap_ImageL1: EnMAPL1Product_SensorGeo):
        enmap_band_mapgeo = self._get_enmap_band_for_matching(enmap_ImageL1)  # in the projection of the reference image
        ref_gA = GeoArray(self.cfg.path_reference_image)

        CRL = COREG_LOCAL(self.cfg.path_reference_image,
                          enmap_band_mapgeo,
                          grid_res=40,
                          max_shift=5,
                          nodata=(ref_gA.nodata, 0),
                          footprint_poly_tgt=reproject_shapelyGeometry(enmap_ImageL1.meta.vnir.ll_mapPoly,
                                                                       4326, enmap_band_mapgeo.epsg),
                          mask_baddata_tgt=self._get_enmap_mask_for_matching(enmap_ImageL1)
                          )
        TPG = CRL.tiepoint_grid
        # CRL.view_CoRegPoints(shapes2plot='vectors', hide_filtered=False, figsize=(20, 20),
        #                      savefigPath='/home/gfz-fe/scheffler/temp/EnPT/Archachon_AROSICS_tiepoints.png')

        valid_tiepoints = TPG.CoRegPoints_table[TPG.CoRegPoints_table.OUTLIER.__eq__(False)].copy()

        return valid_tiepoints

    def _interpolate_tiepoints_into_space(self, tiepoints, outshape, metric='ABS_SHIFT'):
        rows = np.array(tiepoints.Y_IM)
        cols = np.array(tiepoints.X_IM)
        data = np.array(tiepoints[metric])

        from time import time
        t0 = time()

        # https://github.com/agile-geoscience/xlines/blob/master/notebooks/11_Gridding_map_data.ipynb

        from scipy.interpolate import Rbf
        # f = Rbf(cols, rows, data, function='linear')
        # f = Rbf(cols, rows, data)
        # data_full = f(*np.meshgrid(np.arange(outshape[1]),
        #                            np.arange(outshape[0])))

        # rows_lowres = np.arange(0, outshape[0] + 10, 10)
        # cols_lowres = np.arange(0, outshape[1] + 10, 10)
        rows_lowres = np.arange(0, outshape[0] + 5, 5)
        cols_lowres = np.arange(0, outshape[1] + 5, 5)
        f = Rbf(cols, rows, data)
        data_interp_lowres = f(*np.meshgrid(cols_lowres, rows_lowres))

        # https://stackoverflow.com/questions/24978052/interpolation-over-regular-grid-in-python
        # from sklearn.gaussian_process import GaussianProcess
        # gp = GaussianProcess(theta0=0.1, thetaL=.001, thetaU=1., nugget=0.01)
        # gp.fit(X=np.column_stack([rr[vals], cc[vals]]), y=M[vals])
        # rr_cc_as_cols = np.column_stack([rr.flatten(), cc.flatten()])
        # interpolated = gp.predict(rr_cc_as_cols).reshape(M.shape)

        from scipy.interpolate import RegularGridInterpolator
        RGI = RegularGridInterpolator(points=[cols_lowres, rows_lowres],
                                      values=data_interp_lowres,
                                      method='linear')
        rows_full = np.arange(outshape[0])
        cols_full = np.arange(outshape[1])
        data_full = RGI(np.dstack(np.meshgrid(cols_full, rows_full, indexing='ij')))

        print('interpolation runtime: %.2fs' % (time() - t0))

        from matplotlib import pyplot as plt
        plt.figure()
        im = plt.imshow(data_full)
        plt.colorbar(im)
        plt.scatter(cols, rows, c=data, edgecolors='black')
        plt.title(metric)
        plt.show()

        return data_full

    def optimize_geolayer(self,
                          enmap_ImageL1: EnMAPL1Product_SensorGeo):
        if self.cfg.enable_absolute_coreg:


            tiepoints = self._compute_tie_points(enmap_ImageL1)
            xshift_map = self._interpolate_tiepoints_into_space(tiepoints,
                                                                (1800, 1800),  # FIXME hardcoded
                                                                # metric='ABS_SHIFT'
                                                                metric='X_SHIFT_M')
            yshift_map = self._interpolate_tiepoints_into_space(tiepoints,
                                                                (1800, 1800),  # FIXME hardcoded
                                                                # metric='ABS_SHIFT'
                                                                metric='Y_SHIFT_M')
            a=1