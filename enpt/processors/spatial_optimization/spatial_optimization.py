# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""EnPT spatial optimization module.

Adapts the EnMAP image geometry to a given Sentinel-2 L2A dataset.
Fits the VNIR detector data to the reference image. Corrects for keystone.
"""

__author__ = 'Daniel Scheffler'

import os
from typing import Optional

import numpy as np
from osgeo import gdal  # noqa

from arosics import COREG_LOCAL
from geoarray import GeoArray
from py_tools_ds.geo.coord_trafo import reproject_shapelyGeometry, transform_coordArray, get_utm_zone
from py_tools_ds.geo.projection import EPSG2WKT
from py_tools_ds.processing.progress_mon import ProgressBar

from ...options.config import EnPTConfig
from ...model.images.images_sensorgeo import EnMAPL1Product_SensorGeo
from ..spatial_transform import Geometry_Transformer


class Spatial_Optimizer(object):
    def __init__(self, config: EnPTConfig):
        """Create an instance of Spatial_Optimizer."""
        self.cfg = config
        self._ref_Im: Optional[GeoArray, None] = GeoArray(self.cfg.path_reference_image)

        self._EnMAP_Im: Optional[EnMAPL1Product_SensorGeo, None] = None
        self._EnMAP_band: Optional[GeoArray, None] = None
        self._EnMAP_mask: Optional[GeoArray, None] = None
        self._ref_band_prep: Optional[GeoArray, None] = None
        self._EnMAP_bandIdx = 39  # FIXME hardcoded

    def _get_enmap_band_for_matching(self) \
            -> GeoArray:
        """Return the EnMAP band to be used in co-registration in UTM projection at 15m resolution."""
        self._EnMAP_Im.logger.warning(f'Statically using band {self._EnMAP_bandIdx + 1} for co-registration.')

        enmap_band_sensorgeo = self._EnMAP_Im.vnir.data[:, :, self._EnMAP_bandIdx]

        # transform from sensor to map geometry to make it usable for tie point detection
        # -> co-registration runs at 15m resolution to minimize information loss from resampling
        # -> a UTM projection is used to have shift length in meters
        self._EnMAP_Im.logger.info('Temporarily transforming EnMAP band %d to map geometry for co-registration.'
                                   % (self._EnMAP_bandIdx + 1))
        GT = Geometry_Transformer(lons=self._EnMAP_Im.meta.vnir.lons[:, :, self._EnMAP_bandIdx],
                                  lats=self._EnMAP_Im.meta.vnir.lats[:, :, self._EnMAP_bandIdx],
                                  backend='gdal',
                                  resamp_alg='bilinear',
                                  nprocs=self.cfg.CPUs)

        self._EnMAP_band = \
            GeoArray(*GT.to_map_geometry(enmap_band_sensorgeo,
                                         tgt_prj=self._get_optimal_utm_epsg(),
                                         tgt_coordgrid=((0, 15), (0, -15)),
                                         src_nodata=self._EnMAP_Im.vnir.data.nodata,
                                         tgt_nodata=0
                                         ),
                     nodata=0)

        return self._EnMAP_band

    def _get_enmap_mask_for_matching(self) \
            -> GeoArray:
        """Return the EnMAP mask to be used in co-registration in UTM projection at 15m resolution."""
        # use the water mask
        enmap_mask_sensorgeo = self._EnMAP_Im.vnir.mask_landwater[:] == 2  # 2 is water here

        # transform from sensor to map geometry to make it usable for tie point detection
        self._EnMAP_Im.logger.info('Temporarily transforming EnMAP water mask to map geometry for co-registration.')
        GT = Geometry_Transformer(lons=self._EnMAP_Im.meta.vnir.lons[:, :, self._EnMAP_bandIdx],
                                  lats=self._EnMAP_Im.meta.vnir.lats[:, :, self._EnMAP_bandIdx],
                                  backend='gdal',
                                  resamp_alg='nearest',
                                  nprocs=self.cfg.CPUs)

        self._EnMAP_mask = \
            GeoArray(*GT.to_map_geometry(enmap_mask_sensorgeo,
                                         tgt_prj=self._get_optimal_utm_epsg(),
                                         tgt_coordgrid=((0, 15), (0, -15)),
                                         src_nodata=0,  # 0=background
                                         tgt_nodata=0
                                         ),
                     nodata=0)

        return self._EnMAP_mask

    def _get_optimal_utm_epsg(self) -> int:
        """Return the EPSG code of the UTM zone that optimally covers the given EnMAP image."""
        x, y = [c.tolist()[0] for c in self._EnMAP_Im.meta.vnir.ll_mapPoly.centroid.xy]

        return int(f'{326 if y > 0 else 327}{get_utm_zone(x)}')

    @staticmethod
    def _get_suitable_nodata_value(arr: np.ndarray):
        """Get a suitable nodata value for the input array, which is not contained in the array data."""
        dtype = str(np.dtype(arr.dtype))
        try:
            # use a suitable nodata value
            nodata = dict(int8=-128, uint8=0, int16=-9999, uint16=9999, int32=-9999,
                          uint32=9999, float32=-9999., float64=-9999.)[dtype]
            if nodata not in arr:
                return nodata
            else:
                # use a suitable alternative nodata value
                alt_nodata = dict(int8=127, uint8=255, int16=32767, uint16=65535, int32=65535,
                                  uint32=65535, float32=65535, float64=65535)[dtype]
                if alt_nodata not in arr:
                    return alt_nodata
        except AttributeError:
            return None

    def _get_reference_band_for_matching(self) \
            -> GeoArray:
        """Return the reference image band to be used in co-registration in UTM projection at 15m resolution."""
        if self._EnMAP_band is None:
            raise RuntimeError('To prepare the reference image band, '
                               'the EnMAP band for matching needs to be computed first.')

        self._EnMAP_Im.logger.info('Preparing reference image for co-registration.')

        try:
            # get input nodata value if there is one defined in the metadata
            src_nodata = GeoArray(self.cfg.path_reference_image)._nodata  # noqa

            # only use the first band of the reference image for co-registration
            # TODO: select the most similar wavelength if CWLs are contained in the metadata
            #       and provide a user option to specify the band
            src_ds = gdal.Translate('', self.cfg.path_reference_image, format='MEM', bandList=[1])

            # set up the target dataset and coordinate grid
            driver = gdal.GetDriverByName('MEM')
            rows, cols = self._EnMAP_band.shape
            dst_ds = driver.Create('', cols, rows, 1, gdal.GDT_Float32)
            dst_ds.SetProjection(self._EnMAP_band.prj)
            dst_ds.SetGeoTransform(self._EnMAP_band.gt)
            if src_nodata is not None:
                dst_ds.GetRasterBand(1).SetNoDataValue(int(src_nodata))
                dst_ds.GetRasterBand(1).Fill(int(src_nodata))

            # reproject the needed subset from the reference image to a UTM 15m grid (like self._EnMAP_band)
            cb = ProgressBar(prefix='Warping progress    ') if not self.cfg.disable_progress_bars else None
            gdal.SetConfigOption('GDAL_NUM_THREADS', f'{self.cfg.CPUs}')
            gdal.ReprojectImage(
                src_ds,
                dst_ds,
                src_ds.GetProjection(),
                self._EnMAP_band.prj, gdal.GRA_Cubic,
                callback=cb
            )
            out_data = dst_ds.GetRasterBand(1).ReadAsArray()

        finally:
            del src_ds
            del dst_ds
            gdal.SetConfigOption('GDAL_NUM_THREADS', None)

        self._ref_band_prep = GeoArray(out_data, self._EnMAP_band.gt, self._EnMAP_band.prj, nodata=src_nodata)

        # try to set a meaningful nodata value if it cannot be auto-detected
        # -> only needed for AROSICS where setting a nodata value avoids warnings
        if self._ref_band_prep.nodata is None:
            self._ref_band_prep.nodata = self._get_suitable_nodata_value(out_data)

        return self._ref_band_prep

    def _compute_tie_points(self):
        # avoid to run RANSAC within AROSICS if more than 50 lines of a gapfill image were appended
        # (since the RPC coefficients are computed for the main image, there may be decreasing geometry accuracy with
        #  increasing distance from the main image)
        if os.path.exists(self.cfg.path_l1b_enmap_image_gapfill) and \
                (self.cfg.n_lines_to_append is None or self.cfg.n_lines_to_append > 50):
            tieP_filter_level = 2
        else:
            tieP_filter_level = 3

        # compute tie points within AROSICS
        CRL = COREG_LOCAL(self._ref_band_prep,
                          self._EnMAP_band,
                          grid_res=40,
                          max_shift=10,  # 5 EnMAP pixels (co-registration is running at 15m UTM grid)
                          nodata=(self._ref_Im.nodata, 0),
                          footprint_poly_tgt=reproject_shapelyGeometry(self._EnMAP_Im.meta.vnir.ll_mapPoly,
                                                                       4326, self._EnMAP_band.epsg),
                          mask_baddata_tgt=self._EnMAP_mask,
                          tieP_filter_level=tieP_filter_level,
                          progress=self.cfg.disable_progress_bars is False
                          )
        TPG = CRL.tiepoint_grid
        # CRL.view_CoRegPoints(shapes2plot='vectors', hide_filtered=False, figsize=(20, 20),
        #                      savefigPath='/home/gfz-fe/scheffler/temp/EnPT/Archachon_AROSICS_tiepoints.png')

        valid_tiepoints = TPG.CoRegPoints_table[TPG.CoRegPoints_table.OUTLIER.__eq__(False)].copy()

        return valid_tiepoints

    @staticmethod
    def _interpolate_tiepoints_into_space(tiepoints, outshape, metric='ABS_SHIFT'):
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
                                      values=data_interp_lowres.T,  # must be in shape [x, y]
                                      method='linear',
                                      bounds_error=False)
        rows_full = np.arange(outshape[0])
        cols_full = np.arange(outshape[1])
        data_full = RGI(np.dstack(np.meshgrid(cols_full, rows_full)))

        print('interpolation runtime: %.2fs' % (time() - t0))

        # from matplotlib import pyplot as plt
        # plt.figure()
        # im = plt.imshow(data_full)
        # plt.colorbar(im)
        # plt.scatter(cols, rows, c=data, edgecolors='black')
        # plt.title(metric)
        # plt.show()

        return data_full

    def optimize_geolayer(self,
                          enmap_ImageL1: EnMAPL1Product_SensorGeo):
        self._EnMAP_Im = enmap_ImageL1
        self._get_enmap_band_for_matching()
        self._get_enmap_mask_for_matching()
        self._get_reference_band_for_matching()

        enmap_ImageL1.logger.info('Computing tie points between the EnMAP image and the given spatial reference image.')
        tiepoints = self._compute_tie_points()

        enmap_ImageL1.logger.info('Generating misregistration array.')
        xshift_map = self._interpolate_tiepoints_into_space(tiepoints,
                                                            self._EnMAP_band.shape,
                                                            metric='X_SHIFT_M')
        yshift_map = self._interpolate_tiepoints_into_space(tiepoints,
                                                            self._EnMAP_band.shape,
                                                            metric='Y_SHIFT_M')

        ULx, ULy = self._EnMAP_band.box.boxMapXY[0]
        xgsd, ygsd = self._EnMAP_band.xgsd, self._EnMAP_band.ygsd
        rows, cols = self._EnMAP_band.shape
        xgrid_map, ygrid_map = np.meshgrid(np.arange(ULx, ULx + cols * xgsd, xgsd),
                                           np.arange(ULy, ULy - rows * ygsd, -ygsd))

        xgrid_map_coreg = xgrid_map + xshift_map
        ygrid_map_coreg = ygrid_map + yshift_map

        # transform map to sensor geometry
        enmap_ImageL1.logger.info('Transforming spatial optimization results back to sensor geometry.')
        lons_band = self._EnMAP_Im.meta.vnir.lons[:, :, self._EnMAP_bandIdx]
        lats_band = self._EnMAP_Im.meta.vnir.lats[:, :, self._EnMAP_bandIdx]
        GT = Geometry_Transformer(lons=lons_band,
                                  lats=lats_band,
                                  backend='gdal',
                                  resamp_alg='bilinear',
                                  nprocs=self.cfg.CPUs)

        geolayer_sensorgeo = \
            GT.to_sensor_geometry(GeoArray(np.dstack([xgrid_map_coreg,
                                                      ygrid_map_coreg]),
                                           geotransform=self._EnMAP_band.gt,
                                           projection=self._EnMAP_band.prj),
                                  tgt_nodata=0)

        enmap_ImageL1.logger.info('Applying results of spatial optimization to geolayer.')
        lons_coreg, lats_coreg = transform_coordArray(prj_src=self._ref_band_prep.prj,
                                                      prj_tgt=EPSG2WKT(4326),
                                                      Xarr=geolayer_sensorgeo[:, :, 0],
                                                      Yarr=geolayer_sensorgeo[:, :, 1])

        diffs_lons_coreg = lons_band - lons_coreg
        diffs_lats_coreg = lats_band - lats_coreg

        # enmap_ImageL1.meta.vnir.lons -= diffs_lons_coreg[:, :, np.newaxis]
        # enmap_ImageL1.meta.vnir.lats -= diffs_lats_coreg[:, :, np.newaxis]
        # enmap_ImageL1.meta.swir.lons -= diffs_lons_coreg[:, :, np.newaxis]
        # enmap_ImageL1.meta.swir.lats -= diffs_lats_coreg[:, :, np.newaxis]
        enmap_ImageL1.meta.vnir.lons = enmap_ImageL1.meta.vnir.lons - diffs_lons_coreg[:, :, np.newaxis]
        enmap_ImageL1.meta.vnir.lats = enmap_ImageL1.meta.vnir.lats - diffs_lats_coreg[:, :, np.newaxis]
        enmap_ImageL1.meta.swir.lons = enmap_ImageL1.meta.swir.lons - diffs_lons_coreg[:, :, np.newaxis]
        enmap_ImageL1.meta.swir.lats = enmap_ImageL1.meta.swir.lats - diffs_lats_coreg[:, :, np.newaxis]

        return enmap_ImageL1
