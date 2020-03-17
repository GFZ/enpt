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

"""EnPT module 'spatial transform', containing everything related to spatial transformations."""

from typing import Union, Tuple, List  # noqa: F401
from multiprocessing import Pool, cpu_count
from collections import OrderedDict
import numpy as np
from scipy.interpolate import griddata as interpolate_griddata, interp1d
from geoarray import GeoArray

from sensormapgeo.sensormapgeo import \
    SensorMapGeometryTransformer, \
    SensorMapGeometryTransformer3D, \
    AreaDefinition
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict, prj_equal, EPSG2WKT, WKT2EPSG, proj4_to_WKT
from py_tools_ds.geo.coord_grid import find_nearest
from py_tools_ds.geo.coord_trafo import transform_any_prj, transform_coordArray

from ...options.config import enmap_coordinate_grid

__author__ = 'Daniel Scheffler'


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
                        path_or_geoarray_sensorgeo: Union[str, GeoArray, np.ndarray],
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

            # compute target resolution and extent (according to EnMAP grid)
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


class Geometry_Transformer_3D(SensorMapGeometryTransformer3D):
    # use Sentinel-2 grid (30m grid with origin at 0/0)
    # EnMAP geolayer contains pixel center coordinate

    def to_sensor_geometry(self,
                           path_or_geoarray_mapgeo: Union[str, GeoArray],
                           src_prj: Union[str, int] = None,
                           src_extent: Tuple[float, float, float, float] = None):
        data_mapgeo = GeoArray(path_or_geoarray_mapgeo)

        if not data_mapgeo.is_map_geo:
            raise RuntimeError('The dataset to be transformed into sensor geometry already represents sensor geometry.')

        return super(Geometry_Transformer_3D, self).to_sensor_geometry(
            data_mapgeo[:],
            src_prj=src_prj or data_mapgeo.prj,
            src_extent=src_extent or list(np.array(data_mapgeo.box.boundsMap)[[0, 2, 1, 3]]))

    def to_map_geometry(self,
                        path_or_geoarray_sensorgeo: Union[str, GeoArray, np.ndarray],
                        tgt_prj:  Union[str, int] = None,
                        tgt_extent: Tuple[float, float, float, float] = None,
                        tgt_res: Tuple[float, float] = None
                        ) -> Tuple[np.ndarray, tuple, str]:
        data_sensorgeo = GeoArray(path_or_geoarray_sensorgeo)

        if data_sensorgeo.is_map_geo:
            raise RuntimeError('The dataset to be transformed into map geometry already represents map geometry.')

        if not tgt_prj:
            raise ValueError(tgt_prj, 'Target projection must be given if area_definition is not given.')

        # compute target resolution and extent (according to EnMAP grid)
        proj4dict = proj4_to_dict(get_proj4info(proj=tgt_prj))

        if 'units' in proj4dict and proj4dict['units'] == 'm':
            if not tgt_res:
                tgt_res = (np.ptp(enmap_coordinate_grid['x']), np.ptp(enmap_coordinate_grid['x']))

            if not tgt_extent:
                # use the extent computed by compute_output_shape and move it to the EnMAP coordinate grid
                tgt_epsg = WKT2EPSG(proj4_to_WKT(get_proj4info(proj=tgt_prj)))
                common_extent = self._get_common_target_extent(tgt_epsg)

                tgt_extent = move_extent_to_EnMAP_grid(tuple(common_extent))

        out_data, out_gt, out_prj = \
            super(Geometry_Transformer_3D, self).to_map_geometry(data_sensorgeo[:], tgt_prj=tgt_prj,
                                                                 tgt_extent=tgt_extent, tgt_res=tgt_res)

        return out_data, out_gt, out_prj


class VNIR_SWIR_SensorGeometryTransformer(object):
    """Class to transform between EnMAP VNIR and SWIR sensor geometry."""

    def __init__(self,
                 lons_vnir: np.ndarray,
                 lats_vnir: np.ndarray,
                 lons_swir: np.ndarray,
                 lats_swir: np.ndarray,
                 prj_vnir: Union[str, int],
                 prj_swir: Union[str, int],
                 res_vnir: Tuple[float, float],
                 res_swir: Tuple[float, float],
                 resamp_alg: str = 'nearest',
                 **gt_opts) -> None:
        """Get an instance of VNIR_SWIR_SensorGeometryTransformer.

        :param lons_vnir:   2D VNIR longitude array corresponding to the sensor geometry arrays passed later
        :param lats_vnir:   2D VNIR latitude array corresponding to the sensor geometry arrays passed later
        :param lons_swir:   2D SWIR longitude array corresponding to the sensor geometry arrays passed later
        :param lats_swir:   2D SWIR latitude array corresponding to the sensor geometry arrays passed later
        :param prj_vnir:    projection of the VNIR if it would be transformed to map geometry (WKT string or EPSG code)
        :param prj_swir:    projection of the SWIR if it would be transformed to map geometry (WKT string or EPSG code)
        :param res_vnir:    pixel dimensions of the VNIR if it would be transformed to map geometry (X, Y)
        :param res_swir:    pixel dimensions of the SWIR if it would be transformed to map geometry (X, Y)
        :param resamp_alg:  resampling algorithm ('nearest', 'bilinear', 'gauss', 'custom')
        :param gt_opts:     additional options to be passed to the Geometric_Transformer class,
                            e.g., 'fill_value', 'radius_of_influence', ...
        """
        [self._validate_lonlat_ndim(ll) for ll in [lons_vnir, lats_vnir, lons_swir, lats_swir]]

        self.vnir_meta = dict(
            lons=lons_vnir,
            lats=lats_vnir,
            prj=prj_vnir,
            res=res_vnir
        )
        self.swir_meta = dict(
            lons=lons_swir,
            lats=lats_swir,
            prj=prj_swir,
            res=res_swir
        )
        self.resamp_alg = resamp_alg
        self.gt_opts = gt_opts

    @staticmethod
    def _validate_lonlat_ndim(coord_array):
        if coord_array.ndim == 3:
            raise RuntimeError("3D longitude/latitude array are not supported because they model keystone effects "
                               "which cannot be transferred between VNIR and SWIR.")
        else:
            assert coord_array.ndim == 2, 'Only 2D coordinate arrays are supported.'

    def transform_sensorgeo_VNIR_to_SWIR(self, data_vnirsensorgeo: np.ndarray) -> np.ndarray:
        """Transform any array in VNIR sensor geometry to SWIR sensor geometry to remove geometric shifts.

        :param data_vnirsensorgeo:      input array in VNIR sensor geometry
        :return:    input array resampled to SWIR sensor geometry
        """
        return self._transform_sensorgeo(data_vnirsensorgeo, inputgeo='vnir')

    def transform_sensorgeo_SWIR_to_VNIR(self, data_swirsensorgeo: np.ndarray) -> np.ndarray:
        """Transform any array in SWIR sensor geometry to VNIR sensor geometry to remove geometric shifts.

        :param data_swirsensorgeo:      input array in SWIR sensor geometry
        :return:    input array resampled to VNIR sensor geometry
        """
        return self._transform_sensorgeo(data_swirsensorgeo, inputgeo='swir')

    def _transform_sensorgeo(self,
                             data2transform: np.ndarray,
                             inputgeo: str = 'vnir') -> np.ndarray:
        # TODO: Avoid the resampling here, maybe by replacing the lon/lat arrays by image coordinates for the source
        #       geometry and by image coordinate differences for the target geometry. Maybe also the proj string for
        #       local coordinate systems helps (see SensorMapGeometryTransformer class).

        if inputgeo not in ['vnir', 'swir']:
            raise ValueError(inputgeo)

        src, tgt = (self.vnir_meta, self.swir_meta) if inputgeo == 'vnir' else (self.swir_meta, self.vnir_meta)

        # temporarily transform the input sensor geometry array to map geometry
        GT_src = Geometry_Transformer(lons=src['lons'], lats=src['lats'], resamp_alg=self.resamp_alg, **self.gt_opts)
        gA_mapgeo = GeoArray(*GT_src.to_map_geometry(data2transform, tgt_prj=src['prj']))

        # generate the target sensor geometry array (target lons/lats define the target swath definition)
        GT_tgt = Geometry_Transformer(lons=tgt['lons'], lats=tgt['lats'], resamp_alg=self.resamp_alg, **self.gt_opts)
        tgt_data_sensorgeo = GT_tgt.to_sensor_geometry(gA_mapgeo)

        return tgt_data_sensorgeo


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


class RPC_Geolayer_Generator(object):
    """
    Class for creating pixel-wise longitude/latitude arrays based on rational polynomial coefficients (RPC).
    """
    def __init__(self,
                 rpc_coeffs: dict,
                 dem: Union[str, GeoArray],
                 enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                 enmapIm_dims_sensorgeo: Tuple[int, int]):
        """Get an instance of RPC_Geolayer_Generator.

        :param rpc_coeffs:              RPC coefficients for a single EnMAP band
        :param dem:                     digital elevation model in map geometry (file path or instance of GeoArray)
        :param enmapIm_cornerCoords:    corner coordinates as tuple of lon/lat pairs
        :param enmapIm_dims_sensorgeo:  dimensions of the EnMAP image in sensor geometry (rows, colunms)
        """
        self.row_off = rpc_coeffs['row_off']
        self.col_off = rpc_coeffs['col_off']
        self.lat_off = rpc_coeffs['lat_off']
        self.lon_off = rpc_coeffs['long_off']
        self.height_off = rpc_coeffs['height_off']
        self.row_scale = rpc_coeffs['row_scale']
        self.col_scale = rpc_coeffs['col_scale']
        self.lat_scale = rpc_coeffs['lat_scale']
        self.lon_scale = rpc_coeffs['long_scale']
        self.height_scale = rpc_coeffs['height_scale']
        self.row_num_coeffs = rpc_coeffs['row_num_coeffs']
        self.row_den_coeffs = rpc_coeffs['row_den_coeffs']
        self.col_num_coeffs = rpc_coeffs['col_num_coeffs']
        self.col_den_coeffs = rpc_coeffs['col_den_coeffs']

        self.dem = GeoArray(dem)
        self.enmapIm_cornerCoords = enmapIm_cornerCoords
        self.enmapIm_dims_sensorgeo = enmapIm_dims_sensorgeo

    def _normalize_map_coordinates(self,
                                   lon: np.ndarray,
                                   lat: np.ndarray,
                                   height: np.ndarray) -> (np.ndarray, np.ndarray, np.ndarray):
        """Normalize map coordinates to [-1, +1] to improve numerical precision.

        :param lon:     longitude array
        :param lat:     latitude array
        :param height:  elevation array
        """
        if not lon.shape == lat.shape == height.shape:
            raise ValueError((lon.shape, lat.shape, height.shape),
                             'Longitude, latitude and height arrays are expected to have the same dimensions.')

        lon_norm = (lon - self.lon_off) / self.lon_scale  # longitude
        lat_norm = (lat - self.lat_off) / self.lat_scale  # latitude
        height_norm = (height - self.height_off) / self.height_scale  # elevation

        msg = 'Coordinate normalization yields significantly out-of-range values for %s. ' \
              'Check the coordinates and RPC coefficients.'
        # for llh, name in zip([lon_norm, lat_norm, height_norm], ['longitudes', 'latitudes', 'heights']):
        for llh, name in zip([lon_norm, lat_norm, height_norm], ['longitudes', 'latitudes']):
            if llh.min() < -1.1 or llh.max() > 1.1:
                raise RuntimeError((llh.min(), llh.max()), msg % name)

        return lon_norm, lat_norm, height_norm

    def _compute_normalized_image_coordinates(self,
                                              lon_norm: np.ndarray,
                                              lat_norm: np.ndarray,
                                              height_norm: np.ndarray) -> (np.ndarray, np.ndarray):
        """Compute normalized sensor geometry coordinates for the given lon/lat/height positions.

        :param lon_norm:    normalized longitudes
        :param lat_norm:    normalized latitudes
        :param height_norm: normalized elevations
        :return:
        """
        P = lat_norm.flatten()
        L = lon_norm.flatten()
        H = height_norm.flatten()

        u = np.zeros((P.size, 20))

        u_data = (ui for ui in [
            1,
            L,
            P,
            H,
            L * P,
            L * H,
            P * H,
            L ** 2,
            P ** 2,
            H ** 2,
            P * L * H,
            L ** 3,
            L * P ** 2,
            L * H ** 2,
            L ** 2 * P,
            P ** 3,
            P * H ** 2,
            L ** 2 * H,
            P ** 2 * H,
            H ** 3
        ])

        for i, ud in enumerate(u_data):
            u[:, i] = ud

        num_row_norm = np.sum(self.row_num_coeffs * u, axis=1)
        den_row_norm = np.sum(self.row_den_coeffs * u, axis=1)
        num_col_norm = np.sum(self.col_num_coeffs * u, axis=1)
        den_col_norm = np.sum(self.col_den_coeffs * u, axis=1)

        row_norm = num_row_norm / den_row_norm
        col_norm = num_col_norm / den_col_norm

        return row_norm, col_norm

    def _denormalize_image_coordinates(self,
                                       row_norm: np.ndarray,
                                       col_norm: np.ndarray) -> (np.ndarray, np.ndarray):
        """De-normalize norrmalized sensor geometry coordinates to get valid image coordinates.

        :param row_norm:    normalized rows
        :param col_norm:    normalized columns
        :return:    de-normalized rows array,  de-normalized columns array,
        """
        rows = row_norm * self.row_scale + self.row_off
        cols = col_norm * self.col_scale + self.col_off

        return rows, cols

    def transform_LonLatHeight_to_RowCol(self,
                                         lon: np.ndarray,
                                         lat: np.ndarray,
                                         height: np.ndarray) -> (np.ndarray, np.ndarray):
        """Get sensor geometry image coordinates for the given 3D map coordinate positions using RPC coefficients.

        :param lon:     longitude array
        :param lat:     latitude array
        :param height:  elevation array
        :return:    rows array, columns array
        """
        # TODO add reshaping

        lon_norm, lat_norm, height_norm = \
            self._normalize_map_coordinates(lon=lon, lat=lat, height=height)
        row_norm, col_norm = \
            self._compute_normalized_image_coordinates(lon_norm=lon_norm, lat_norm=lat_norm, height_norm=height_norm)
        rows, cols = \
            self._denormalize_image_coordinates(row_norm=row_norm, col_norm=col_norm)

        return rows, cols

    @staticmethod
    def _fill_nans_at_corners(arr: np.ndarray, along_axis: int = 0) -> np.ndarray:
        if not arr.ndim == 2:
            raise ValueError(arr.ndim, '2D numpy array expected.')
        if along_axis not in [0, 1]:
            raise ValueError(along_axis, "The 'axis' parameter must be set to 0 or 1")

        kw = dict(kind='linear', fill_value='extrapolate')

        if along_axis == 0:
            # extrapolate in top/bottom direction
            cols_with_nan = np.arange(arr.shape[1])[~np.all(np.isnan(arr), axis=0)]

            for col in cols_with_nan:  # FIXME iterating over columns is slow
                data = arr[:, col]
                idx_goodvals = np.argwhere(~np.isnan(data)).flatten()
                arr[:, col] = interp1d(idx_goodvals, data[idx_goodvals], **kw)(range(data.size))
        else:
            # extrapolate in left/right direction
            rows_with_nan = np.arange(arr.shape[0])[~np.all(np.isnan(arr), axis=1)]

            for row in rows_with_nan:  # FIXME iterating over rows is slow
                data = arr[row, :]
                idx_goodvals = np.argwhere(~np.isnan(data)).flatten()
                arr[row, :] = interp1d(idx_goodvals, data[idx_goodvals], **kw)(range(data.size))

        return arr

    def compute_geolayer(self) -> (np.ndarray, np.ndarray):
        """Compute pixel-wise lon/lat arrays based on RPC coefficients, corner coordinates and image dimensions.

        :return: (2D longitude array, 2D latitude array)
        """
        # transform corner coordinates of EnMAP image to UTM
        grid_utm_epsg = get_UTMEPSG_from_LonLat(*get_center_coord(self.enmapIm_cornerCoords))
        cornerCoordsUTM = np.array([transform_any_prj(4326, grid_utm_epsg, lon, lat)
                                    for lon, lat in self.enmapIm_cornerCoords])
        xmin, xmax = cornerCoordsUTM[:, 0].min(), cornerCoordsUTM[:, 0].max()
        ymin, ymax = cornerCoordsUTM[:, 1].min(), cornerCoordsUTM[:, 1].max()

        # get UTM bounding box and move it to the EnMAP grid
        xmin, ymin, xmax, ymax = move_extent_to_EnMAP_grid((xmin, ymin, xmax, ymax))

        # set up a regular grid of UTM points with a specific mesh width
        meshwidth = 300  # 10 EnMAP pixels
        y_grid_utm, x_grid_utm = np.meshgrid(np.arange(ymax, ymin - meshwidth, -meshwidth),
                                             np.arange(xmin, xmax + meshwidth, meshwidth),
                                             indexing='ij')

        # transform UTM grid to DEM projection
        x_grid_demPrj, y_grid_demPrj = (x_grid_utm, y_grid_utm) if prj_equal(grid_utm_epsg, self.dem.epsg) else \
            transform_coordArray(EPSG2WKT(grid_utm_epsg), EPSG2WKT(self.dem.epsg), x_grid_utm, y_grid_utm)

        # retrieve corresponding heights from DEM
        # -> resample DEM to EnMAP grid?
        xy_pairs_demPrj = np.vstack([x_grid_demPrj.flatten(), y_grid_demPrj.flatten()]).T
        heights = self.dem.read_pointData(xy_pairs_demPrj).flatten()

        # transform UTM points to lon/lat
        lon_grid, lat_grid = transform_coordArray(EPSG2WKT(grid_utm_epsg), EPSG2WKT(4326), x_grid_utm, y_grid_utm)
        lons = lon_grid.flatten()
        lats = lat_grid.flatten()

        # compute floating point EnMAP image coordinates for the selected UTM points
        rows, cols = self.transform_LonLatHeight_to_RowCol(lon=lons, lat=lats, height=heights)

        # interpolate lon/lats to get lon/lat coordinates integer image coordinates of EnMAP image
        rows_im, cols_im = self.enmapIm_dims_sensorgeo
        out_rows_grid, out_cols_grid = np.meshgrid(range(rows_im), range(cols_im), indexing='ij')
        out_xy_pairs = np.vstack([out_cols_grid.flatten(), out_rows_grid.flatten()]).T

        lons_interp = interpolate_griddata(np.array((cols, rows)).T, lons, out_xy_pairs,
                                           method='linear').reshape(*out_cols_grid.shape)
        lats_interp = interpolate_griddata(np.array((cols, rows)).T, lats, out_xy_pairs,
                                           method='linear').reshape(*out_rows_grid.shape)

        # lons_interp / lats_interp may contain NaN values in case xmin, ymin, xmax, ymax has been set too small
        # to account for RPC transformation errors
        # => fix that by extrapolation at NaN value positions

        lons_interp = self._fill_nans_at_corners(lons_interp, along_axis=1)  # extrapolate in left/right direction
        lats_interp = self._fill_nans_at_corners(lats_interp, along_axis=1)

        # return a geolayer in the exact dimensions like the EnMAP detector array
        return lons_interp, lats_interp

    def __call__(self, *args, **kwargs):
        return self.compute_geolayer()


global_dem_sensorgeo = None  # type: GeoArray


def mp_initializer_for_RPC_3D_Geolayer_Generator(dem_sensorgeo):
    """Declare global variables needed for RPC_3D_Geolayer_Generator._compute_geolayer_for_band()

    :param dem_sensorgeo:   DEM in sensor geometry
    """
    global global_dem_sensorgeo
    global_dem_sensorgeo = dem_sensorgeo


class RPC_3D_Geolayer_Generator(object):
    """
    Class for creating band- AND pixel-wise longitude/latitude arrays based on rational polynomial coefficients (RPC).
    """
    def __init__(self,
                 rpc_coeffs_per_band: dict,
                 dem: Union[str, GeoArray],
                 enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                 enmapIm_dims_sensorgeo: Tuple[int, int],
                 CPUs: int = None):
        """Get an instance of RPC_3D_Geolayer_Generator.

        :param rpc_coeffs_per_band:     RPC coefficients for a single EnMAP band ({'band_1': <rpc_coeffs_dict>,
                                                                                   'band_2': <rpc_coeffs_dict>,
                                                                                   ...}
        :param dem:                     digital elevation model in map geometry (file path or instance of GeoArray)
        :param enmapIm_cornerCoords:    corner coordinates as tuple of lon/lat pairs
        :param enmapIm_dims_sensorgeo:  dimensions of the EnMAP image in sensor geometry (rows, colunms)
        :param CPUs:                    number of CPUs to use
        """
        self.rpc_coeffs_per_band = OrderedDict(sorted(rpc_coeffs_per_band.items()))
        self.dem = dem
        self.enmapIm_cornerCoords = enmapIm_cornerCoords
        self.enmapIm_dims_sensorgeo = enmapIm_dims_sensorgeo
        self.CPUs = CPUs or cpu_count()

        # get validated DEM in map geometry
        # self.logger.debug('Verifying DEM...')  # TODO
        from ..dem_preprocessing import DEM_Processor
        self.dem = DEM_Processor(dem_path_geoarray=dem,
                                 enmapIm_cornerCoords=enmapIm_cornerCoords).dem
        # TODO clip DEM to needed area
        self.dem.to_mem()

    @staticmethod
    def _compute_geolayer_for_band(rpc_coeffs, enmapIm_cornerCoords, enmapIm_dims_sensorgeo, band_idx):
        lons, lats = \
            RPC_Geolayer_Generator(rpc_coeffs=rpc_coeffs,
                                   dem=global_dem_sensorgeo,
                                   enmapIm_cornerCoords=enmapIm_cornerCoords,
                                   enmapIm_dims_sensorgeo=enmapIm_dims_sensorgeo) \
            .compute_geolayer()

        return lons, lats, band_idx

    def compute_geolayer(self):
        rows, cols = self.enmapIm_dims_sensorgeo
        bands = len(self.rpc_coeffs_per_band)
        lons = np.empty((rows, cols, bands), dtype=np.float)
        lats = np.empty((rows, cols, bands), dtype=np.float)

        band_inds = list(range(len(self.rpc_coeffs_per_band)))
        rpc_coeffs_list = list(self.rpc_coeffs_per_band.values())

        if self.CPUs > 1:
            # multiprocessing
            args = [(coeffs, self.enmapIm_cornerCoords, self.enmapIm_dims_sensorgeo, idx)
                    for coeffs, idx in zip(rpc_coeffs_list, band_inds)]

            # FIXME: pickling back large lon/lat arrays to the main process may be an issue on small machines
            # NOTE: With the small test dataset pickling has only a small effect on processing time.
            with Pool(self.CPUs,
                      initializer=mp_initializer_for_RPC_3D_Geolayer_Generator,
                      initargs=(self.dem,)) as pool:
                results = pool.starmap(self._compute_geolayer_for_band, args)

            for res in results:
                band_lons, band_lats, band_idx = res
                lons[:, :, band_idx] = band_lons
                lats[:, :, band_idx] = band_lats

        else:
            # singleprocessing
            global global_dem_sensorgeo
            global_dem_sensorgeo = self.dem

            for band_idx in band_inds:
                lons[:, :, band_idx], lats[:, :, band_idx] = \
                    self._compute_geolayer_for_band(rpc_coeffs=rpc_coeffs_list[band_idx],
                                                    enmapIm_cornerCoords=self.enmapIm_cornerCoords,
                                                    enmapIm_dims_sensorgeo=self.enmapIm_dims_sensorgeo,
                                                    band_idx=band_idx)[:2]

        return lons, lats


def compute_mapCoords_within_sensorGeoDims(sensorgeoCoords_YX: List[Tuple[float, float]],
                                           rpc_coeffs: dict,
                                           dem: Union[str, GeoArray],
                                           enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                                           enmapIm_dims_sensorgeo: Tuple[int, int],
                                           ) -> List[Tuple[float, float]]:
    """Compute map coordinates for a given image coordinate pair of an EnMAP image in sensor geometry.

    :param sensorgeoCoords_YX:      list of requested sensor geometry positions [(row, column), (row, column), ...]
    :param rpc_coeffs:              RPC coefficients describing the relation between sensor and map geometry
    :param dem:                     digital elevation model in MAP geometry
    :param enmapIm_cornerCoords:    MAP coordinates of the EnMAP image
    :param enmapIm_dims_sensorgeo:  dimensions of the sensor geometry EnMAP image (rows, columns)
    :return:
    """
    # compute coordinate array
    RPCGG = RPC_Geolayer_Generator(rpc_coeffs=rpc_coeffs,
                                   dem=dem,
                                   enmapIm_cornerCoords=enmapIm_cornerCoords,  # order does not matter
                                   enmapIm_dims_sensorgeo=enmapIm_dims_sensorgeo)
    lons, lats = RPCGG.compute_geolayer()

    # extract the new corner coordinate from the coordinate arrays computed via RPCs
    rows, cols = enmapIm_dims_sensorgeo

    ul, ur, ll, lr = enmapIm_cornerCoords

    lonlats = []
    for imYX in sensorgeoCoords_YX:
        lonlat = \
            ul if imYX == (0, 0) else \
            ur if imYX == (0, cols - 1) else \
            ll if imYX == (rows - 1, 0) else \
            lr if imYX == (rows - 1, cols - 1) else \
            (lons[imYX], lats[imYX])

        lonlats.append(lonlat)

    return lonlats


def get_UTMEPSG_from_LonLat(lon: float, lat: float) -> int:
    zoneNr = int(1 + (lon + 180.0) / 6.0)
    isNorth = lat >= 0

    return int('326' + str(zoneNr)) if isNorth else int('327' + str(zoneNr))


def get_center_coord(cornerCoordsXY):
    # FIXME center coord is not equal to center of bounding box
    cornerCoordsXY = np.array(cornerCoordsXY)
    x_center = float(np.mean([cornerCoordsXY[:, 0].min(), cornerCoordsXY[:, 0].max()]))
    y_center = float(np.mean([cornerCoordsXY[:, 1].min(), cornerCoordsXY[:, 1].max()]))

    return x_center, y_center


def get_UTMEPSG_from_LonLat_cornersXY(lons: List[float], lats: List[float]):
    return get_UTMEPSG_from_LonLat(*get_center_coord(list(zip(lons, lats))))
