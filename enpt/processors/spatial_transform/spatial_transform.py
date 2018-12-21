# -*- coding: utf-8 -*-
"""EnPT module 'spatial transform', containing everything related to spatial transformations."""

import numpy as np
from typing import Union, Tuple, List  # noqa: F401
from scipy.interpolate import griddata as interpolate_griddata, interp1d
from geoarray import GeoArray

from py_tools_ds.geo.raster.reproject import SensorMapGeometryTransformer, AreaDefinition
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict, prj_equal, EPSG2WKT
from py_tools_ds.geo.coord_grid import find_nearest
from py_tools_ds.geo.coord_trafo import transform_any_prj, transform_coordArray

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
        for llh, name in zip([lon_norm, lat_norm, height_norm], ['longitudes', 'latitudes', 'heights']):
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
        :return:
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
        :return:
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

        pos_nan = np.argwhere(np.isnan(arr))
        rows_nan, cols_nan = pos_nan[:, 0], pos_nan[:, 1]
        kw = dict(kind='linear', fill_value='extrapolate')

        if along_axis == 0:
            # extrapolate in top/bottom direction
            for col in cols_nan:
                data = arr[:, col]
                idx_goodvals = np.argwhere(~np.isnan(data)).flatten()
                arr[:, col] = interp1d(idx_goodvals, data[idx_goodvals], **kw)(range(data.size))
        else:
            # extrapolate in left/right direction
            for row in rows_nan:
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
