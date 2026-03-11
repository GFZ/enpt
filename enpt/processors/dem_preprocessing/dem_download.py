# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018–2026 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz.de)
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

"""EnPT downloader module for digital elevation models."""

import math

import numpy as np
from osgeo import gdal

from py_tools_ds.geo.coord_trafo import transform_any_prj

__author__ = 'Daniel Scheffler'


class CopernicusDEMGenerator:
    def __init__(
            self,
            west: float,
            south: float,
            east: float,
            north: float,
            resolution: float = 30,
            product: str = "GLO-30",
            out_format: str = "ENVI"
    ):
        self.west = west
        self.south = south
        self.east = east
        self.north = north
        self.resolution = resolution
        self.product = product
        self.out_format = out_format

    def run(self, path_out: str):
        """Main entry point: download, mosaic, reproject, and save DEM."""
        # Determine the correct UTM zone and EPSG code
        utm_epsg = self._get_utm_epsg((self.west + self.east) / 2,
                                      (self.south + self.north) / 2)
        print(f"Target projection: EPSG:{utm_epsg} (UTM zone)")

        # Determine tile URLs (in WGS84 degrees)
        urls = self._construct_tile_urls()
        print(f"{len(urls)} DEM tiles potentially required...")

        # Build VRT mosaic
        src_ds = self._build_vrt(urls)
        if src_ds is None:
            raise RuntimeError("No DEM tiles found or could not open them.")

        # Reproject to UTM
        arr, gt, prj, nodata = self._reproject_to_utm(src_ds, utm_epsg)

        # Write output
        self._write_dem(path_out, arr, gt, prj, nodata)
        print("Done.")

    @staticmethod
    def _get_utm_epsg(lon: float, lat: float) -> int:
        """Return UTM EPSG code based on lon/lat center."""
        zone = int((lon + 180) / 6) + 1
        return 32600 + zone if lat >= 0 else 32700 + zone

    def _construct_tile_urls(self) -> list[str]:
        """Generate Copernicus DEM COG URLs covering bbox."""
        if self.product == "GLO-30":
            bucket, arcsec = "copernicus-dem-30m.s3.amazonaws.com", "10"
        else:
            bucket, arcsec = "copernicus-dem-90m.s3.amazonaws.com", "30"

        def deg_range(a, b):
            return range(math.floor(a), math.floor(b) + 1)

        urls = []
        for lat in deg_range(self.south, self.north):
            for lon in deg_range(self.west, self.east):
                ns = f"N{abs(lat):02d}_00" if lat >= 0 else f"S{abs(lat):02d}_00"
                ew = f"E{abs(lon):03d}_00" if lon >= 0 else f"W{abs(lon):03d}_00"
                folder = f"Copernicus_DSM_COG_{arcsec}_{ns}_{ew}_DEM"
                fname = f"{folder}.tif"
                urls.append(f"https://{bucket}/{folder}/{fname}")
        return urls

    def _build_vrt(self, urls: list[str]):
        """Create an in-memory VRT mosaic from URLs."""
        vrt_path = "/vsimem/copernicus_mosaic.vrt"
        gdal.BuildVRT(vrt_path, urls)
        return gdal.Open(vrt_path)

    def _reproject_to_utm(
            self, src_ds: gdal.Dataset,
            utm_epsg: int
    ) -> tuple[np.ndarray, tuple, str, float]:
        """Warp DEM mosaic into UTM projection with fixed 30 m pixels."""
        UL_UR_LL_LR_ll = (
            (self.west, self.north),
            (self.east, self.north),
            (self.west, self.south),
            (self.east, self.south),
        )
        UL_UR_LL_LR_prj = [transform_any_prj(4326, utm_epsg, x, y) for x, y in UL_UR_LL_LR_ll]
        xs, ys = zip(*UL_UR_LL_LR_prj)
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)

        width = int((xmax - xmin) / self.resolution)
        height = int((ymax - ymin) / self.resolution)
        gt = (xmin, self.resolution, 0, ymax, 0, -self.resolution)

        with gdal.GetDriverByName("MEM").Create("", width, height, 1, gdal.GDT_Float32) as dst_ds:
            dst_ds.SetGeoTransform(gt)
            dst_ds.SetProjection(f'EPSG:{utm_epsg}')

            print("Reprojecting and resampling DEM...")
            gdal.Warp(dst_ds, src_ds, resampleAlg="bilinear", dstNodata=-9999)

            arr = dst_ds.GetRasterBand(1).ReadAsArray()
            prj_wkt = dst_ds.GetProjection()

        nodata = -9999
        arr[np.isnan(arr)] = nodata

        return arr, gt, prj_wkt, nodata

    def _write_dem(
            self, path_out: str,
            arr: np.ndarray,
            gt: tuple,
            prj: str,
            nodata: float
    ):
        """Write DEM array to disk."""
        driver = gdal.GetDriverByName(self.out_format)
        if driver is None:
            raise RuntimeError(f"Unsupported format: {self.out_format}")
        rows, cols = arr.shape

        with driver.Create(path_out, cols, rows, 1, gdal.GDT_Float32) as ds:
            ds.SetGeoTransform(gt)
            ds.SetProjection(prj)
            band = ds.GetRasterBand(1)
            band.WriteArray(arr)
            band.SetNoDataValue(nodata)
            ds.FlushCache()

        print(f"DEM saved: {path_out}")
