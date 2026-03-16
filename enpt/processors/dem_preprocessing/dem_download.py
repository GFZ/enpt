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
from geoarray import GeoArray
from osgeo import gdal

from py_tools_ds.geo.coord_trafo import transform_any_prj

__author__ = 'Daniel Scheffler'


class CopernicusDEMGenerator:
    def __init__(
        self,
        extent: tuple[float, float, float, float],
        tgt_epsg: int,
        resolution: float = 30,
        product: str = "GLO-30"
    ):
        """
        Parameters
        ----------
        extent
            (xmin, ymin, xmax, ymax) in target projection
        tgt_epsg
            EPSG code of the target projection
        resolution
            Output pixel size in target projection units
        """
        self.xmin, self.ymin, self.xmax, self.ymax = extent
        self.tgt_epsg = tgt_epsg
        self.resolution = resolution
        self.product = product

    def run(self) -> GeoArray:
        """Download, mosaic, reproject, and save DEM."""
        print(f"Target projection: EPSG:{self.tgt_epsg}")

        # convert target extent → WGS84
        wgs_bbox = self._extent_to_wgs84()

        # determine tile URLs
        urls = self._construct_tile_urls(wgs_bbox)
        print(f"{len(urls)} DEM tiles potentially required...")

        # build VRT mosaic
        src_ds = self._build_vrt(urls)
        if src_ds is None:
            raise RuntimeError("No DEM tiles found or could not open them.")

        # warp to requested grid
        arr, gt, prj, nodata = self._warp_to_target(src_ds)

        return GeoArray(arr, geotransform=gt, projection=prj, nodata=nodata)

    def _extent_to_wgs84(self) -> tuple[float, float, float, float]:
        """Transform target extent to WGS84."""
        corners = [
            (self.xmin, self.ymin),
            (self.xmin, self.ymax),
            (self.xmax, self.ymin),
            (self.xmax, self.ymax),
        ]

        ll = [transform_any_prj(self.tgt_epsg, 4326, x, y) for x, y in corners]

        xs, ys = zip(*ll)

        return min(xs), min(ys), max(xs), max(ys)

    def _construct_tile_urls(self, wgs_bbox: tuple[float, float, float, float]) -> list[str]:
        """Generate Copernicus DEM URLs covering WGS84 bbox."""
        west, south, east, north = wgs_bbox

        if self.product == "GLO-30":
            bucket, arcsec = "copernicus-dem-30m.s3.amazonaws.com", "10"
        else:
            bucket, arcsec = "copernicus-dem-90m.s3.amazonaws.com", "30"

        def deg_range(a, b):
            return range(math.floor(a), math.floor(b) + 1)

        urls = []

        for lat in deg_range(south, north):
            for lon in deg_range(west, east):
                ns = f"N{abs(lat):02d}_00" if lat >= 0 else f"S{abs(lat):02d}_00"
                ew = f"E{abs(lon):03d}_00" if lon >= 0 else f"W{abs(lon):03d}_00"

                folder = f"Copernicus_DSM_COG_{arcsec}_{ns}_{ew}_DEM"
                fname = f"{folder}.tif"

                urls.append(f"https://{bucket}/{folder}/{fname}")

        return urls

    def _build_vrt(self, urls: list[str]) -> gdal.Dataset:
        """Create in-memory VRT mosaic."""
        vrt_path = "/vsimem/copernicus_mosaic.vrt"
        gdal.BuildVRT(vrt_path, urls)
        return gdal.Open(vrt_path)

    def _warp_to_target(
            self,
            src_ds: gdal.Dataset,
            dst_nodata: int = -9999
    ) -> tuple[np.ndarray, tuple, str, float]:
        """Subset or warp DEM depending on target projection."""
        if self.tgt_epsg == 4326:
            print("Subsetting DEM (native Copernicus grid retained)...")
            rsp_alg = "nearest"
            xres = yres = None
        else:
            print("Reprojecting and resampling DEM...")
            rsp_alg = "bilinear"
            xres = yres = self.resolution

        with gdal.Warp(
                "",
                src_ds,
                format="MEM",
                dstSRS=f"EPSG:{self.tgt_epsg}",
                outputBounds=(self.xmin, self.ymin, self.xmax, self.ymax),
                resampleAlg=rsp_alg,
                xRes=xres,
                yRes=yres,
                dstNodata=dst_nodata
            ) as ds:
            arr = ds.GetRasterBand(1).ReadAsArray()
            gt = ds.GetGeoTransform()
            prj_wkt = ds.GetProjection()

        arr[np.isnan(arr)] = dst_nodata

        return arr, gt, prj_wkt, dst_nodata
