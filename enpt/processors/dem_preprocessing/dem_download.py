# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018–2025 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
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

from typing import Union, Tuple  # noqa: F401
import math

import numpy as np
from osgeo import gdal, osr

from py_tools_ds.geo.coord_trafo import transform_any_prj

__author__ = 'Daniel Scheffler'


class DEM_Downloader(object):
    def __init__(self,
                 enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                 progress: bool = False):
        self.enmapIm_cornerCoords = enmapIm_cornerCoords
        self.progress = progress

    @staticmethod
    def _integer_degree_range(self, min_deg: float, max_deg: float):
        return list(range(math.floor(min_deg), math.floor(max_deg) + 1))

    @staticmethod
    def _format_tile_name(self, lat_deg: int, lon_deg: int, arcsec_str: str):
        ns = f"N{abs(lat_deg):02d}_00" if lat_deg >= 0 else f"S{abs(lat_deg):02d}_00"
        ew = f"E{abs(lon_deg):03d}_00" if lon_deg >= 0 else f"W{abs(lon_deg):03d}_00"
        folder = f"Copernicus_DSM_COG_{arcsec_str}_{ns}_{ew}_DEM"
        filename = f"Copernicus_DSM_COG_{arcsec_str}_{ns}_{ew}_DEM.tif"
        return folder, filename

    @staticmethod
    def _construct_cog_urls_for_bbox(self, west, south, east, north, product="GLO-30"):
        if product == "GLO-30":
            bucket = "copernicus-dem-30m.s3.amazonaws.com"
            arcsec_str = "10"
        else:
            bucket = "copernicus-dem-90m.s3.amazonaws.com"
            arcsec_str = "30"

        urls = []
        for latd in self._integer_degree_range(south, north):
            for lond in self._integer_degree_range(west, east):
                folder, filename = self._format_tile_name(latd, lond, arcsec_str)
                url = f"https://{bucket}/{folder}/{filename}"
                urls.append(url)
        return urls

    def _download_dem_from_url(self, url, output_dir):
        pass

    def make_dem_from_header(self, hdr_file, out_prefix, product="GLO-30", out_format="ENVI", topo=False):
        width, height, transform, epsg = parse_map_info(hdr_file)

        west, north = transform[0], transform[3]
        east = west + width * transform[1]
        south = north + height * transform[5]

        wgs_west, wgs_south, wgs_east, wgs_north = transform_bbox_to_wgs84(west, south, east, north, epsg)

        urls = construct_cog_urls_for_bbox(wgs_west, wgs_south, wgs_east, wgs_north, product)
        print(f"{len(urls)} possible DEM tiles...")

        # Download and mosaic all DEM tiles
        vrt_path = "/vsimem/temp_mosaic.vrt"
        gdal.BuildVRT(vrt_path, urls)
        src_ds = gdal.Open(vrt_path)
        if src_ds is None:
            raise RuntimeError("No DEM tiles found or could not open them.")

        # Reproject to target grid
        mem_drv = gdal.GetDriverByName("MEM")
        dst_ds = mem_drv.Create("", width, height, 1, gdal.GDT_Float32)
        dst_ds.SetGeoTransform(transform)

        dst_srs = osr.SpatialReference()
        dst_srs.ImportFromEPSG(epsg)
        dst_ds.SetProjection(dst_srs.ExportToWkt())

        gdal.Warp(dst_ds, src_ds, format="MEM", resampleAlg="bilinear")

        dem_arr = dst_ds.GetRasterBand(1).ReadAsArray()
        nodata = np.finfo(np.float32).min
        dem_arr[np.isnan(dem_arr)] = nodata

        # Write output file
        driver = gdal.GetDriverByName(out_format)
        if driver is None:
            raise RuntimeError(f"Unsupported format: {out_format}")

        out_ds = driver.Create(out_prefix, width, height, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(transform)
        out_ds.SetProjection(dst_srs.ExportToWkt())
        out_ds.GetRasterBand(1).WriteArray(dem_arr)
        out_ds.GetRasterBand(1).SetNoDataValue(nodata)
        out_ds = None


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
    def _get_utm_epsg(lon, lat):
        """Return UTM EPSG code based on lon/lat center."""
        zone = int((lon + 180) / 6) + 1
        return 32600 + zone if lat >= 0 else 32700 + zone

    def _construct_tile_urls(self):
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

    def _build_vrt(self, urls):
        """Create an in-memory VRT mosaic from URLs."""
        vrt_path = "/vsimem/copernicus_mosaic.vrt"
        gdal.BuildVRT(vrt_path, urls)
        return gdal.Open(vrt_path)

    def _reproject_to_utm(self, src_ds, utm_epsg):
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
            gdal.Warp(dst_ds, src_ds, format="MEM", dstSRS=f'EPSG:{utm_epsg}',
                      resampleAlg="bilinear", outputBounds=(xmin, ymin, xmax, ymax))

            arr = dst_ds.GetRasterBand(1).ReadAsArray()
            prj_wkt = dst_ds.GetProjection()

        nodata = 0
        arr[np.isnan(arr)] = nodata

        return arr, gt, prj_wkt, nodata

    def _write_dem(self, path_out, arr, gt, prj, nodata):
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
