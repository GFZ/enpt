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

"""EnPT pre-processing module for digital elevation models."""

from typing import Union, Tuple  # noqa: F401
from multiprocessing import cpu_count
import numpy as np

from geoarray import GeoArray
from py_tools_ds.geo.projection import get_proj4info, proj4_to_dict
from py_tools_ds.geo.coord_trafo import reproject_shapelyGeometry, transform_any_prj
from py_tools_ds.geo.vector.topology import get_footprint_polygon, get_overlap_polygon

from ..spatial_transform import Geometry_Transformer

__author__ = 'Daniel Scheffler'


class DEM_Processor(object):
    def __init__(self, dem_path_geoarray: Union[str, GeoArray],
                 enmapIm_cornerCoords: Tuple[Tuple[float, float]],
                 CPUs: int = None):
        self.dem = GeoArray(dem_path_geoarray)
        self.enmapIm_cornerCoords = enmapIm_cornerCoords
        self.CPUs = CPUs or cpu_count()

        self._set_nodata_if_not_provided()
        self._validate_input()

    def _validate_input(self):
        # check geocoding of DEM
        if not self.dem.is_map_geo:
            raise ValueError((self.dem.gt, self.dem.prj),
                             'The provided digital elevation model has no valid geo-coding or projection.')

        # check if provided as WGS-84
        proj4dict = proj4_to_dict(get_proj4info(proj=self.dem.prj))
        if 'datum' not in proj4dict or proj4dict['datum'] != 'WGS84':
            raise ValueError(proj4dict,
                             "The digital elevation model must be provided with 'WGS84' as geographic datum.")

        # check overlap
        dem_ll_mapPoly = reproject_shapelyGeometry(self.dem.footprint_poly, prj_src=self.dem.epsg, prj_tgt=4326)
        enmapIm_ll_mapPoly = get_footprint_polygon(self.enmapIm_cornerCoords, fix_invalid=True)
        overlap_dict = get_overlap_polygon(dem_ll_mapPoly, enmapIm_ll_mapPoly)
        overlap_perc = round(overlap_dict['overlap percentage'], 4)

        if overlap_perc < 100:
            # compute minimal extent in user provided projection
            cornersXY = np.array([transform_any_prj(4326, self.dem.epsg, x, y) for x, y in self.enmapIm_cornerCoords])
            xmin, xmax = cornersXY[:, 0].min(), cornersXY[:, 0].max()
            ymin, ymax = cornersXY[:, 1].min(), cornersXY[:, 1].max()

            raise ValueError('The provided digital elevation model covers %.1f percent of the EnMAP image. It must '
                             'cover it completely. The minimal needed extent in the provided projection is: \n'
                             'xmin: %s, xmax: %s, ymin: %s, ymax: %s.' % (overlap_perc, xmin, xmax, ymin, ymax))

    def _set_nodata_if_not_provided(self):
        # noinspection PyProtectedMember
        if self.dem._nodata is None:
            self.dem.nodata = -9999

    def fill_gaps(self):
        pass

    def compute_slopes(self):
        # compute on map geometry (as provided)
        pass

    def compute_aspect(self):
        # compute on map geometry (as provided)
        pass

    def to_sensor_geometry(self,
                           lons: np.ndarray,
                           lats: np.ndarray):
        GT = Geometry_Transformer(lons=lons, lats=lats, nprocs=self.CPUs)
        data_sensorgeo = GT.to_sensor_geometry(self.dem)

        return GeoArray(data_sensorgeo)
