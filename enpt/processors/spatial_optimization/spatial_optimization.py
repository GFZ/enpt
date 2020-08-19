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

from arosics import COREG_LOCAL
from geoarray import GeoArray

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

        # transform from sensor to map geometry to make ithe band usable for tie point detection
        enmap_ImageL1.logger.info('Temporarily transforming EnMAP band %d to map geometry for co-registration.'
                                  % (bandidx + 1))
        GT = Geometry_Transformer(lons=enmap_ImageL1.meta.vnir.lons[:, :, bandidx],
                                  lats=enmap_ImageL1.meta.vnir.lats[:, :, bandidx],
                                  fill_value=0,
                                  resamp_alg='gauss',
                                  radius_of_influence=30,
                                  nprocs=self.cfg.CPUs)

        ref_gA = GeoArray(self.cfg.path_reference_image)
        ref_ULx, ref_ULy = ref_gA.box.boxMapXY[0]

        enmap_band_mapgeo = \
            GeoArray(*GT.to_map_geometry(enmap_band_sensorgeo,
                                         tgt_prj=ref_gA.prj,  # TODO correct?
                                         tgt_coordgrid=((ref_ULx, ref_ULx + ref_gA.xgsd),
                                                        (ref_ULy, ref_ULy + ref_gA.ygsd))
                                         ),
                     nodata=0)

        return enmap_band_mapgeo

    def _compute_tie_points(self,
                            enmap_ImageL1: EnMAPL1Product_SensorGeo):
        enmap_band_mapgeo = self._get_enmap_band_for_matching(enmap_ImageL1)

        CRL = COREG_LOCAL(self.cfg.path_reference_image,
                          enmap_band_mapgeo,
                          grid_res=50)
        TPG = CRL.tiepoint_grid
        CRL.view_CoRegPoints(shapes2plot='vectors', hide_filtered=False)

        # filter out tie points over water
        pass  # TODO
        a = 1  # FIXME

    def optimize_geolayer(self,
                          enmap_ImageL1: EnMAPL1Product_SensorGeo):
        if self.cfg.enable_absolute_coreg:
            self._compute_tie_points(enmap_ImageL1)
