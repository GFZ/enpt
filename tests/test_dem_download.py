#!/usr/bin/env python
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

"""
test_dem_preprocessing
----------------------

Tests for `processors.dem_preprocessing.dem_download` module.
"""

from unittest import TestCase

import pytest
from geoarray import GeoArray

from enpt.processors.dem_preprocessing import CopernicusDEMGenerator

__author__ = 'Daniel Scheffler'


class Test_CopernicusDEMGenerator(TestCase):
    @staticmethod
    def _validate(dem: GeoArray, tgt_epsg: int):
        assert dem.size
        assert dem.is_map_geo
        assert dem.is_inmem
        assert dem.epsg == tgt_epsg
        assert dem[:].std() > 0

    def test_run_lonlat_glo30(self):
        dem = CopernicusDEMGenerator(
            extent=(11.0, 47.0, 11.3, 47.2),
            tgt_epsg=4326,
            product="GLO-30",
        ).run()
        self._validate(dem, 4326)

    def test_run_utm_glo30(self):
        dem = CopernicusDEMGenerator(
            extent=(636690.0, 4940340.0, 666600.0, 4950210.0),
            tgt_epsg=32630,
            product="GLO-30",
        ).run()
        self._validate(dem, 32630)

    def test_run_lonlat_glo90(self):
        dem = CopernicusDEMGenerator(
            extent=(11.0, 47.0, 11.3, 47.2),
            tgt_epsg=4326,
            product="GLO-90",
        ).run()
        self._validate(dem, 4326)

    def test_run_utm_glo90(self):
        dem = CopernicusDEMGenerator(
            extent=(636690.0, 4940340.0, 666600.0, 4950210.0),
            tgt_epsg=32630,
            product="GLO-90",
        ).run()
        self._validate(dem, 32630)

    def test_wrong_product(self):
        with pytest.raises(ValueError):
            CopernicusDEMGenerator(
                extent=(636690.0, 4940340.0, 666600.0, 4950210.0),
                tgt_epsg=32630,
                product="GLO-120")


if __name__ == '__main__':
    pytest.main()
