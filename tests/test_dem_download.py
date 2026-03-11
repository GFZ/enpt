#!/usr/bin/env python
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

"""
test_dem_preprocessing
----------------------

Tests for `processors.dem_preprocessing.dem_download` module.
"""

from unittest import TestCase
from tempfile import TemporaryDirectory
import os

import pytest
from geoarray import GeoArray

from enpt.processors.dem_preprocessing import CopernicusDEMGenerator

__author__ = 'Daniel Scheffler'


class Test_CopernicusDEMGenerator(TestCase):
    def test_CopernicusDEMGenerator(self):
        with TemporaryDirectory() as td:
            demgen = CopernicusDEMGenerator(
                west=11.0, south=47.0, east=11.3, north=47.2,
                product="GLO-30", out_format="GTiff"
            )
            demgen.run(os.path.join(td, "output_dem.tif"))
            dem = GeoArray(os.path.join(td, "output_dem.tif"))[:]

            assert os.path.exists(os.path.join(td, "output_dem.tif"))
            assert dem.std() > 0

if __name__ == '__main__':
    pytest.main()
