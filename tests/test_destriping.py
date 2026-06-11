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
test_destriping
---------------

Tests for `processors.destriping.destriping` module.
"""

from unittest import TestCase
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import pytest
from geoarray import GeoArray

from enpt.processors.destriping.destriping import Destriper
from enpt.options.config import config_for_testing, EnPTConfig
from enpt.io.reader import L1B_Reader

__author__ = 'Daniel Scheffler'


class Test_Destriping(TestCase):
    def test_destriping(self):
        cfg = EnPTConfig(**config_for_testing)

        with ZipFile(cfg.path_l1b_enmap_image, "r") as zf, \
             TemporaryDirectory(cfg.working_dir) as td:
            zf.extractall(td)
            L1_obj = L1B_Reader(config=cfg).read_inputdata(
                root_dir_main=td,
                compute_snr=False)
            swir_sub = L1_obj.swir.data.get_subset(zslice=slice(38, 40))

        # TODO: prepare a subset that has horizontal stripes

        # swir_sub = GeoArray('/home/gfz-fe/scheffler/temp/EnPT/destriping/ENMAP_DT0000009666_SWIR_sub.bsq')

        a = 1
        dst, diff = (
            Destriper().destripe(
                array=swir_sub[:],
                sensor='swir',
                high_freq=True,
                low_freq=False,
                spatial_domain=True,
                spectral_domain=False,
                mode='stripes',
                along_track_direction=True
            ))
        dst_gA = GeoArray(dst)
        dst_gA.show(band=1)
        diff_gA = GeoArray(diff)
        diff_gA.show(band=1)
        a = 1


if __name__ == '__main__':
    pytest.main()
