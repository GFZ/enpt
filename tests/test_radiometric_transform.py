#!/usr/bin/env python
#  -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
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

import os
from unittest import TestCase
import tempfile
import zipfile

from enpt.processors.radiometric_transform import Radiometric_Transformer
from enpt.options.config import EnPTConfig, config_for_testing

__author__ = 'Daniel Scheffler'


class Test_Radiometric_Transformer(TestCase):

    def setUp(self):
        """Set up the needed test data"""
        self.cfg = EnPTConfig(**config_for_testing)
        self.pathList_testimages = [self.cfg.path_l1b_enmap_image,
                                    self.cfg.path_l1b_enmap_image_gapfill]
        self.RT = Radiometric_Transformer(config=self.cfg)

    def test_transform_TOARad2TOARef(self):
        from enpt.io.reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_SensorGeo

        for l1b_file in self.pathList_testimages:
            with tempfile.TemporaryDirectory() as tmpdir:
                print("Tmp dir: %s" % tmpdir)
                with zipfile.ZipFile(l1b_file, "r") as zf:
                    zf.extractall(tmpdir)

                    root_dir = os.path.join(tmpdir, os.listdir(tmpdir)[0])

                    # create EnPT Level 1 image
                    L1_obj = L1B_Reader(config=self.cfg).read_inputdata(root_dir)

                    # input assertions
                    assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)

                    # run transformation
                    L1_obj = self.RT.transform_TOARad2TOARef(L1_obj)  # for now only test if its runnable without error

            # output assertions
            assert isinstance(L1_obj, EnMAPL1Product_SensorGeo)
            assert L1_obj.vnir.detector_meta.unitcode == 'TOARef'
            assert L1_obj.swir.detector_meta.unitcode == 'TOARef'


if __name__ == '__main__':
    import pytest
    pytest.main()
