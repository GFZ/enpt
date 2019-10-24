#!/usr/bin/env python
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
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
test_l1b_reader
---------------

Tests for `l1b_reader` module.
"""

import unittest
import os
from os import path
import tempfile
import zipfile
import shutil

from enpt.io.reader import L1B_Reader
from enpt.model.images import EnMAPL1Product_SensorGeo
from enpt.options.config import EnPTConfig, config_for_testing, config_for_testing_dlr

__author__ = 'Daniel Scheffler'


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.config = EnPTConfig(**config_for_testing)
        self.pathList_testimages = [self.config.path_l1b_enmap_image,
                                    self.config.path_l1b_enmap_image_gapfill]
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        shutil.rmtree(self.config.output_dir)

    def test_read_and_save_inputdata(self):
        print("")
        print("################################################")
        print("#                                              #")
        print("# Test reading EnMAP Level-1B products from SG #")
        print("#                                              #")
        print("################################################")
        print("")
        print("")

        print("================================================")
        print("Unzip Test data files")
        print("================================================")
        print("")
        for l1b_file in self.pathList_testimages:
            with zipfile.ZipFile(l1b_file, "r") as zf:
                zf.extractall(self.tmpdir)
        prods = [os.path.join(self.tmpdir, os.path.basename(self.pathList_testimages[0]).split(".zip")[0]),
                 os.path.join(self.tmpdir, os.path.basename(self.pathList_testimages[1]).split(".zip")[0])]
        print("Done!")
        print("")
        print("")

        print("================================================")
        print("Create the L1B_Reader new instance")
        print("================================================")
        print("")
        rd = L1B_Reader(config=self.config)
        print("Done!")
        print("")
        print("")

        # TEST FOR ONE IMAGE ONLY
        print("=======================================================================================================")
        print("Test: Read and write ONE image only")
        print("=======================================================================================================")
        print("")
        for prod in prods:
            # for l1b_file in self.pathList_testimages:
            #     with zipfile.ZipFile(l1b_file, "r") as zf:
            #
            #         zf.extractall(self.tmpdir)
            #
            #         prod = os.path.join(self.tmpdir, os.path.basename(l1b_file).split(".zip")[0])
            print("-------------------------------------")
            print("Test with %s" % os.path.basename(prod))
            print("-------------------------------------")
            print("Tmp dir: %s" % self.tmpdir)
            print("")
            print(" * Without SNR ")
            print("")
            L1_obj = rd.read_inputdata(prod, compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNone(L1_obj.swir.detector_meta.snr)
            root_dir_written_L1_data = L1_obj.save(path.join(self.tmpdir, "no_snr"))

            # read self written L1 data
            L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

            print("")
            print(" * With SNR ")
            print("")
            L1_obj = rd.read_inputdata(prod)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
            root_dir_written_L1_data = L1_obj.save(path.join(self.tmpdir, "with_snr"))
            L1_obj = rd.read_inputdata(root_dir_written_L1_data)
            self.assertIsNotNone(L1_obj, EnMAPL1Product_SensorGeo)

        # TEST FOR ONE IMAGE ONLY
        print("=======================================================================================================")
        print("Test: read, join and write 2 images if possible!")
        print("=======================================================================================================")
        print("")

        for k_prod1, k_prod2 in ((0, 1), (1, 0)):
            for n_lines in (-1, 10, 50, 80, 100, 150):  # TODO isolate the test for different number of lines
                tempdir = tempfile.mkdtemp(dir=self.config.working_dir)
                if n_lines is -1:
                    n_lines = "all"
                print("-----------------------------------------------------------------------------------------------")
                print("Test with %s and %s, with: %s lines" % (os.path.basename(prods[k_prod1]),
                                                               os.path.basename(prods[k_prod2]),
                                                               n_lines))
                print("-----------------------------------------------------------------------------------------------")

                if n_lines == "all":
                    n_lines = None

                print("")
                print(" * Without SNR")
                print("")
                L1_obj = rd.read_inputdata(prods[k_prod1], prods[k_prod2], n_lines, compute_snr=False)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNone(L1_obj.swir.detector_meta.snr)
                root_dir_written_L1_data = L1_obj.save(path.join(tempdir, "no_snr"))
                L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                print("")
                print(" * With SNR")
                print("")
                L1_obj = rd.read_inputdata(prods[k_prod1], prods[k_prod2], n_lines)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
                root_dir_written_L1_data = L1_obj.save(path.join(tempdir, "with_snr"))
                L1_obj = rd.read_inputdata(root_dir_written_L1_data)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                print("")
                print("")
                shutil.rmtree(tempdir)

        return


class Test_L1B_Reader_DLR(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.config = EnPTConfig(**config_for_testing_dlr)
        self.pathList_testimages = [self.config.path_l1b_enmap_image,
                                    self.config.path_l1b_enmap_image_gapfill]
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        shutil.rmtree(self.config.output_dir)

    def test_read_inputdata(self):
        with zipfile.ZipFile(self.pathList_testimages[0], "r") as zf:
            zf.extractall(self.tmpdir)

        RD = L1B_Reader(config=self.config)

        L1_obj = RD.read_inputdata(self.tmpdir, compute_snr=False)
        L1_obj.save(path.join(self.config.output_dir, "no_snr"))


if __name__ == "__main__":
    unittest.main()
