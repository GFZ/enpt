#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_l1b_reader
---------------

Tests for `l1b_reader` module.
"""

import unittest
from glob import glob
import os
from os import path
import sys
import tempfile
import zipfile
# from datetime import datetime
import shutil

from enpt.options.config import EnPTConfig
from . import config_for_testing


class Test_L1B_Reader(unittest.TestCase):
    """Tests for L1B_Reader class.."""

    def setUp(self):
        self.pathList_testimages = glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip"))
        self.config = EnPTConfig(**config_for_testing)
        self.tmpdir = tempfile.mkdtemp(dir=self.config.working_dir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)
        # shutil.rmtree(self.config.output_dir)

    def test_read_and_save_inputdata(self):
        from enpt.io.reader import L1B_Reader
        from enpt.model.images import EnMAPL1Product_SensorGeo


        print("")
        print("################################################")
        print("#                                              #")
        print("# Test reading EnMAP Level-1B products from SG #")
        print("#                                              #")
        print("################################################")
        print("")
        print("")

        print("================================================")
        print("Create the L1B_Reader new instance")
        print("================================================")
        print("")
        print("")
        rd = L1B_Reader(config=self.config)

        prods = (config_for_testing['path_l1b_enmap_image'],
                 config_for_testing['path_l1b_enmap_image_ext'],
                 config_for_testing['path_l1b_enmap_image_ext2'])

        # TEST FOR ONE IMAGE ONLY
        print("=======================================================================================================")
        print("Test: Read ONE image only")
        print("=======================================================================================================")
        print("")

        for prod in prods:
            print("-------------------------------------")
            print("Test with %s" % os.path.basename(prod))
            print("-------------------------------------")
            print("")
            print(" * Without SNR ")
            print("")
            L1_obj = rd.read_inputdata(prod, compute_snr=False)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNone(L1_obj.swir.detector_meta.snr)
        #     root_dir_written_L1_data = L1_obj.save(path.join(self.tmpdir, "no_snr"))
        #     L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
        #     self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            print("")
            print(" * With SNR ")
            print("")
            L1_obj = rd.read_inputdata(prod)
            self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
            self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
        #     root_dir_written_L1_data = L1_obj.save(path.join(self.config.output_dir, "with_snr"))
        #     L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
        #     self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
            print("")
            print("")

        print("=======================================================================================================")
        print("Test: Write ONE image only")
        print("=======================================================================================================")
        print("")
        print("-----------------------------------------------------------------------------------------------")
        print("Test with %s" % os.path.basename(prods[0]))
        print("-----------------------------------------------------------------------------------------------")
        print("")
        print(" * With SNR ")
        print("")
        L1_obj = rd.read_inputdata(prods[0])
        root_dir_written_L1_data = L1_obj.save(path.join(self.config.output_dir, "with_snr"))
        L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
        self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
        shutil.rmtree(path.join(self.config.output_dir, "with_snr"))
        print("")
        print("")


        # TEST FOR TWO IMAGES
        print("=======================================================================================================")
        print("Test: Read TWO images")
        print("=======================================================================================================")
        print("")
        for prod1, prod2 in ((0, 1),
                             (1, 0),
                             (0, 2),
                             (2, 0),
                             (1, 2),
                             (2, 1)):
            for n_lines in (10, 50, 100, 500, 1000, 1500):
                print("-----------------------------------------------------------------------------------------------")
                print("Test with %s and %s, with: %s lines" % (os.path.basename(prods[prod1]),
                                                               os.path.basename(prods[prod2]),
                                                               n_lines))
                print("-----------------------------------------------------------------------------------------------")
                print("")
                print(" * Without SNR")
                print("")
                L1_obj = rd.read_inputdata(prods[prod1], prods[prod2], n_lines, compute_snr=False)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNone(L1_obj.swir.detector_meta.snr)
                # root_dir_written_L1_data = L1_obj.save(path.join(self.tmpdir, "no_snr"))
                # L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
                # self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                print("")
                print(" * with SNR ")
                print("")
                L1_obj = rd.read_inputdata(prods[prod1], prods[prod2], n_lines)
                self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                self.assertIsNotNone(L1_obj.vnir.detector_meta.snr)
                self.assertIsNotNone(L1_obj.swir.detector_meta.snr)
                # root_dir_written_L1_data = \
                #     L1_obj.save(path.join(self.config.output_dir,
                #                           "with_snr_%s_%s_%s_lines" % (os.path.basename(prods[prod1]),
                #                                                        os.path.basename(prods[prod2]),
                #                                                        n_lines))
                #                 )
                # L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
                # self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
                print("")
                print("")

        print("=======================================================================================================")
        print("Test: Write TWO images ")
        print("=======================================================================================================")
        print("")
        print("-----------------------------------------------------------------------------------------------")
        print("Test with %s and %s, with: %s lines" % (os.path.basename(prods[0]),os.path.basename(prods[1]), 500))
        print("-----------------------------------------------------------------------------------------------")
        print("")
        print(" * With SNR ")
        print("")
        L1_obj = rd.read_inputdata(prods[0], prods[1], 500)
        root_dir_written_L1_data = L1_obj.save(path.join(self.config.output_dir, "with_snr"))
        L1_obj = rd.read_inputdata(root_dir_written_L1_data, compute_snr=False)
        self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)
        shutil.rmtree(path.join(self.config.output_dir, "with_snr"))
        print("")
        print("")


        return


if __name__ == "__main__":
    unittest.main()
