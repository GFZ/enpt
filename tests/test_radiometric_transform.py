# -*- coding: utf-8 -*-

import os
import unittest
from glob import glob
import tempfile
import zipfile
from datetime import datetime

from enpt.processors.radiometric_transform import TOARad2TOARef_Transformer
from enpt.options.config import EnPTConfig

from . import config_for_testing


class Radiometric_Transformer_Tester(unittest.TestCase):

    def setUp(self):
        """Set up the needed test data"""

        self.cfg = EnPTConfig(**config_for_testing)
        self.pathList_testimages = glob(os.path.join(os.path.dirname(__file__), "data", "EnMAP_Level_1B", "*.zip"))
        self.RT = TOARad2TOARef_Transformer(config=self.cfg)

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
                    L1_obj = L1B_Reader(config=self.cfg)\
                        .read_inputdata(root_dir, observation_time=datetime(2015, 12, 7, 10))

                    # input assertions
                    self.assertIsInstance(L1_obj, EnMAPL1Product_SensorGeo)

                    # run transformation
                    output = self.RT.transform_dummy(L1_obj)  # for now only test if its runnable without error

            # output assertions
            self.assertIsInstance(output, EnMAPL1Product_SensorGeo)
