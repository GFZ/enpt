# -*- coding: utf-8 -*-

import os
import unittest

from enpt.model.images import EnMAP_L1B
from enpt.processors.radiometric_transform import TOARad2TOARef_Transformer



class Radiometric_Transformer_Tester(unittest.TestCase):

    def setUp(self):
        """Set up the needed test data"""

        path_testimage = os.path.abspath('./tests/data/Landsat8__500x500x2.bsq')
        self.EIm = EnMAP_L1B(pathImage=path_testimage)
        self.RT = TOARad2TOARef_Transformer(None, None)


    def test_transform_TOARad2TOARef(self):
        # input assertions
        self.assertIsInstance(self.EIm, EnMAP_L1B)

        # run
        output = self.RT.transform(self.EIm) # for now only test if its runnable without error

        # output assertions
        self.assertIsInstance(output, EnMAP_L1B)

