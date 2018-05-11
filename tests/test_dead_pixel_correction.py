# -*- coding: utf-8 -*-

import unittest

import numpy as np
from geoarray import GeoArray

from enpt.processors.dead_pixel_correction import Dead_Pixel_Corrector


class Test_Dead_Pixel_Corrector(unittest.TestCase):

    def setUp(self):
        """Set up the needed test data"""

        # self.cfg = EnPTConfig(**config_for_testing)
        # self.RT = Dead_Pixel_Corrector(config=self.cfg)
        self.DPC = Dead_Pixel_Corrector(algorithm='spectral', interp='linear')

        # create test data
        self.im = np.random.randint(0, 10, (50, 1000, 88), np.int16)  # VNIR size
        self.deadpixelmap = np.zeros((self.im.shape[2], self.im.shape[1]))

        for band, column in \
            [[0, 2],  # first band
             [1, 2],  # first band, same column
             [50, 4],  # 2 adjacent bands
             [51, 4],  # 2 adjacent bands
             [60, 20],  # single dead column
             [87, 50],  # second last band, same column
             [86, 50],  # last band, same column
             [87, 2]]:  # single dead column, last band
            self.im[:, column, band] = 0
            self.deadpixelmap[band, column] = 1

    def test_correct(self):
        corrected = self.DPC.correct(self.im, self.deadpixelmap)

        # output assertions
        self.assertIsInstance(corrected, (GeoArray, np.ndarray))
        # self.assertNotEqual(np.mean(self.im[:, 2, 0]), 0)  # first band
        # self.assertNotEqual(np.mean(self.im[:, 2, 1]), 0)  # first band, same column
        self.assertNotEqual(np.mean(self.im[:, 4, 50]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 4, 10]), 0)  # 2 adjacent bands
        self.assertNotEqual(np.mean(self.im[:, 20, 60]), 0)  # single dead column
        # self.assertNotEqual(np.mean(self.im[:, 50, 86]), 0)  # second last band, same column
        # self.assertNotEqual(np.mean(self.im[:, 50, 87]), 0)  # last band, same column
        # self.assertNotEqual(np.mean(self.im[:, 2, 87]), 0)  # single dead column, last band
