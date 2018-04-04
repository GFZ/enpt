# -*- coding: utf-8 -*-

"""
test_srf
--------

Tests for `model.srf` module.
"""

from unittest import TestCase

from enpt.model.srf import SRF


class Test_SRF(TestCase):
    def setUp(self):
        pass

    def test_from_cwl_fwhm(self):
        srf = SRF.from_cwl_fwhm(cwls=[800, 1000], fwhms=[10, 20])
        self.assertIsInstance(srf, SRF)
