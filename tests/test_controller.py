# -*- coding: utf-8 -*-

"""
test_controller
---------------

Tests for `execution.controller` module.
"""

from unittest import TestCase
import shutil

from enpt.execution.controller import EnPT_Controller

from . import config_for_testing


class Test_EnPT_Controller(TestCase):
    def setUp(self):
        self.CTR = EnPT_Controller(**config_for_testing)

    def tearDown(self):
        shutil.rmtree(self.CTR.cfg.output_dir)

    def test_run_all_processors(self):
        self.CTR.run_all_processors()
