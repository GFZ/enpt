#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_controller
---------------

Tests for `execution.controller` module.
"""

from unittest import TestCase, main
import shutil

from enpt.execution.controller import EnPT_Controller
from enpt.options.config import config_for_testing


class Test_EnPT_Controller(TestCase):
    def setUp(self):
        self.CTR = EnPT_Controller(**config_for_testing)

    def tearDown(self):
        # NOTE: ignore_errors deletes the folder, regardless of whether it contains read-only files
        shutil.rmtree(self.CTR.cfg.output_dir, ignore_errors=True)

    def test_run_all_processors(self):
        self.CTR.run_all_processors()


if __name__ == '__main__':
    main()
