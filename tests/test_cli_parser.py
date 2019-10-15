#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tools - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Daniel Scheffler (GFZ Potsdam, daniel.scheffler@gfz-potsdam.de)
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
test_cli_parser
---------------

Tests for enpt.bin.enpt_cli.py
"""

from unittest import TestCase
import os
from runpy import run_path
from multiprocessing import cpu_count

import enpt


path_run_enpt = os.path.abspath(os.path.join(enpt.__path__[0], '..', 'bin', 'enpt_cli.py'))


class Test_CLIParser(TestCase):
    baseargs = []

    def setUp(self):
        self.parser_run = run_path(path_run_enpt)['get_enpt_argparser']()
        self.get_config = run_path(path_run_enpt)['get_config']

    def test_param_acceptance(self):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--CPUs', '10'])
        config = self.get_config(parsed_args)
        self.assertEqual(config.CPUs, 10)

        # test if parameter fallbacks are working ('CPUs' has a fallback)
        parsed_args = self.parser_run.parse_args(self.baseargs)
        config = self.get_config(parsed_args)

        self.assertNotIsInstance(config.CPUs, str)
        self.assertEqual(config.CPUs, cpu_count())

    def test_json_opts(self):
        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": 10}}'])
        config = self.get_config(parsed_args)
        self.assertEqual(config.CPUs, 10)

        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": "None"}}'])
        config = self.get_config(parsed_args)
        self.assertEqual(config.CPUs, cpu_count())
