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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
test_cli_parser
---------------

Tests for enpt.bin.enpt_cli.py
"""

from unittest import TestCase
from argparse import ArgumentError
import os
from runpy import run_path
from multiprocessing import cpu_count

import enpt

__author__ = 'Daniel Scheffler'


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

    def test_param_list(self):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--target_coord_grid', '0', '30', '0', '30'])
        self.get_config(parsed_args)  # we don't check the result here as EnPT_Config generates a dict from it

        try:
            self.parser_run.parse_args(self.baseargs + ['--target_coord_grid', '0', '30'])
        except SystemExit as e:
            assert isinstance(e.__context__, ArgumentError)
        else:
            raise ValueError("Exception not raised")

    def test_param_boolean(self):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', 'True'])
        config = self.get_config(parsed_args)
        self.assertIsInstance(config.enable_ac, bool)
        self.assertEqual(config.enable_ac, True)

        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', 'false'])
        config = self.get_config(parsed_args)
        self.assertIsInstance(config.enable_ac, bool)
        self.assertEqual(config.enable_ac, False)

        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', '0'])
        config = self.get_config(parsed_args)
        self.assertIsInstance(config.enable_ac, bool)
        self.assertEqual(config.enable_ac, False)

        try:
            self.parser_run.parse_args(self.baseargs + ['--enable_ac', 'treu'])
        except SystemExit as e:
            assert isinstance(e.__context__, ArgumentError)
        else:
            raise ValueError("Exception not raised")

    def test_json_opts(self):
        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": 10}}'])
        config = self.get_config(parsed_args)
        self.assertEqual(config.CPUs, 10)

        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": "None"}}'])
        config = self.get_config(parsed_args)
        self.assertEqual(config.CPUs, cpu_count())


if __name__ == '__main__':
    import nose2
    nose2.main()
