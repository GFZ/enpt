#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
test_cli_parser
---------------

Tests for enpt.bin.cli.py
"""

from unittest import TestCase
import os
from runpy import run_path
from multiprocessing import cpu_count
from io import StringIO
from unittest.mock import patch

import pytest
import enpt

__author__ = 'Daniel Scheffler'


path_run_enpt = os.path.abspath(os.path.join(enpt.__path__[0], 'cli.py'))


class Test_CLIParser(TestCase):
    baseargs = []

    def setUp(self):
        self.parser_run = run_path(path_run_enpt)['get_enpt_argparser']()
        self.get_config = run_path(path_run_enpt)['get_config']

    def test_param_acceptance(self):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--CPUs', '10'])
        config = self.get_config(parsed_args)
        assert config.CPUs == 10

        # test if parameter fallbacks are working ('CPUs' has a fallback)
        parsed_args = self.parser_run.parse_args(self.baseargs)
        config = self.get_config(parsed_args)

        assert not isinstance(config.CPUs, str)
        assert config.CPUs == cpu_count()

    @patch('sys.stderr', new_callable=StringIO)  # catch STDERR so it does not pollute the test output
    def test_param_list(self, mock_stderr):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--target_coord_grid', '0', '30', '0', '30'])
        self.get_config(parsed_args)  # we don't check the result here as EnPT_Config generates a dict from it

        with (pytest.raises(SystemExit)):
            self.parser_run.parse_args(self.baseargs + ['--target_coord_grid', '0', '30'])

        assert '-tgtgrid/--target_coord_grid: expected 4 arguments' in mock_stderr.getvalue()

    @patch('sys.stderr', new_callable=StringIO)  # catch STDERR so it does not pollute the test output
    def test_param_boolean(self, mock_stderr):
        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', 'True'])
        config = self.get_config(parsed_args)
        assert isinstance(config.enable_ac, bool)
        assert config.enable_ac is True

        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', 'false'])
        config = self.get_config(parsed_args)
        assert isinstance(config.enable_ac, bool)
        assert config.enable_ac is False

        parsed_args = self.parser_run.parse_args(self.baseargs +
                                                 ['--enable_ac', '0'])
        config = self.get_config(parsed_args)
        assert isinstance(config.enable_ac, bool)
        assert config.enable_ac is False

        with (pytest.raises(SystemExit)):
            self.parser_run.parse_args(self.baseargs + ['--enable_ac', 'treu'])

        assert '--enable_ac: Boolean value expected.' in mock_stderr.getvalue()

    def test_json_opts(self):
        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": 10}}'])
        config = self.get_config(parsed_args)
        assert config.CPUs == 10

        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['--json_config', '{"general_opts": {"CPUs": "None"}}'])
        config = self.get_config(parsed_args)
        assert config.CPUs == cpu_count()

    def test_tgtprj(self):
        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['-tgtprj', 'Geographic'])
        config = self.get_config(parsed_args)
        CTR = enpt.EnPT_Controller(config)
        assert CTR.cfg.target_projection_type == 'Geographic'
        assert CTR.cfg.target_epsg == 4326

        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['-tgtprj', 'UTM'])
        config = self.get_config(parsed_args)
        CTR = enpt.EnPT_Controller(config)
        assert CTR.cfg.target_projection_type == 'UTM'

    def test_tgtepsg(self):
        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['-tgtepsg', '32617'])
        config = self.get_config(parsed_args)
        CTR = enpt.EnPT_Controller(config)
        assert CTR.cfg.target_projection_type == 'UTM'
        assert CTR.cfg.target_epsg == 32617

        parsed_args = self.parser_run.parse_args(
            self.baseargs + ['-tgtepsg', '4326'])
        config = self.get_config(parsed_args)
        CTR = enpt.EnPT_Controller(config)
        assert CTR.cfg.target_projection_type == 'Geographic'
        assert CTR.cfg.target_epsg == 4326


if __name__ == '__main__':
    pytest.main()
