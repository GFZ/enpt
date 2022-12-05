#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2022 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
test_config
-----------

Tests for `options.config` module.
"""
import os
from json import \
    dumps, \
    JSONDecodeError

from unittest import TestCase

from enpt.options.config import \
    get_options, \
    path_options_default, \
    EnPTConfig, \
    EnPTValidator, \
    enpt_schema_config_output

__author__ = 'Daniel Scheffler'


class Test_get_options(TestCase):
    def test_target_is_file_no_validation(self):
        opts_dict = get_options(os.path.join(path_options_default), validation=False)
        self.assertIsInstance(opts_dict, dict)

    def test_target_is_file_validation(self):
        opts_dict = get_options(os.path.join(path_options_default))
        self.assertIsInstance(opts_dict, dict)


class Test_EnPTConfig(TestCase):
    def test_plain_args(self):
        cfg = EnPTConfig(CPUs=10)
        self.assertIsInstance(cfg, EnPTConfig)
        self.assertTrue(cfg.CPUs == 10)

    def test_jsonconfig_str_allfine(self):
        cfg = '{"a": 1 /*comment*/, "b":2}'
        cfg = EnPTConfig(json_config=cfg)
        self.assertIsInstance(cfg, EnPTConfig)

    def test_jsonconfig_str_nojson(self):
        cfg = 'dict(a=1 /*comment*/, b=2)'
        with self.assertRaises(ValueError):
            EnPTConfig(json_config=cfg)

    def test_jsonconfig_str_badcomment(self):
        cfg = '{"a": 1 /comment*/, "b":2}'
        with self.assertWarns(UserWarning), self.assertRaises(JSONDecodeError):
            EnPTConfig(json_config=cfg)

    def test_jsonconfig_str_undecodable_val(self):
        cfg = '{"a": None /comment*/, "b":2}'
        with self.assertWarns(UserWarning), self.assertRaises(JSONDecodeError):
            EnPTConfig(json_config=cfg)

    def test_jsonconfig_str_schema_violation(self):
        cfg = '{"general_opts": {"CPUs": "badvalue"}}'
        with self.assertRaises(ValueError):
            EnPTConfig(json_config=cfg)

    def test_jsonconfig_file(self):
        cfg = os.path.join(path_options_default)
        cfg = EnPTConfig(json_config=cfg)
        self.assertIsInstance(cfg, EnPTConfig)

    def test_jsonconfig_param_acceptance(self):
        cfg = EnPTConfig(json_config='{"general_opts": {"CPUs": 10}}')
        self.assertIsInstance(cfg, EnPTConfig)
        self.assertTrue(cfg.CPUs == 10)

    def test_to_jsonable_dict(self):
        cfg = EnPTConfig()
        jsonable_dict = cfg.to_jsonable_dict()
        self.assertIsInstance(cfg.to_jsonable_dict(), dict)

        # test if dict is jsonable
        dumps(jsonable_dict)

    def test_to_dict_validity(self):
        cfg = EnPTConfig()
        params = cfg.to_dict()
        self.assertIsInstance(cfg.to_jsonable_dict(), dict)

        # check validity
        EnPTValidator(allow_unknown=True, schema=enpt_schema_config_output).validate(params)

    def test_invalid_filepath(self):
        with self.assertRaises(FileNotFoundError):
            EnPTConfig(path_l1b_enmap_image='/path/to/not/existing/image.tif')


if __name__ == '__main__':
    import pytest
    pytest.main()
