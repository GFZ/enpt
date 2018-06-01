#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_config
-----------

Tests for `options.config` module.
"""
import os
from json import \
    dumps, \
    JSONDecodeError

from unittest import TestCase, main

from enpt.options.config import \
    get_options, \
    path_options_default, \
    EnPTConfig, \
    EnPTValidator, \
    enpt_schema_config_output


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


if __name__ == '__main__':
    main()
