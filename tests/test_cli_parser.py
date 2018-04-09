#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_cli_parser
---------------

Tests for enpt.bin.enpt_cli.py
"""

from unittest import TestCase
import os
from runpy import run_path
from multiprocessing import cpu_count

from enpt import __path__


path_run_enpt = os.path.abspath(os.path.join(__path__[0], '..', 'bin', 'enpt_cli.py'))


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
