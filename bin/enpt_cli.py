#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""EnPT console argument parser."""

import argparse

from enpt import __version__
from enpt.options.config import EnPTConfig
from enpt.execution.controller import EnPT_controller


def get_enpt_argparser():
    """Return argument parser for enpt_cli.py program."""

    ##########################################################
    # CONFIGURE MAIN PARSER FOR THE EnPT PREPROCESSING CHAIN #
    ##########################################################

    parser = argparse.ArgumentParser(
        prog='enpt_cli.py',
        description='=' * 70 + '\n' + 'EnMAP Processing Tools console argument parser. ',
        epilog="use '>>> enpt_cli.py -h' for detailed documentation and usage hints.")

    add = parser.add_argument
    add('--version', action='version', version=__version__)

    add('-jc', '--json_config', nargs='?', type=str,
        help='file path of a JSON file containing options. See here for an example: '
             'https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/'
             'EnPT/blob/master/enpt/options/options_default.json')
    add('--CPUs', type=int, default=None,
        help='number of CPU cores to be used for processing (default: "None" -> use all available')

    # link parser to run function
    parser.set_defaults(func=run_job)

    return parser


def parsedArgs_to_user_opts(cli_args):
    # type: (argparse.Namespace) -> dict
    """Convert argparse Namespace object to dictionary of explicitly given parameters.

    NOTE:   All options that have not been given explicitly (None values) are removed. Reason: EnPTConfig prefers
            directly passed arguments against those that are passed withi a JSON config file.
            So, e.g., if CPUs=None (default), the 'CPUs' parameter given within a JSON config file would be overridden.

            => only override JSON configuration if parameters are explicitly given (e.g., CPUs is set to 10)
            => if json_opts are given: default options are overridden with the options in the JSON config.

    :param cli_args:    options as parsed by the argparse.ArgumentParser
    """

    # convert argparse Namespace object to dictionary
    opts = {k: v for k, v in vars(cli_args).items() if not k.startswith('_') and k != 'func'}

    # remove those options that have not been given explicitly (None values)
    user_opts = dict()
    for k, v in opts.items():
        # values are None if they are not given by the user -> don't pass to set_config
        if v is None:
            continue
        else:
            user_opts.update({k: v})

    return user_opts


def get_config(cli_args: argparse.Namespace):
    return EnPTConfig(**parsedArgs_to_user_opts(cli_args))


def run_job(config: EnPTConfig):
    EnPT_controller(config)


if __name__ == '__main__':
    parsed_args = get_enpt_argparser().parse_args()
    parsed_args.func(parsed_args)
