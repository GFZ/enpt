#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""EnPT console argument parser."""

import argparse

from enpt import __version__
from enpt.options.config import EnPTConfig
from enpt.execution.controller import EnPT_Controller


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

    # NOTE: don't define any defaults here for parameters that are passed to EnPTConfig!
    #       -> otherwise, we cannot distinguish between explicity given parameters and default values
    #       => see docs in parsedArgs_to_user_opts() for explanation
    add('-jc', '--json_config', nargs='?', type=str,
        help='file path of a JSON file containing options. See here for an example: '
             'https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/'
             'EnPT/blob/master/enpt/options/options_default.json')
    add('--CPUs', type=int, default=None,
        help='number of CPU cores to be used for processing (default: "None" -> use all available')
    add('-im', '--path_l1b_enmap_image', type=str, default=None,
        help='input path of the EnMAP L1B image to be processed '
             '(zip-archive or root directory; must be given if not contained in --json-config.)')
    add('-imgap', '--path_l1b_enmap_image_gapfill', type=str,  default=None,
        help='input path of an adjacent EnMAP L1B image to be used for gap-filling (zip-archive or root directory)')
    add('-od', '--output_dir', type=str, default=None,
        help='output directory where processed data and log files are saved')
    add('-wd', '--working_dir', type=str, default=None,
        help='directory to be used for temporary files')
    add('-nla', '--n_lines_to_append', type=int, default=None,
        help='number of lines to be added to the main image [if None, use the whole imgap]. Requires --imgap to be set')
    add('--path_earthSunDist', type=str, default=None,
        help='input path of the earth sun distance model')
    add('--path_solar_irr', type=str,  default=None,
        help='input path of the solar irradiance model')
    add('--scale_factor_toa_ref', type=int, default=None,
        help='scale factor to be applied to TOA reflectance result')
    add('--enable_keystone_correction', type=int, default=False)
    add('--enable_vnir_swir_coreg', type=int, default=False)
    add('--path_reference_image', type=str, default=None)
    add('--sicor_cache_dir', type=str, default=None)
    add('--auto_download_ecmwf', type=bool, default=False)
    add('--enable_cloud_screening', type=bool, default=False)
    add('--scale_factor_boa_ref', type=int, default=10000)
    add('--run_smile_P', type=bool, default=False)
    add('--run_deadpix_P', type=bool, default=True)
    add('--deadpix_P_algorithm', type=str, default="spectral")
    add('--deadpix_P_interp', type=str, default="linear")
    add('--ortho_resampAlg', type=int, default=1)

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
    CTR = EnPT_Controller(config)
    CTR.run_all_processors()


if __name__ == '__main__':
    parsed_args = get_enpt_argparser().parse_args()
    parsed_args.func(get_config(parsed_args))

    print('\nready.')
