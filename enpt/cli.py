#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz-potsdam.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz-potsdam.de),
# Stephane Guillaso (GFZ Potsdam, stephane.guillaso@gfz-potsdam.de)
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

"""EnPT console argument parser."""

import argparse
import sys

from enpt import __version__
from enpt.options.config import EnPTConfig
from enpt.execution.controller import EnPT_Controller

__author__ = 'Daniel Scheffler'


def get_enpt_argparser():
    """Return argument parser for the 'enpt' program."""

    ##########################################################
    # CONFIGURE MAIN PARSER FOR THE EnPT PREPROCESSING CHAIN #
    ##########################################################

    parser = argparse.ArgumentParser(
        prog='enpt',
        description='=' * 70 + '\n' + 'EnMAP Processing Tool console argument parser. ',
        epilog="use '>>> enpt -h' for detailed documentation and usage hints.")

    add = parser.add_argument
    add('--version', action='version', version=__version__)

    # NOTE: don't define any defaults here for parameters that are passed to EnPTConfig!
    #       -> otherwise, we cannot distinguish between explicity given parameters and default values
    #       => see docs in parsedArgs_to_user_opts() for explanation
    add('-jc', '--json_config', nargs='?', type=str,
        help='file path of a JSON file containing options. See here for an example: '
             'https://git.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/'
             'EnPT/blob/master/enpt/options/options_default.json')
    add('--CPUs', type=int, default=None,
        help='number of CPU cores to be used for processing (default: "None" -> use all available)')
    add('-im', '--path_l1b_enmap_image', type=str, default=None,
        help='input path of the EnMAP L1B image to be processed '
             '(zip-archive or root directory; must be given if not contained in --json-config.)')
    add('-imgap', '--path_l1b_enmap_image_gapfill', type=str,  default=None,
        help='input path of an adjacent EnMAP L1B image to be used for gap-filling (zip-archive or root directory)')
    add('-dem', '--path_dem', type=str,  default=None,
        help='input path of digital elevation model in map or sensor geometry; GDAL compatible file format '
             '(must cover the EnMAP L1B data completely if given in map geometry or must have the same pixel '
             'dimensions like the EnMAP L1B data if given in sensor geometry)')
    add('-dummyfmt', '--is_dummy_dataformat', type=_str2bool, default=False, nargs='?', const=True,
        help='Set to true in case of the preliminary, GFZ-internal dataformat as used for the Alpine test dataset. '
             '(default: False. Note: This will be removed in future.)')
    add('-ele', '--average_elevation', type=int, default=0,
        help='average elevation in meters above sea level; may be provided if no DEM is available; '
             'ignored if DEM is given')
    add('-od', '--output_dir', type=str, default=None,
        help='output directory where processed data and log files are saved')
    add('-of', '--output_format', type=str, default='GTiff',
        help="file format of all raster output files ('GTiff': GeoTIFF, 'ENVI':  ENVI BSQ; default: 'ENVI')")
    add('-ointlv', '--output_interleave', type=str, default='pixel',
        help="raster data interleaving type ('band', 'line', 'pixel'; default: 'pixel')")
    add('-wd', '--working_dir', type=str, default=None,
        help='directory to be used for temporary files')
    add('-nla', '--n_lines_to_append', type=int, default=None,
        help='number of lines to be added to the main image [if None, use the whole imgap]. Requires --imgap to be set')
    add('-dbb', '--drop_bad_bands', type=_str2bool, default=True,
        help='if set to True (default), the water absorption bands between 1358 and 1453 nm as well as between 1814 '
             'and 1961 nm are excluded from processing and will not be contained in the L2A product')
    add('-dpb', '--disable_progress_bars', type=_str2bool, default=False, nargs='?', const=True,
        help='whether to disable all progress bars during processing')
    add('--path_earthSunDist', type=str, default=None,
        help='input path of the earth sun distance model')
    add('--path_solar_irr', type=str,  default=None,
        help='input path of the solar irradiance model')
    add('--scale_factor_toa_ref', type=int, default=None,
        help='scale factor to be applied to TOA reflectance result')
    add('--enable_keystone_correction', type=_str2bool, default=False, nargs='?', const=True,
        help='Enable keystone correction')
    add('--enable_vnir_swir_coreg', type=_str2bool, default=False, nargs='?', const=True,
        help='Enable VNIR/SWIR co-registration')
    add('--path_reference_image', type=str, default=None,
        help='Reference image for co-registration.')
    add('--polymer_root', type=str, default=None,
        help='Polymer root directory (that contains the subdirectory for ancillary data)')
    add('--enable_ac', type=_str2bool, default=True, nargs='?', const=True,
        help="Enable atmospheric correction using SICOR algorithm (default: True). If False, the L2A output contains "
             "top-of-atmosphere reflectance")
    add('--mode_ac', type=str, default=None, nargs='?',
        help="3 modes to determine which atmospheric correction is applied at which surfaces (default: land): "
             "('land', water', 'combined')")
    add('--auto_download_ecmwf', type=_str2bool, default=True, nargs='?', const=True,
        help='Automatically download ECMWF AUX data when running Polymer atmospheric correction for water surfaces')
    add('--scale_factor_boa_ref', type=int, default=10000,
        help='Scale factor to be applied to BOA reflectance result')
    add('--threads', type=int, default=-1,
        help='Number of threads in ACwater Polymer: 0 for single thread; < 0 for as many as there are CPUs; '
             'and > 0 gives the number of threads')
    add('--blocksize', type=int, default=100,
        help='Block size in ACwater Polymer')
    add('--run_smile_P', type=_str2bool, default=False, nargs='?', const=True,
        help='Enable extra smile detection and correction (provider smile coefficients are ignored)')
    add('--run_deadpix_P', type=_str2bool, default=True, nargs='?', const=True,
        help='Enable dead pixel correction')
    add('--deadpix_P_algorithm', type=str, default="spectral",
        help="Algorithm for dead pixel correction ('spectral' or 'spatial')")
    add('--deadpix_P_interp_spectral', type=str, default="linear",
        help="Spectral interpolation algorithm to be used during dead pixel correction "
             "('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')")
    add('--deadpix_P_interp_spatial', type=str, default="linear",
        help="Spatial interpolation algorithm to be used during dead pixel correction "
             "('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')")
    add('--ortho_resampAlg', type=str, default='bilinear',
        help="Ortho-rectification resampling algorithm ('nearest', 'bilinear', 'gauss')")
    add('--vswir_overlap_algorithm', type=str, default='swir_only',
        help="Algorithm specifying how to deal with the spectral bands in the VNIR/SWIR spectral overlap region "
             "('order_by_wvl', 'average', 'vnir_only', 'swir_only')")
    add('-tgtprj', '--target_projection_type', type=str, default='UTM',
        help="Projection type of the raster output files ('UTM', 'Geographic') (default: 'UTM')")
    add('-tgtepsg', '--target_epsg', type=int, default=None,
        help="Custom EPSG code of the target projection (overrides target_projection_type)")
    add('-tgtgrid', '--target_coord_grid', nargs=4, type=float, default=None,
        help="Custom target coordinate grid where is output is resampled to ([x0, x1, y0, y1], e.g., [0, 30, 0, 30])")

    # link parser to run function
    parser.set_defaults(func=run_job)

    return parser


def parsedArgs_to_user_opts(cli_args: argparse.Namespace) -> dict:
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

    # activate absolute coreg if a reference image is given (the argparser does not separate these two options)
    if opts['path_reference_image']:
        opts['enable_absolute_coreg'] = True

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


def _str2bool(v):
    """Convert string parameter to bool.

    From: https://stackoverflow.com/a/43357954/2952871

    :param v:
    :return:
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(parsed_args: argparse.Namespace = None) -> int:
    """Run the argument parser and forward the arguments to the linked functions.

    :param parsed_args:     argparse.Namespace instance of already parsed arguments
                            (allows to call main() from test_cli_parser.py while passing specific arguments)

    :return:  exitcode (0: all fine, >=1: failed)
    """
    if not parsed_args:
        parsed_args = get_enpt_argparser().parse_args()  # type: argparse.Namespace

    parsed_args.func(get_config(parsed_args))

    print('\nready.')

    return 0


if __name__ == '__main__':
    sys.exit(main())  # pragma: no cover
