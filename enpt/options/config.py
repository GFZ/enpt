# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2021 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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

"""EnPT configuration module.

Provides the configuration that is later passed to individual submodules.
"""

import os
import json
from json import JSONDecodeError
import datetime
import pkgutil
import warnings
from pprint import pformat

from jsmin import jsmin
from cerberus import Validator
from collections import OrderedDict
from collections.abc import Mapping
import numpy as np
from multiprocessing import cpu_count

from .options_schema import \
    enpt_schema_input, \
    enpt_schema_config_output, \
    parameter_mapping, \
    get_param_from_json_config
from ..version import \
    __version__, \
    __versionalias__

__author__ = 'Daniel Scheffler'


path_enptlib = os.path.dirname(pkgutil.get_loader("enpt").path)
path_options_default = os.path.join(path_enptlib, 'options', 'options_default.json')

try:
    # from acwater.acwater import polymer_ac_enmap
    path_polymer = os.path.abspath(os.path.join(os.path.dirname(pkgutil.get_loader("polymer").path), os.pardir))
except AttributeError:
    path_polymer = ''

config_for_testing_water = dict(
    path_l1b_enmap_image=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
                     # Arcachon
                     'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-730.zip'

                     # Arcachon full tile 2
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z.zip'
                     )),
    # path_l1b_enmap_image_gapfill=os.path.abspath(
    #     os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
    #                  'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-730.zip')),
    path_dem=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data',
                     'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__tile2'
                     '__DEM_ASTER.bsq')),
    log_level='DEBUG',
    output_dir=os.path.join(path_enptlib, '..', 'tests', 'data', 'test_outputs'),
    disable_progress_bars=False,
    is_dummy_dataformat=False,
    auto_download_ecmwf=True,
    average_elevation=0,
    deadpix_P_algorithm='spectral',
    deadpix_P_interp_spatial='linear',
    deadpix_P_interp_spectral='linear',
    enable_keystone_correction=False,
    enable_vnir_swir_coreg=False,
    n_lines_to_append=None,
    ortho_resampAlg='gauss',
    run_deadpix_P=True,
    run_smile_P=False,
    scale_factor_boa_ref=10000,
    scale_factor_toa_ref=10000,
    enable_ac=True,
    mode_ac='combined',
    polymer_root=path_polymer,
    threads=-1,
    blocksize=100,
    vswir_overlap_algorithm='swir_only',
    CPUs=16
)


config_for_testing = dict(
    path_l1b_enmap_image=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B', 'AlpineTest1_CWV2_SM0.zip')),
    path_l1b_enmap_image_gapfill=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B', 'AlpineTest2_CWV2_SM0.zip')),
    path_dem=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'dem_map_geo.bsq')),
    log_level='DEBUG',
    output_dir=os.path.join(path_enptlib,  '..', 'tests', 'data', 'test_outputs'),
    n_lines_to_append=50,
    disable_progress_bars=True,
    is_dummy_dataformat=True,
    enable_ac=False,
    ortho_resampAlg='bilinear',
    CPUs=16
)


config_for_testing_dlr = dict(
    path_l1b_enmap_image=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
                     # Alps
                     # 'ENMAP01-____L1B-DT000000987_20130205T105307Z_001_V000101_20190426T143700Z__rows0-99.zip'

                     # Alps full
                     # 'ENMAP01-____L1B-DT000000987_20130205T105307Z_001_V000101_20190426T143700Z.zip'

                     # Arcachon
                     'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-799.zip'

                     # Arcachon 1000x30
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-730.zip'

                     # Arcachon full tile 2
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z.zip'

                     # Arcachon full tile 3, reprocessed 05/2020
                     # 'ENMAP01-____L1B-DT000400126_20170218T110119Z_003_V000204_20200508T124425Z.zip'

                     # Arcachon tile 3 (full), downloaded from enmap.org
                     # 'L1B_Arcachon_3__enmap.org.zip',
                     )),
    # path_l1b_enmap_image_gapfill=os.path.abspath(
    #     os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
    #                  # Alps
    #                  'ENMAP01-____L1B-DT000000987_20130205T105307Z_001_V000101_20190426T143700Z__rows100-199.zip'
    #
    #                  # Arcachon
    #                  # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows800-899.zip'
    #                  )),
    path_dem=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data',
                     # Alps
                     # 'DLR_L2A_DEM_UTM32.bsq'

                     # Arcachon tile 2 ASTER DEM (02/2020)
                     'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__tile2__DEM_ASTER.bsq'

                     # Arcachon tile 3 ASTER DEM (05/2020)
                     # 'ENMAP01-____L1B-DT000400126_20170218T110119Z_003_V000204_20200508T124425Z__tile3__DEM_ASTER.bsq'
                     # '15_DEM_UTM__with_prj.tif'
                     )),
    log_level='DEBUG',
    output_dir=os.path.join(path_enptlib,  '..', 'tests', 'data', 'test_outputs'),
    n_lines_to_append=50,
    disable_progress_bars=False,
    is_dummy_dataformat=False,
    # output_format='ENVI',
    # output_interleave='band',
    # target_projection_type='Geographic',
    # target_epsg=32632,
    # target_coord_grid=[-1.37950, -1.37923, 44.60710, 44.60737],
    enable_absolute_coreg=True,
    path_reference_image=os.path.join(path_enptlib, '..', 'tests', 'data', 'T30TXQ_20170218T110111_B05__sub.tif'),
    enable_ac=True,
    mode_ac='land',
    CPUs=32,
    ortho_resampAlg='gauss',
    vswir_overlap_algorithm='swir_only'
)


enmap_coordinate_grid_utm = dict(x=np.array([0, 30]),
                                 y=np.array([0, 30]))
enmap_xres, enmap_yres = np.ptp(enmap_coordinate_grid_utm['x']), np.ptp(enmap_coordinate_grid_utm['y'])


class EnPTConfig(object):
    def __init__(self, json_config='', **user_opts):
        """Create a job configuration.

        :arg json_config:
             path to JSON file containing configuration parameters or a string in JSON format

        :key CPUs:
             number of CPU cores to be used for processing (default: "None" -> use all available)

        :key path_l1b_enmap_image:
            input path of the EnMAP L1B image to be processed
            (zip-archive or root directory; must be given if not contained in --json-config.)

        :key path_l1b_enmap_image_gapfill:
            input path of an adjacent EnMAP L1B image to be used for gap-filling (zip-archive or root directory)

        :key path_dem:
            input path of digital elevation model in map or sensor geometry; GDAL compatible file format (must cover
            the EnMAP L1B data completely if given in map geometry or must have the same pixel dimensions like the
            EnMAP L1B data if given in sensor geometry)

        :key average_elevation:
            average elevation in meters above sea level; may be provided if no DEM is available; ignored if DEM is given

        :key output_dir:
            output directory where processed data and log files are saved

        :key output_format:
            file format of all raster output files ('GTiff': GeoTIFF, 'ENVI':  ENVI BSQ; default: 'ENVI')

        :key output_interleave:
            raster data interleaving type (default: 'pixel')
            - 'band': band-sequential (BSQ),
            - 'line': data interleaved-by-line (BIL; only usable for ENVI output format),
            - 'pixel' data interleaved-by-pixel (BIP)

        :key working_dir:
            directory to be used for temporary files

        :key n_lines_to_append:
            number of lines to be added to the main image [if None, use the whole imgap].
            Requires 'path_l1b_enmap_image_gapfill' to be set.

        :key drop_bad_bands:
            if set to True (default), the water absorption bands between 1358 and 1453 nm as well
            as between 1814 and 1961 nm are excluded from processing and will not be contained in the L2A product

        :key disable_progress_bars:
            whether to disable all progress bars during processing

        :key path_earthSunDist:
             input path of the earth sun distance model

        :key path_solar_irr:
            input path of the solar irradiance model

        :key scale_factor_toa_ref:
            scale factor to be applied to TOA reflectance result

        :key enable_keystone_correction:
            Enable keystone correction

        :key enable_vnir_swir_coreg:
            Enable VNIR/SWIR co-registration

        :key enable_absolute_coreg:
            Enable the co-registration of the EnMAP image to the reference image given with 'path_reference_image'

        :key path_reference_image:
            Reference image for co-registration.

        :key polymer root:
            Polymer root directory (that contains the subdirectory for ancillary data).

        :key enable_ac:
            Enable atmospheric correction using SICOR algorithm (default: True).
            If False, the L2A output contains top-of-atmosphere reflectance.

        :key mode_ac:
            3 modes to determine which atmospheric correction is applied at which surfaces (default: land):
            - 'land': SICOR (developed for land surfaces is applied to land AND water surfaces
            - 'water': POLYMER (developed for water surfaces) is applied to water only
                       (land surfaces are no included in the L2A product)
            - 'combined': SICOR is applied to land and POLYMER is applied to water surfaces;
                          NOTE that this may result in edge effects, e.g., at coastlines

        :key auto_download_ecmwf:
            Automatically download ECMWF AUX data when running Polymer atmospheric correction for water surfaces

        :key scale_factor_boa_ref:
            Scale factor to be applied to BOA reflectance result

        :key threads:
            number of threads for multiprocessing of blocks (see bellow):
            - 'threads = 0': for single thread
            - 'threads < 0': for as many threads as there are CPUs
            - 'threads > 0': gives the number of threads

        :key blocksize:
            block size for multiprocessing

        :key run_smile_P:
            Enable extra smile detection and correction (provider smile coefficients are ignored)

        :key run_deadpix_P:
            Enable dead pixel correction

        :key deadpix_P_algorithm:
            Algorithm for dead pixel correction ('spectral' or 'spatial')

        :key deadpix_P_interp_spectral:
            Spectral interpolation algorithm to be used during dead pixel correction
             ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')

        :key deadpix_P_interp_spatial:
            Spatial interpolation algorithm to be used during dead pixel correction
             ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')

        :key ortho_resampAlg:
            Ortho-rectification resampling algorithm ('nearest', 'bilinear', 'gauss')

        :key target_projection_type:
            Projection type of the raster output files ('UTM', 'Geographic') (default: 'UTM')

        :key target_epsg:
            Custom EPSG code of the target projection (overrides target_projection_type)

        :key target_coord_grid:
            Custom target coordinate grid where the output is resampled to ([x0, x1, y0, y1], e.g., [0, 30, 0, 30])
        """
        # fixed attributes
        self.version = __version__
        self.versionalias = __versionalias__

        #######################
        # POPULATE PARAMETERS #
        #######################

        # args
        self.json_config = json_config
        self.kwargs = user_opts

        # get validated options dict from JSON-options
        self.json_opts_fused_valid = self.get_json_opts(validate=True)

        gp = self.get_parameter

        ###################
        # general options #
        ###################

        self.is_dummy_dataformat = gp('is_dummy_dataformat')
        if 'is_dlr_dataformat' in user_opts:
            warnings.warn("The 'is_dlr_dataformat' flag is deprecated and will not exist in future. "
                          "Please set 'is_dummy_dataformat' to False instead.", DeprecationWarning)
            self.is_dummy_dataformat = user_opts['is_dlr_dataformat'] is False

        self.CPUs = gp('CPUs', fallback=cpu_count())
        self.log_level = gp('log_level')
        self.create_logfile = gp('create_logfile')
        self.path_l1b_enmap_image = self.absPath(gp('path_l1b_enmap_image'))
        self.path_l1b_enmap_image_gapfill = self.absPath(gp('path_l1b_enmap_image_gapfill'))
        self.path_dem = self.absPath(gp('path_dem'))
        self.average_elevation = gp('average_elevation')
        self.path_l1b_snr_model = self.absPath(gp('path_l1b_snr_model'))
        self.working_dir = self.absPath(gp('working_dir')) or None
        self.n_lines_to_append = gp('n_lines_to_append')
        self.drop_bad_bands = gp('drop_bad_bands')
        self.disable_progress_bars = gp('disable_progress_bars')

        ##################
        # output options #
        ##################

        self.output_dir = self.absPath(gp('output_dir', fallback=os.path.abspath(os.path.curdir)))
        self.output_format = gp('output_format')
        self.output_interleave = gp('output_interleave')

        ###########################
        # processor configuration #
        ###########################

        # toa_ref
        self.path_earthSunDist = self.absPath(gp('path_earthSunDist'))
        self.path_solar_irr = self.absPath(gp('path_solar_irr'))
        self.scale_factor_toa_ref = gp('scale_factor_toa_ref')

        # geometry
        self.enable_keystone_correction = gp('enable_keystone_correction')
        self.enable_vnir_swir_coreg = gp('enable_vnir_swir_coreg')
        self.enable_absolute_coreg = gp('enable_absolute_coreg')
        self.path_reference_image = gp('path_reference_image')

        # atmospheric_correction
        self.polymer_root = gp('polymer_root')
        self.enable_ac = gp('enable_ac')
        self.mode_ac = gp('mode_ac')
        self.auto_download_ecmwf = gp('auto_download_ecmwf')
        self.scale_factor_boa_ref = gp('scale_factor_boa_ref')
        self.threads = gp('threads')
        self.blocksize = gp('blocksize')

        # smile
        self.run_smile_P = gp('run_smile_P')

        # dead_pixel
        self.run_deadpix_P = gp('run_deadpix_P')
        self.deadpix_P_algorithm = gp('deadpix_P_algorithm')
        self.deadpix_P_interp_spectral = gp('deadpix_P_interp_spectral')
        self.deadpix_P_interp_spatial = gp('deadpix_P_interp_spatial')

        # orthorectification / VSWIR fusion
        self.ortho_resampAlg = gp('ortho_resampAlg')
        self.vswir_overlap_algorithm = gp('vswir_overlap_algorithm')
        self.target_projection_type = gp('target_projection_type')
        self.target_epsg = gp('target_epsg')
        grid = gp('target_coord_grid')
        self.target_coord_grid = dict(x=np.array(grid[:2]), y=np.array(grid[2:])) if grid else None

        #########################
        # validate final config #
        #########################

        EnPTValidator(allow_unknown=True, schema=enpt_schema_config_output).validate(self.to_dict())

        # check if given paths point to existing files
        if os.getenv('IS_ENPT_GUI_TEST') != "1":
            paths = {k: v for k, v in self.__dict__.items() if k.startswith('path_')}
            for k, fp in paths.items():
                if fp and not os.path.isfile(fp):
                    raise FileNotFoundError("The file path provided at the '%s' parameter does not point "
                                            "to an existing file (%s)." % (k, fp))

        if not self.path_dem:
            warnings.warn('No digital elevation model provided. Note that this may cause uncertainties, e.g., '
                          'in the atmospheric correction and the orthorectification.', RuntimeWarning, stacklevel=2)

        # check invalid interleave
        if self.output_interleave == 'line' and self.output_format == 'GTiff':
            warnings.warn("The interleaving type 'line' is not supported by the GTiff output format. Using 'pixel'.",
                          UserWarning)
            self.output_interleave = 'pixel'

        # override target_projection_type if target_epsg is given
        if self.target_epsg:
            self.target_projection_type = \
                'Geographic' if self.target_epsg == 4326 else \
                'UTM' if len(str(self.target_epsg)) == 5 and str(self.target_epsg)[:3] in ['326', '327'] else \
                'NA'
        if self.target_projection_type == 'Geographic':
            self.target_epsg = 4326

        # set target coordinate grid to the UTM EnMAP grid if no other grid is provided and target projection is UTM
        self.target_coord_grid = \
            self.target_coord_grid if self.target_coord_grid else \
            enmap_coordinate_grid_utm if self.target_projection_type == 'UTM' else None

        # bug warning regarding holes in bilinear resampling output
        if self.target_projection_type == 'Geographic' and self.ortho_resampAlg == 'bilinear':
            warnings.warn("There is currently a bug that causes holes in the bilinear resampling results if the "
                          "target projection is 'Geographic'. It is recommended to use 'nearest' or 'gauss' instead.",
                          UserWarning)

    @staticmethod
    def absPath(path):
        return path if not path or os.path.isabs(path) else os.path.abspath(os.path.join(path_enptlib, path))

    def get_parameter(self, key_user_opts, fallback=None):
        # 1. priority: parameters that have directly passed to EnPTConfig within user_opts
        if key_user_opts in self.kwargs:
            return self.kwargs[key_user_opts]

        # 2. priority: default options, overridden by eventually provided json_config
        else:
            param = get_param_from_json_config(key_user_opts, self.json_opts_fused_valid)
            if not param:
                if fallback:
                    return fallback
            return param

    def get_json_opts(self, validate=True):
        """Get a dictionary of EnPT config parameters.

        NOTE: Reads the default options from options_default.json and updates the values with those from database.
        """
        def update_dict(d, u):
            for k, v in u.items():
                if isinstance(v, Mapping):
                    d[k] = update_dict(d.get(k, {}), v)
                else:
                    d[k] = v
            return d

        # read options_default.json
        default_options = get_options(path_options_default, validation=validate)

        ###############################################################################################################
        # if json config is provided (via python bindings or CLI parser -> override all options with that json config #
        ###############################################################################################################

        if self.json_config:
            if self.json_config.startswith("{"):
                try:
                    params_dict = json.loads(jsmin(self.json_config))
                except JSONDecodeError:
                    warnings.warn('The given JSON options string could not be decoded. '
                                  'JSON decoder failed with the following error:')
                    raise
            elif os.path.isfile(self.json_config):
                try:
                    with open(self.json_config, 'r') as inF:
                        params_dict = json.loads(jsmin(inF.read()))
                except JSONDecodeError:
                    warnings.warn('The given JSON options file %s could not be decoded. '
                                  'JSON decoder failed with the following error:' % self.json_config)
                    raise

            else:
                raise ValueError("The parameter 'json_config' must be a JSON formatted string or a JSON file on disk.")

            # convert values to useful data types and update the default values
            params_dict = json_to_python(params_dict)
            update_dict(default_options, params_dict)

        if validate:
            EnPTValidator(allow_unknown=True, schema=enpt_schema_input).validate(default_options)

        json_options = default_options
        return json_options

    def to_dict(self):
        """Generate a dictionary in the same structure like the one in options_default.json from the current config."""

        def nested_set(dic, keys, value):
            for k in keys[:-1]:
                dic = dic.setdefault(k, {})
            dic[keys[-1]] = value

        outdict = dict()
        for key_user_opts, subkeys in parameter_mapping.items():
            nested_set(outdict, subkeys, getattr(self, key_user_opts))

        return outdict

    def to_jsonable_dict(self):
        return python_to_json(self.to_dict())

    def save(self, path_outfile):
        """Save the JobConfig instance to a JSON file in the same structure like the one in options_default.json.

        :param path_outfile:    path of the output JSON file
        """
        with open(path_outfile, 'w') as outF:
            json.dump(self.to_jsonable_dict(), outF, skipkeys=False, indent=4)

    def __repr__(self):
        return pformat(self.to_dict())


def json_to_python(value):
    def is_number(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    if type(value) is dict:
        return {json_to_python(k): json_to_python(v) for k, v in value.items()}
    elif type(value) is list:
        return [json_to_python(v) for v in value]
    else:
        if value == "None":
            return None
        if value == "slice(None, None, None)":
            return slice(None)
        if value is True or value == "true":
            return True
        if value is False or value == "false":
            return False
        if is_number(value):
            try:
                if str(int(value)) != str(float(value)):
                    return int(value)
                else:
                    return float(value)
            except ValueError:
                return float(value)
        else:
            return value


def python_to_json(value):
    if type(value) in [dict, OrderedDict]:
        return {python_to_json(k): python_to_json(v) for k, v in value.items()}
    elif type(value) is list:
        return [python_to_json(v) for v in value]
    elif type(value) is np.ndarray:
        return [python_to_json(v) for v in value.tolist()]
    else:
        if value is None:
            return "None"
        if value is slice(None):
            return "slice(None, None, None)"
        if value is True:
            return "true"
        if value is False:
            return "false"
        if type(value) is datetime.datetime:
            return datetime.datetime.strftime(value, '%Y-%m-%d %H:%M:%S.%f%z')
        else:
            return value


class EnPTValidator(Validator):
    def __init__(self, *args, **kwargs):
        """Get an instance of EnPTValidator.

        :param args:    Arguments to be passed to cerberus.Validator
        :param kwargs:  Keyword arguments to be passed to cerberus.Validator
        """
        super(EnPTValidator, self).__init__(*args, **kwargs)

    def validate(self, document2validate, **kwargs):
        if super(EnPTValidator, self).validate(document=document2validate, **kwargs) is False:
            raise ValueError("Options is malformed: %s" % str(self.errors))


def get_options(target: str, validation: bool = True):
    """Return dictionary with all options.

    :param target:      if path to file, then json is used to load, otherwise the default template is used
    :param validation:  True / False, whether to validate options read from files or not
    :return: dictionary with options
    """
    if os.path.isfile(target):
        with open(target, "r") as fl:
            options = json_to_python(json.loads(jsmin(fl.read())))

        if validation is True:
            EnPTValidator(allow_unknown=True, schema=enpt_schema_input).validate(options)

        return options
    else:
        raise FileNotFoundError("Options file not found at file path %s." % target)
