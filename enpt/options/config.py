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
from collections import OrderedDict, Mapping
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
    enable_ice_retrieval=False,
    CPUs=16
)


config_for_testing_dlr = dict(
    path_l1b_enmap_image=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
                     # Alps
                     'ENMAP01-____L1B-DT000000987_20130205T105307Z_001_V000101_20190426T143700Z__rows0-99.zip'

                     # Arcachon
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-799.zip'

                     # Arcachon 1000x30
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows700-730.zip'
                     )),
    path_l1b_enmap_image_gapfill=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B',
                     # Alps
                     'ENMAP01-____L1B-DT000000987_20130205T105307Z_001_V000101_20190426T143700Z__rows100-199.zip'

                     # Arcachon
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__rows800-899.zip'
                     )),
    path_dem=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data',
                     # Alps
                     'DLR_L2A_DEM_UTM32.bsq'

                     # Arcachon
                     # 'ENMAP01-____L1B-DT000400126_20170218T110115Z_002_V000204_20200206T182719Z__DEM_ASTER.bsq'
                     )),
    log_level='DEBUG',
    output_dir=os.path.join(path_enptlib,  '..', 'tests', 'data', 'test_outputs'),
    n_lines_to_append=50,
    disable_progress_bars=False,
    is_dummy_dataformat=False,
    enable_ac=True,
    enable_ice_retrieval=False,
    CPUs=1,
    ortho_resampAlg='gauss',
    vswir_overlap_algorithm='swir_only'
)


enmap_coordinate_grid = dict(x=np.array([0, 30]),
                             y=np.array([0, 30]))
enmap_xres, enmap_yres = np.ptp(enmap_coordinate_grid['x']), np.ptp(enmap_coordinate_grid['y'])
assert enmap_xres == enmap_yres, 'Unequal X/Y resolution of the output grid!'


class EnPTConfig(object):
    def __init__(self, json_config='', **user_opts):
        """Create a job configuration.

        :arg json_config:
             path to JSON file containing configuration parameters or a string in JSON format

        :key CPUs:
             number of CPU cores to be used for processing (default: "None" -> use all available

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

        :key working_dir:
            directory to be used for temporary files

        :key n_lines_to_append:
            number of lines to be added to the main image [if None, use the whole imgap].
            Requires 'path_l1b_enmap_image_gapfill' to be set.

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

        :key path_reference_image:
            Reference image for co-registration.

        :key enable_ac:
            Enable atmospheric correction using SICOR algorithm (default: True).
            If False, the L2A output contains top-of-atmosphere reflectance.

        :key auto_download_ecmwf:
            Automatically download ECMWF data for atmospheric correction

        :key enable_ice_retrieval:
            Enable ice retrieval (default); increases accuracy of water vapour retrieval

        :key enable_cloud_screening:
            Enable cloud screening during atmospheric correction

        :key scale_factor_boa_ref:
            Scale factor to be applied to BOA reflectance result

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
            warnings.warn("The 'is_dlr_dataformat' flag is deprectated and will not exist in future. "
                          "Please set 'is_dummy_dataformat' to False instead.", DeprecationWarning)
            self.is_dummy_dataformat = user_opts['is_dlr_dataformat'] is False

        self.CPUs = gp('CPUs', fallback=cpu_count())
        self.log_level = gp('log_level')
        self.create_logfile = gp('create_logfile')
        self.path_l1b_enmap_image = self.absPath(gp('path_l1b_enmap_image'))
        self.path_l1b_enmap_image_gapfill = self.absPath(gp('path_l1b_enmap_image_gapfill'))
        self.path_dem = self.absPath(gp('path_dem'))
        self.average_elevation = self.absPath(gp('average_elevation'))
        self.path_l1b_snr_model = self.absPath(gp('path_l1b_snr_model'))
        self.working_dir = self.absPath(gp('working_dir')) or None
        self.n_lines_to_append = gp('n_lines_to_append')
        self.disable_progress_bars = gp('disable_progress_bars')

        ##################
        # output options #
        ##################

        self.output_dir = self.absPath(gp('output_dir', fallback=os.path.abspath(os.path.curdir)))

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
        self.path_reference_image = gp('path_reference_image')

        # atmospheric_correction
        self.enable_ac = gp('enable_ac')
        self.auto_download_ecmwf = gp('auto_download_ecmwf')
        self.enable_ice_retrieval = gp('enable_ice_retrieval')
        self.enable_cloud_screening = gp('enable_cloud_screening')
        self.scale_factor_boa_ref = gp('scale_factor_boa_ref')

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

        #########################
        # validate final config #
        #########################

        EnPTValidator(allow_unknown=True, schema=enpt_schema_config_output).validate(self.to_dict())

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
        if value in [True, "true"]:
            return True
        if value in [False, "false"]:
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
        """

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
