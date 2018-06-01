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

import sicor

from .options_schema import \
    enpt_schema_input, \
    enpt_schema_config_output, \
    parameter_mapping, \
    get_param_from_json_config
from ..version import \
    __version__, \
    __versionalias__


path_enptlib = os.path.dirname(pkgutil.get_loader("enpt").path)
path_options_default = os.path.join(path_enptlib, 'options', 'options_default.json')


config_for_testing = dict(
    path_l1b_enmap_image=os.path.abspath(
        os.path.join(path_enptlib, '..', 'tests', 'data', 'EnMAP_Level_1B', 'AlpineTest1_CWV2_SM0.zip')),
    log_level='DEBUG',
    output_dir=os.path.join(path_enptlib,  '..', 'tests', 'data', 'test_outputs')
)


class EnPTConfig(object):
    def __init__(self, json_config='', **user_opts):
        """Create a job configuration.

        :param json_config  path to JSON file containing configuration parameters or a string in JSON format
        :param user_opts    keyword arguments
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

        self.CPUs = gp('CPUs', fallback=cpu_count())
        self.log_level = gp('log_level')
        self.create_logfile = gp('create_logfile')
        self.path_l1b_enmap_image = self.absPath(gp('path_l1b_enmap_image'))
        self.path_l1b_enmap_image_gapfill = self.absPath(gp('path_l1b_enmap_image_gapfill'))
        self.path_l1b_snr_model = self.absPath(gp('path_l1b_snr_model'))
        self.working_dir = self.absPath(gp('working_dir')) or None

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
        self.sicor_cache_dir = gp('sicor_cache_dir', fallback=sicor.__path__[0])
        self.auto_download_ecmwf = gp('auto_download_ecmwf')
        self.enable_cloud_screening = gp('enable_cloud_screening')
        self.scale_factor_boa_ref = gp('scale_factor_boa_ref'),

        # smile
        self.run_smile_P = gp('run_smile_P')

        # dead_pixel
        self.run_deadpix_P = gp('run_deadpix_P')
        self.deadpix_P_algorithm = gp('deadpix_P_algorithm')
        self.deadpix_P_interp = gp('deadpix_P_interp')

        # orthorectification
        self.ortho_resampAlg = gp('ortho_resampAlg')

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


def get_options(target: str, validation: bool=True):
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
