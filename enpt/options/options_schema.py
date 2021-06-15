# -*- coding: utf-8

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

"""Definition of EnPT options schema (as used by cerberus library)."""


enpt_schema_input = dict(

    general_opts=dict(
        type='dict', required=False,
        schema=dict(
            CPUs=dict(type='integer', required=False, nullable=True),
            log_level=dict(type='string', required=False, allowed=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
            create_logfile=dict(type='boolean', required=False),
            path_l1b_enmap_image=dict(type='string', required=False),
            path_l1b_enmap_image_gapfill=dict(type='string', required=False),
            path_dem=dict(type='string', required=False),
            is_dummy_dataformat=dict(type='boolean', required=False),
            average_elevation=dict(type='integer', required=False),
            path_l1b_snr_model=dict(type='string', required=False),
            working_dir=dict(type='string', required=False, nullable=True),
            n_lines_to_append=dict(type='integer', required=False, nullable=True, min=0),
            drop_bad_bands=dict(type='boolean', required=False),
            disable_progress_bars=dict(type='boolean', required=False, nullable=True),
        )),

    output=dict(
        type='dict', required=False,
        schema=dict(
            output_dir=dict(type='string', required=False),
            output_format=dict(type='string', required=False, allowed=['GTiff', 'ENVI']),
            output_interleave=dict(type='string', required=False, allowed=['band', 'line', 'pixel'])
        )),

    processors=dict(
        type='dict', required=False,
        schema=dict(

            toa_ref=dict(
                type='dict', required=False,
                schema=dict(
                    path_earthSunDist=dict(type='string', required=False),
                    path_solar_irr=dict(type='string', required=False),
                    scale_factor_toa_ref=dict(type='integer', required=False, min=1),
                )),

            geometry=dict(
                type='dict', required=False,
                schema=dict(
                    enable_keystone_correction=dict(type='boolean', required=False),
                    enable_vnir_swir_coreg=dict(type='boolean', required=False),
                    enable_absolute_coreg=dict(type='boolean', required=False),
                    path_reference_image=dict(type='string', required=False),
                )),

            atmospheric_correction=dict(
                type='dict', required=False,
                schema=dict(
                    polymer_root=dict(type='string', required=False),
                    enable_ac=dict(type='boolean', required=False),
                    mode_ac=dict(type='string', required=False, allowed=['land', 'water', 'combined']),
                    auto_download_ecmwf=dict(type='boolean', required=False),
                    scale_factor_boa_ref=dict(type='integer', required=False, min=1),
                    threads=dict(type='integer', required=False),
                    blocksize=dict(type='integer', required=False),

                )),

            smile=dict(
                type='dict', required=False,
                schema=dict(
                    run_processor=dict(type='boolean', required=False),
                )),

            dead_pixel=dict(
                type='dict', required=False,
                schema=dict(
                    run_processor=dict(type='boolean', required=False),
                    algorithm=dict(type='string', required=False, allowed=['spectral', 'spatial']),
                    interp_method_spectral=dict(type='string', required=False,
                                                allowed=['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']),
                    interp_method_spatial=dict(type='string', required=False,
                                               allowed=['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']),
                )),

            orthorectification=dict(
                type='dict', required=False,
                schema=dict(
                    resamp_alg=dict(type='string', required=False, allowed=['nearest', 'bilinear', 'gauss']),
                    vswir_overlap_algorithm=dict(type='string', required=False,
                                                 allowed=['order_by_wvl', 'average', 'vnir_only', 'swir_only']),
                    target_projection_type=dict(type='string', required=False, allowed=['UTM', 'Geographic']),
                    target_epsg=dict(type='integer', required=False, nullable=True, min=0, forbidden=[0]),
                    target_coord_grid=dict(type='list', required=False, nullable=True, minlength=4, maxlength=4)
                ))
        ))
)


parameter_mapping = dict(
    # general opts
    CPUs=('general_opts', 'CPUs'),
    log_level=('general_opts', 'log_level'),
    create_logfile=('general_opts', 'create_logfile'),
    path_l1b_enmap_image=('general_opts', 'path_l1b_enmap_image'),
    path_l1b_enmap_image_gapfill=('general_opts', 'path_l1b_enmap_image_gapfill'),
    path_dem=('general_opts', 'path_dem'),
    is_dummy_dataformat=('general_opts', 'is_dummy_dataformat'),
    average_elevation=('general_opts', 'average_elevation'),
    path_l1b_snr_model=('general_opts', 'path_l1b_snr_model'),
    working_dir=('general_opts', 'working_dir'),
    n_lines_to_append=('general_opts', 'n_lines_to_append'),
    drop_bad_bands=('general_opts', 'drop_bad_bands'),
    disable_progress_bars=('general_opts', 'disable_progress_bars'),

    # output
    output_dir=('output', 'output_dir'),
    output_format=('output', 'output_format'),
    output_interleave=('output', 'output_interleave'),

    # processors > toa_ref
    path_earthSunDist=('processors', 'toa_ref', 'path_earthSunDist'),
    path_solar_irr=('processors', 'toa_ref', 'path_solar_irr'),
    scale_factor_toa_ref=('processors', 'toa_ref', 'scale_factor_toa_ref'),

    # processors > geometry
    enable_keystone_correction=('processors', 'geometry', 'enable_keystone_correction'),
    enable_vnir_swir_coreg=('processors', 'geometry', 'enable_vnir_swir_coreg'),
    enable_absolute_coreg=('processors', 'geometry', 'enable_absolute_coreg'),
    path_reference_image=('processors', 'geometry', 'path_reference_image'),

    # processors > atmospheric_correction
    polymer_root=('processors', 'atmospheric_correction', 'polymer_root'),
    enable_ac=('processors', 'atmospheric_correction', 'enable_ac'),
    mode_ac=('processors', 'atmospheric_correction', 'mode_ac'),
    auto_download_ecmwf=('processors', 'atmospheric_correction', 'auto_download_ecmwf'),
    scale_factor_boa_ref=('processors', 'atmospheric_correction', 'scale_factor_boa_ref'),
    threads=('processors', 'atmospheric_correction', 'threads'),
    blocksize=('processors', 'atmospheric_correction', 'blocksize'),

    # processors > smile
    run_smile_P=('processors', 'smile', 'run_processor'),

    # processors > dead_pixel
    run_deadpix_P=('processors', 'dead_pixel', 'run_processor'),
    deadpix_P_algorithm=('processors', 'dead_pixel', 'algorithm'),
    deadpix_P_interp_spectral=('processors', 'dead_pixel', 'interp_method_spectral'),
    deadpix_P_interp_spatial=('processors', 'dead_pixel', 'interp_method_spatial'),

    # processors > orthorectification
    ortho_resampAlg=('processors', 'orthorectification', 'resamp_alg'),
    vswir_overlap_algorithm=('processors', 'orthorectification', 'vswir_overlap_algorithm'),
    target_projection_type=('processors', 'orthorectification', 'target_projection_type'),
    target_epsg=('processors', 'orthorectification', 'target_epsg'),
    target_coord_grid=('processors', 'orthorectification', 'target_coord_grid'),
)


def get_updated_schema(source_schema, key2update, new_value):
    def deep_update(schema, key2upd, new_val):
        """Return true if update, else false."""
        for key in schema:
            if key == key2upd:
                schema[key] = new_val
            elif isinstance(schema[key], dict):
                deep_update(schema[key], key2upd, new_val)

        return schema

    from copy import deepcopy
    tgt_schema = deepcopy(source_schema)

    tgt_schema['processors']['schema']['orthorectification']['schema']['target_coord_grid'] = \
        dict(type='dict', required=False, nullable=True)

    return deep_update(tgt_schema, key2update, new_value)


enpt_schema_config_output = get_updated_schema(enpt_schema_input, key2update='required', new_value=True)


def get_param_from_json_config(paramname, json_config):
    keymap = parameter_mapping[paramname]  # tuple

    dict2search = json_config
    for i, k in enumerate(keymap):
        if i < len(keymap) - 1:
            # not the last element of the tuple -> contains a sub-dictionary
            dict2search = dict2search[k]
        elif isinstance(k, list):
            return [dict2search[sk] for sk in k]
        else:
            return dict2search[k]
