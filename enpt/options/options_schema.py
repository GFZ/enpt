"""Definition of EnPT options schema (as used by cerberus library)."""
enpt_schema_input = dict(

    general_opts=dict(
        type='dict', required=False,
        schema=dict(
            CPUs=dict(type='integer', required=False, nullable=True),
            log_level=dict(type='string', required=False, allowed=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
            path_l1b_enmap_image=dict(type='string', required=False),
            path_l1b_enmap_image_gapfill=dict(type='string', required=False),
            path_l1b_snr_model=dict(type='string', required=False),
        )),

    processors=dict(
        type='dict', required=False,
        schema=dict(

            toa_ref=dict(
                type='dict', required=False,
                schema=dict(
                    path_earthSunDist=dict(type='string', required=False),
                    path_solar_irr=dict(type='string', required=False),
                )),

            geometry=dict(
                type='dict', required=False,
                schema=dict(
                    enable_keystone_correction=dict(type='boolean', required=False),
                    enable_vnir_swir_coreg=dict(type='boolean', required=False),
                    path_reference_image=dict(type='string', required=False),
                )),

            atmospheric_correction=dict(
                type='dict', required=False,
                schema=dict(
                    auto_download_ecmwf=dict(type='boolean', required=False),
                    enable_cloud_screening=dict(type='boolean', required=False),
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
                )),

            orthorectification=dict(
                type='dict', required=False,
                schema=dict(
                    resamp_alg=dict(type='integer', required=False),
                ))
        ))
)


parameter_mapping = dict(
    # general opts
    CPUs=('general_opts', 'CPUs'),
    log_level=('general_opts', 'log_level'),
    path_l1b_enmap_image=('general_opts', 'path_l1b_enmap_image'),
    path_l1b_enmap_image_gapfill=('general_opts', 'path_l1b_enmap_image_gapfill'),
    path_l1b_snr_model=('general_opts', 'path_l1b_snr_model'),

    # processors > toa_ref
    path_earthSunDist=('processors', 'toa_ref', 'path_earthSunDist'),
    path_solar_irr=('processors', 'toa_ref', 'path_solar_irr'),

    # processors > geometry
    enable_keystone_correction=('processors', 'geometry', 'enable_keystone_correction'),
    enable_vnir_swir_coreg=('processors', 'geometry', 'enable_vnir_swir_coreg'),
    path_reference_image=('processors', 'geometry', 'path_reference_image'),

    # processors > atmospheric_correction
    auto_download_ecmwf=('processors', 'atmospheric_correction', 'auto_download_ecmwf'),
    enable_cloud_screening=('processors', 'atmospheric_correction', 'enable_cloud_screening'),

    # processors > smile
    run_smile_P=('processors', 'smile', 'run_processor'),

    # processors > dead_pixel
    run_deadpix_P=('processors', 'dead_pixel', 'run_processor'),

    # processors > orthorectification
    ortho_resampAlg=('processors', 'orthorectification', 'resamp_alg'),
)


def get_updated_schema(source_schema, key2update, new_value):
    def deep_update(schema, key2upd, new_val):
        """Return true if update, else false"""

        for key in schema:
            if key == key2upd:
                schema[key] = new_val
            elif isinstance(schema[key], dict):
                deep_update(schema[key], key2upd, new_val)

        return schema

    from copy import deepcopy
    tgt_schema = deepcopy(source_schema)
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
