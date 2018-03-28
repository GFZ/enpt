"""Definition of EnPT options schema (as used by cerberus library)."""
enpt_schema_input = dict(

    general_opts=dict(
        type='dict', required=False,
        schema=dict(
            CPUs=dict(type='integer', required=False, nullable=True),
            log_level=dict(type='string', required=False, allowed=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
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
                )),

            atmospheric_correction=dict(
                type='dict', required=False,
                schema=dict(
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
