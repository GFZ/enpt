Usage
=====

Usage of the Python API
***********************

To start run whole EnPT processing pipeline via the Python API::

    from enpt.execution.controller import EnPT_Controller

    config_minimal = dict(
        path_l1b_enmap_image='/path/ENMAP*L1B*.zip',
        path_dem='/path/to/overlapping/DEM.bsq'
    )
    CTR = EnPT_Controller(**config_minimal)
    CTR.run_all_processors()

Further configuration parameters are documented here_.
Note that the class `EnPTConfig` takes the same keyword arguments like the `EnPT_Controller` class.

.. _here: http://enmap.gitext.gfz-potsdam.de/GFZ_Tools_EnMAP_BOX/EnPT/doc/enpt.options.html#enpt.options.config.EnPTConfig

You can also pass a JSON-File with your EnPT configuration to the `EnPT_Controller` class. This allows you to easily
copy and reuse configuration files. A template with all possible options and defaults can be found in
`enpt/options/options_default.json`_.

The corresponding Python call looks like this::

    from enpt.execution.controller import EnPT_Controller
    from enpt.options.config import EnPTConfig

    CFG = EnPTConfig(json_config='/path/to/your/config_file.json')
    CTR = EnPT_Controller(config=CFG)
    CTR.run_all_processors()

.. _enpt/options/options_default.json: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/EnPT/blob/master/enpt/options/options_default.json


Command line utilities
**********************

enpt_cli.py
-----------

At the command line, EnPT provides the **enpt_cli.py** command:

.. argparse::
   :filename: ./../bin/enpt_cli.py
   :func: get_enpt_argparser
   :prog: enpt_cli.py


QGIS GUI
********

There is a separate graphical user interface (GUI) for EnPT than can be installed as an EnMAP-Box application in QGIS.
To install it in QGIS, please refer to the separate repository enpt_enmapboxapp_.

Here is screenshot of the current version:

.. image:: img/screenshot_enpt_enmapboxapp_874x1047.png

.. _enpt_enmapboxapp: https://gitext.gfz-potsdam.de/EnMAP/GFZ_Tools_EnMAP_BOX/enpt_enmapboxapp