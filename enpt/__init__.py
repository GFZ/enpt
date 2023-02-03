# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2023 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""EnMAP processing tool (EnPT) software package developed by GFZ."""

import sensormapgeo as __smg  # noqa (E402 + F401)  # only to avoid later import error due to static TLS
import os as __os

from .version import __version__, __versionalias__   # noqa (E402 + F401)
from .options.config import EnPTConfig
from .execution.controller import EnPT_Controller

__author__ = """Daniel Scheffler"""
__email__ = 'danschef@gfz-potsdam.de'
__all__ = ['__version__',
           '__versionalias__',
           '__author__',
           '__email__',
           'EnPTConfig',
           'EnPT_Controller'
           ]

# $PROJ_LIB was renamed to $PROJ_DATA in proj=9.1.1, which leads to issues with fiona>=1.8.20,<1.9
# https://github.com/conda-forge/pyproj-feedstock/issues/130
# -> fix it by setting PROJ_DATA
if 'GDAL_DATA' in __os.environ and 'PROJ_DATA' not in __os.environ and 'PROJ_LIB' not in __os.environ:
    __os.environ['PROJ_DATA'] = __os.path.join(__os.path.dirname(__os.environ['GDAL_DATA']), 'proj')
