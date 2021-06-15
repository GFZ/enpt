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

"""EnPT metadata modules. All objects and functions regarding EnMAP metadata are implemented here."""

from .metadata_sensorgeo import EnMAP_Metadata_L1B_Detector_SensorGeo, EnMAP_Metadata_L1B_SensorGeo  # noqa: F401
from .metadata_mapgeo import EnMAP_Metadata_L2A_MapGeo  # noqa: F401


__author__ = 'Daniel Scheffler'


# Define L1B_product_props
L1B_product_props = dict(
    xml_detector_label=dict(
        VNIR='VNIRDetector',
        SWIR='SWIRDetector'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)


L1B_product_props_DLR = dict(
    xml_detector_label=dict(
        VNIR='vnir',
        SWIR='swir'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)


# Define L1B_product_props
L2A_product_props_DLR = dict(
    xml_detector_label=dict(
        VNIR='vnir',
        SWIR='swir'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)
