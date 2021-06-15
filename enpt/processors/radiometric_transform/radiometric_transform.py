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

"""EnPT 'radiometric transform' module.

Contained Transformations:
    - TOA radiance to TOA reflectance
"""

import math
import numpy as np

from ...model.images import EnMAPL1Product_SensorGeo, EnMAP_Detector_SensorGeo  # noqa F401  # flake8 issue
from ...options.config import EnPTConfig

__author__ = 'Daniel Scheffler'


class Radiometric_Transformer(object):
    """Class for performing all kinds of radiometric transformations of EnMAP images."""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of Radiometric_Transformer.

        :param config: EnPT configuration object
        """
        self.cfg = config
        self.solarIrr = config.path_solar_irr  # path of model for solar irradiance
        self.earthSunDist = config.path_earthSunDist  # path of model for earth sun distance

    def transform_TOARad2TOARef(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        """Transform top-of-atmosphere radiance to top-of-atmosphere reflectance (16-bit signed-integer).

        NOTE: The following formula is used:
                toaRef = (scale_factor * math.pi * toaRad * earthSunDist**2) /
                         (solIrr * math.cos(zenithAngleDeg))

        :param enmap_ImageL1:   instance of the class 'EnMAPL1Product_ImGeo'
        :return:
        """
        for detectorName in enmap_ImageL1.detector_attrNames:
            detector: EnMAP_Detector_SensorGeo = getattr(enmap_ImageL1, detectorName)

            enmap_ImageL1.logger.info('Converting TOA radiance to TOA reflectance for %s detector...'
                                      % detector.detector_name)

            # compute TOA reflectance
            constant = \
                self.cfg.scale_factor_toa_ref * math.pi * enmap_ImageL1.meta.earthSunDist ** 2 / \
                (math.cos(math.radians(enmap_ImageL1.meta.geom_sun_zenith)))
            solIrr = detector.detector_meta.solar_irrad.reshape(1, 1, detector.data.bands)
            toaRef = (constant * detector.data[:] / solIrr).astype(np.int16)

            # update EnMAP image
            detector.data = toaRef
            detector.detector_meta.unit = '0-%d' % self.cfg.scale_factor_toa_ref
            detector.detector_meta.unitcode = 'TOARef'

        return enmap_ImageL1
