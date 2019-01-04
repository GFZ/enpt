# -*- coding: utf-8 -*-
"""EnPT 'radiometric transform' module.

Contained Transformations:
    - TOA radiance to TOA reflectance
"""

import math
import numpy as np

from ...model.images import EnMAPL1Product_SensorGeo, EnMAP_Detector_SensorGeo  # noqa F401  # flake8 issue
from ...options.config import EnPTConfig


class Radiometric_Transformer(object):
    """Class for performing all kinds of radiometric transformations of EnMAP images."""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of Radiometric_Transformer."""
        self.cfg = config
        self.solarIrr = config.path_solar_irr  # path of model for solar irradiance
        self.earthSunDist = config.path_earthSunDist  # path of model for earth sun distance

    def transform_TOARad2TOARef(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        """Transform top-of-atmosphere radiance to top-of-atmosphere reflectance.

        NOTE: The following formula is used:
                toaRef = (scale_factor * math.pi * toaRad * earthSunDist**2) /
                         (solIrr * math.cos(zenithAngleDeg))

        :param enmap_ImageL1:   instance of the class 'EnMAPL1Product_ImGeo'
        :return:
        """
        for detectorName in enmap_ImageL1.detector_attrNames:
            detector = getattr(enmap_ImageL1, detectorName)  # type: EnMAP_Detector_SensorGeo

            enmap_ImageL1.logger.info('Converting TOA radiance to TOA reflectance for %s detector...'
                                      % detector.detector_name)

            # compute TOA reflectance
            constant = \
                self.cfg.scale_factor_toa_ref * math.pi * enmap_ImageL1.meta.earthSunDist ** 2 / \
                (math.cos(math.radians(enmap_ImageL1.meta.geom_sun_zenith)))
            solIrr = np.array([detector.detector_meta.solar_irrad[band] for band in detector.detector_meta.srf.bands])\
                .reshape(1, 1, detector.data.bands)
            toaRef = (constant * detector.data[:] / solIrr).astype(np.int16)

            # update EnMAP image
            detector.data = toaRef
            detector.detector_meta.unit = '0-%d' % self.cfg.scale_factor_toa_ref
            detector.detector_meta.unitcode = 'TOARef'

        return enmap_ImageL1
