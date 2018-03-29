# -*- coding: utf-8 -*-
"""EnPT process controller module."""

from ..model.images import EnMAPL1Product_SensorGeo
from ..processors.radiometric_transform import TOARad2TOARef_Transformer
from ..options.config import EnPTConfig


class EnPT_controller(object):
    """Class of EnPT process controller."""

    def __init__(self, config: EnPTConfig=None, **config_kwargs: dict):
        # type: (str) -> None
        """Initialize the Process Controller.

        :param config:          an instance of the EnPTConfig class (overrides config_kwargs)
        :param config_kwargs:   configuration parameters to be passed to EnPTConfig class
        """
        self.cfg = config or EnPTConfig(**config_kwargs)

        self.Images = [EnMAPL1Product_SensorGeo, ]

    def run_all_processors(self):
        """Run all processors at once."""
        self.run_toaRad2toaRef()
        self.run_atmospheric_correction()

    def run_toaRad2toaRef(self):
        """Run conversion from TOA radiance to TOA reflectance."""
        # get a new instance of radiometric transformer
        RT = TOARad2TOARef_Transformer(None, None)

        # run transformation to TOARef
        RT.transform_dummy(self.Images[0])

    def run_atmospheric_correction(self):
        """Run atmospheric correction only."""
        pass
