# -*- coding: utf-8 -*-
"""EnPT process controller module."""

from ..model.images import EnMAPL1Product_ImGeo
from ..processors.radiometric_transform import TOARad2TOARef_Transformer


class EnPT_controller(object):

    """Class of EnPT process controller"""

    def __init__(self, entity_ID):
        # type: (str) -> None
        """Initialize the Process Controller

        :param entity_ID:   ID of the scene to be processed, e.g. 'EM_20170515_12345'
        """
        raise NotImplementedError('The process controller is not yet working.')

        self.Images = [EnMAPL1Product_ImGeo, ]

    def run_all_processors(self):
        """Run all processors at once."""
        self.run_toaRad2toaRef()
        self.run_atmospheric_correction()

    def run_toaRad2toaRef(self):
        """Run conversion from TOA radiance to TOA reflectance."""
        # get a new instance of radiometric transformer
        RT = TOARad2TOARef_Transformer(None, None)

        # run transformation to TOARef
        RT.transform(self.Images[0])

    def run_atmospheric_correction(self):
        """Run atmospheric correction only."""
        pass
