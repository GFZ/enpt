# -*- coding: utf-8 -*-

from ..model.images import EnMAP_L1B
from ..processors.radiometric_transform import TOARad2TOARef_Transformer
from geomultisens.processing.process_controller import process_controller




class EnPT_controller(process_controller):

    def __init__(self, entity_ID):
        # type: (str) -> None

        """EnPT Process Controller

        :param entity_ID:   ID of the scene to be processed, e.g. 'EM_20170515_12345'
        """

        raise NotImplementedError('The process controller is not yet working.')

        self.Images = [EnMAP_L1B, ]


    def run_all_processors(self):
        """
        Run all processors at once.
        """

        self.run_toaRad2toaRef()
        self.run_atmospheric_correction()


    def run_toaRad2toaRef(self):
        """
        Run conversion from TOA radiance to TOA reflectance.
        """

        # get a new instance of radiometric transformer
        RT = TOARad2TOARef_Transformer(None, None)

        # run transformation to TOARef
        RT.transform(self.Images[0])


    def run_atmospheric_correction(self):
        """
        Run atmospheric correction only.
        """

        pass