# -*- coding: utf-8 -*-

from ..processors.radiometric_transform import TOARad2TOARef_Transformer
from GeoMultiSens_dev.processing.process_controller import process_controller




class EnPT_controller(process_controller):

    def __init__(self, entity_ID):
        # type: (str) -> self

        """EnPT Process Controller

        :param entity_ID:   ID of the scene to be processed, e.g. 'EM_20170515_12345'
        """

        pass


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
        RT = TOARad2TOARef_Transformer()

        # run transformation to TOARef
        RT.transform()


    def run_atmospheric_correction(self):
        """
        Run atmospheric correction only.
        """

        pass