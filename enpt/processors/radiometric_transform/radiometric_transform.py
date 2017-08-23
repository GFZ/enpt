# -*- coding: utf-8 -*-

"""
transforms TOA radiance to TOA reflectance
"""


from enpt.model.images import EnMAPL1Product_SensorGeo


class TOARad2TOARef_Transformer(object):

    def __init__(self, solarIrr, sunEarthDist):
        # type: (dict, dict) -> None
        """Class for all kinds of radiometric transformations.

        :param solarIrr:        model for solar irradiance
        :param sunEarthDist:    model for sun earth distance
        """

        self.solarIrr = solarIrr
        self.sunEarthDist = sunEarthDist

    def transform(self, enmap_ImageL1: EnMAPL1Product_SensorGeo):
        """Transform top-of-atmosphere radiance to top-of-atmosphere reflectance.

        :param enmap_ImageL1: instance of the class 'EnMAPL1Product_ImGeo'
        :return:
        """
        # TODO: Implement that function! The code below is only dummy code.
        enmap_ImageL1.logger.info('I am the TOARef transformer and I am logging something! \n')

        # get original data for VNIR
        rows, cols, _ = enmap_ImageL1.vnir.data.shape
        original_data = enmap_ImageL1.vnir.data[int(rows / 2):int((rows / 2) + 5),
                                                int(cols / 2):int((cols / 2) + 5), 0]

        # do something
        result = original_data * 2

        # give some output
        enmap_ImageL1.logger.info('This is a small subset of the input image: \n%s \n' % repr(original_data))
        enmap_ImageL1.logger.info('This is what comes out after multiplying it with 2: \n%s \n' % repr(result))

        # merge the results into the EnMAP image
        enmap_ImageL1.vnir.data = result

        return enmap_ImageL1
