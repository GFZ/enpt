# -*- coding: utf-8 -*-
"""EnPT 'dead pixel correction' module.

Performs the Dead Pixel Correction using a linear interpolation in spectral dimension.
"""
from typing import Union
import logging
from tqdm import tqdm

import numpy as np
from scipy.interpolate import griddata
from geoarray import GeoArray


class Dead_Pixel_Corrector(object):
    def __init__(self, algorithm: str='spectral', interp: str='linear', filter_halfwidth: int=1,
                 logger: logging.Logger=None):
        """Get an instance of Dead_Pixel_Corrector.

        :param algorithm:           algorithm how to correct dead pixels
                                    'spectral': interpolate in the spectral domain
                                    'spatial':  interpolate in the spatial domain
        :param interp:              interpolation algorithm ('linear', 'bilinear', 'cubic', 'spline')
        :param filter_halfwidth:    half width of interpolation filter
                                    (determines the number of adjacant pixels to be respected in interpation)
        :param logger:
        """
        self.algorithm = algorithm
        self.interp_alg = interp
        self.fhw = filter_halfwidth
        self.logger = logger or logging.getLogger()

    @staticmethod
    def _validate_inputs(image2correct: GeoArray, deadcolumn_map: GeoArray):
        if (image2correct.bands, image2correct.columns) != deadcolumn_map.shape:
            raise ValueError('The given image to be corrected (shape: %s) requires a dead column map with shape '
                             '(%s, %s). Received %s.'
                             % (image2correct.shape, image2correct.bands, image2correct.columns, deadcolumn_map.shape))

    def correct(self, image2correct: Union[np.ndarray, GeoArray], deadcolumn_map: Union[np.ndarray, GeoArray]):
        """

        NOTE: dead columns in the first or the last band are unmodified.

        :param image2correct:
        :param deadcolumn_map:
        :return:
        """
        image2correct = GeoArray(image2correct) if not isinstance(image2correct, GeoArray) else image2correct
        deadcolumn_map = GeoArray(deadcolumn_map) if not isinstance(deadcolumn_map, GeoArray) else deadcolumn_map

        # validate
        self._validate_inputs(image2correct, deadcolumn_map)

        # correct
        if self.algorithm == 'spectral':

            for band, column in tqdm(np.argwhere(deadcolumn_map), disable=False):  # set disable=True to mute progress
                self.logger.debug('Correcting dead column %s in band %s.' % (column, band))

                #################
                # Interpolation #
                #################

                # first or last bands (number depends on filter half width)
                if band in \
                        list(range(self.fhw)) + list(range(image2correct.bands-1, image2correct.bands-self.fhw-1, -1)):
                    pass  # TODO

                # any other band
                else:
                    column_data_orig = image2correct[:, column, band-self.fhw:band+self.fhw+1]
                    column_data_float = column_data_orig.astype(np.float)
                    column_data_float[:, 1] = np.NaN

                    x, y = np.indices(column_data_float.shape)
                    interp = np.array(column_data_float)

                    interp[np.isnan(interp)] = \
                        griddata(np.array([x[~np.isnan(column_data_float)],
                                           y[~np.isnan(column_data_float)]]).T,  # points we know
                                 column_data_float[~np.isnan(column_data_float)],  # values we know
                                 np.array([x[np.isnan(column_data_float)],
                                           y[np.isnan(column_data_float)]]).T,  # points to interpolate
                                 method=self.interp_alg)

                    # copy corrected columns to image2correct
                    interp[np.isnan(interp)] = column_data_orig[np.isnan(interp)]
                    image2correct[:, column, band-self.fhw:band+self.fhw+1] = interp.astype(image2correct.dtype)

        else:
            raise NotImplementedError("Currently only the algorithm 'spectral' is implemented.")

        return image2correct
