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
    def __init__(self, algorithm: str = 'spectral', interp: str = 'linear', filter_halfwidth: int = 1,
                 logger: logging.Logger = None):
        """Get an instance of Dead_Pixel_Corrector.

        :param algorithm:           algorithm how to correct dead pixels
                                    'spectral': interpolate in the spectral domain
                                    'spatial':  interpolate in the spatial domain
        :param interp:              interpolation algorithm ('linear', 'bilinear', 'cubic', 'spline')
        :param filter_halfwidth:    half width of interpolation filter
                                    (determines the number of adjacant pixels to be respected in interpolation)
        :param logger:
        """
        self.algorithm = algorithm
        self.interp_alg = interp
        self.fhw = filter_halfwidth
        self.logger = logger or logging.getLogger()

    @staticmethod
    def _validate_inputs(image2correct: GeoArray, deadpixel_map: GeoArray):
        if deadpixel_map.ndim == 2:
            if (image2correct.bands, image2correct.columns) != deadpixel_map.shape:
                raise ValueError('The given image to be corrected (shape: %s) requires a dead column map with shape '
                                 '(%s, %s). Received %s.'
                                 % (image2correct.shape, image2correct.bands,
                                    image2correct.columns, deadpixel_map.shape))
        elif deadpixel_map.ndim == 3:
            if image2correct.shape != deadpixel_map.shape:
                raise ValueError('The given image to be corrected (shape: %s) requires a dead pixel map with equal '
                                 'shape. Received %s.' % (image2correct.shape, deadpixel_map.shape))
        else:
            raise ValueError('Unexpected number of dimensions of dead column map.')

    def _correct_using_2D_deadpixelmap(self,
                                       image2correct: GeoArray,
                                       deadcolumn_map: GeoArray,
                                       progress=False):
        # TODO speed this up
        """

        NOTE: dead columns in the first or the last band are unmodified.

        :param image2correct:
        :param deadcolumn_map:
        :param progress:        whether to show progress bars
        :return:
        """
        assert deadcolumn_map.ndim == 2, "2D dead column map expected."

        #################
        # Interpolation #
        #################

        if self.algorithm == 'spectral':
            # set bands where no spectral interpolation is possible -> spatial interpolation
            band_nbrs_spatial_interp = \
                list(range(self.fhw)) + list(range(image2correct.bands - 1, image2correct.bands - self.fhw - 1, -1))

            # spatial interpolation (fallback) #
            ####################################

            # NOTE: this is done first, because spectral interpolation needs the information on the outermost pixels

            for band, column in np.argwhere(deadcolumn_map):

                # only first or last bands (number depends on filter half width)
                if band in band_nbrs_spatial_interp:
                    if column in [0, image2correct.shape[1] - 1]:
                        # currently, dead pixels in the outermost bands at the outermost columns are not corrected
                        self.logger.warning('Unable to correct dead column %s in band %s.' % (column, band))
                    else:
                        self.logger.debug('Correcting dead column %s in band %s.' % (column, band))

                        band_data_orig = image2correct[:, column - 1:column + 2, band]  # target and adjacent columns
                        band_data_float = band_data_orig.astype(np.float)
                        band_data_float[:, 1] = np.NaN

                        x, y = np.indices(band_data_float.shape)
                        interp = np.array(band_data_float)

                        interp[np.isnan(interp)] = \
                            griddata(np.array([x[~np.isnan(band_data_float)],
                                               y[~np.isnan(band_data_float)]]).T,  # points we know
                                     band_data_float[~np.isnan(band_data_float)],  # values we know
                                     np.array([x[np.isnan(band_data_float)],
                                               y[np.isnan(band_data_float)]]).T,  # points to interpolate
                                     method=self.interp_alg)

                        # copy corrected columns to image2correct
                        interp[np.isnan(interp)] = band_data_orig[np.isnan(interp)]
                        image2correct[:, column - 1:column + 2, band] = interp.astype(image2correct.dtype)

            # spectral interpolation #
            #########################

            for band, column in tqdm(np.argwhere(deadcolumn_map), disable=False if progress else True):
                if band in band_nbrs_spatial_interp:
                    continue  # already interpolated spatially above

                # any other band
                else:
                    self.logger.debug('Correcting dead column %s in band %s.' % (column, band))

                    column_data_orig = image2correct[:, column, band - self.fhw:band + self.fhw + 1]
                    column_data_float = column_data_orig.astype(np.float)
                    column_data_float[:, self.fhw] = np.NaN

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
                    image2correct[:, column, band - self.fhw:band + self.fhw + 1] = interp.astype(image2correct.dtype)

        else:
            raise NotImplementedError("Currently only the algorithm 'spectral' is implemented.")

        return image2correct

    def _correct_using_3D_deadpixelmap(self,
                                       image2correct: GeoArray,
                                       deadpixel_map: GeoArray,
                                       progress=False):
        # TODO implement that
        assert deadpixel_map.ndim == 3, "3D dead pixel map expected."

        self.logger.warning('Dead pixel correction skipped. Currently a dead pixel correction based on a 3D dead pixel '
                            'map has not been implemented.')
        return image2correct

    def correct(self, image2correct: Union[np.ndarray, GeoArray], deadpixel_map: Union[np.ndarray, GeoArray],
                progress=False):
        image2correct = GeoArray(image2correct) if not isinstance(image2correct, GeoArray) else image2correct

        self._validate_inputs(image2correct, deadpixel_map)

        if 1 in list(np.unique(deadpixel_map)):
            if deadpixel_map.ndim == 2:
                return self._correct_using_2D_deadpixelmap(image2correct, deadpixel_map, progress)
            else:
                return self._correct_using_3D_deadpixelmap(image2correct, deadpixel_map, progress)

        else:
            self.logger.info("Dead pixel correction skipped because dead pixel mask labels no pixels as 'defective'.")
            return image2correct
