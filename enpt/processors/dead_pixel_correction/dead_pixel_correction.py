# -*- coding: utf-8 -*-
"""EnPT 'dead pixel correction' module.

Performs the Dead Pixel Correction using a linear interpolation in spectral dimension.
"""
from typing import Union
from numbers import Number  # noqa: F401
import logging
from tqdm import tqdm

import numpy as np
from scipy.interpolate import griddata, interp1d
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

            # NOTE: this is done first, because spectral interpolation needs the information of the outermost pixels

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

    def _interpolate_nodata_spectrally(self, image2correct: GeoArray, deadpixel_map: GeoArray, progress=False):
        assert deadpixel_map.ndim == 3, "3D dead pixel map expected."
        if deadpixel_map.shape != image2correct.shape:
            raise ValueError("Dead pixel map and image to be corrected must have equal shape.")

        image_corrected = interp_nodata_along_axis(image2correct, axis=2, nodata=deadpixel_map[:],
                                                   method=self.interp_alg, fill_value='extrapolate',
                                                   progress=progress)

        return image_corrected

    def _correct_using_3D_deadpixelmap(self,
                                       image2correct: GeoArray,
                                       deadpixel_map: GeoArray,
                                       progress=False):
        assert deadpixel_map.ndim == 3, "3D dead pixel map expected."
        if deadpixel_map.shape != image2correct.shape:
            raise ValueError("Dead pixel map and image to be corrected must have equal shape.")

        if self.algorithm == 'spectral':
            image_corrected = interp_nodata_along_axis(image2correct, axis=2, nodata=deadpixel_map[:],
                                                       method=self.interp_alg, fill_value='extrapolate')
        else:
            raise NotImplementedError("Currently only the algorithm 'spectral' is implemented.")

        return image_corrected

    def correct(self, image2correct: Union[np.ndarray, GeoArray], deadpixel_map: Union[np.ndarray, GeoArray],
                progress=False):
        image2correct = GeoArray(image2correct) if not isinstance(image2correct, GeoArray) else image2correct

        self._validate_inputs(image2correct, deadpixel_map)

        if 1 in list(np.unique(deadpixel_map)):
            if deadpixel_map.ndim == 2:
                deadcolumn_map = deadpixel_map

                # compute dead pixel percentage
                prop = np.any(deadcolumn_map, axis=0).sum() * image2correct.shape[0] / np.dot(*image2correct.shape[:2])

                # convert 2D deadcolumn_map to 3D deadpixel_map
                B, C = deadcolumn_map.shape
                deadpixel_map = np.empty((image2correct.shape[0], C, B), np.bool)
                deadpixel_map[:, :, :] = deadcolumn_map.T

            else:
                # compute dead pixel percentage
                prop = np.any(deadpixel_map, axis=2).sum() / np.dot(*image2correct.shape[:2])

            self.logger.info('Percentage of defective pixels: %.2f' % (prop * 100))

            # run correction
            if self.algorithm == 'spectral':
                return self._interpolate_nodata_spectrally(image2correct, deadpixel_map, progress=progress)
            else:
                raise NotImplementedError("Currently only the algorithm 'spectral' is implemented.")

        else:
            self.logger.info("Dead pixel correction skipped because dead pixel mask labels no pixels as 'defective'.")
            return image2correct


def interp_nodata_along_axis_2d(data_2d: np.ndarray, axis: int = 0,
                                nodata: Union[np.ndarray, Number] = np.nan,
                                method='linear', fill_value='extrapolate', progress=False):
    if data_2d.ndim != 2:
        raise ValueError('Expected a 2D array. Received a %dD array.' % data_2d.ndim)
    if axis > data_2d.ndim:
        raise ValueError("axis=%d is out of bounds for data with %d dimensions." % (axis, data_2d.ndim))

    data_2d = data_2d if axis == 1 else data_2d.T

    if isinstance(nodata, np.ndarray):
        badmask_full = nodata if axis == 1 else nodata.T

        if badmask_full.shape != data_2d.shape:
            raise ValueError('No-data mask and data must have the same shape.')

    else:
        badmask_full = ~np.isfinite(data_2d) if ~np.isfinite(nodata) else data_2d == nodata

    # call 1D interpolation vectorized
    #   => group the dataset by rows that have nodata at the same column position
    #   => remember the row positions, call the intpolation for these rows at once (vectorized)
    #      and substitute the original data  at the previously grouped row positions
    rowsNrs2interp = np.argwhere(np.any(badmask_full, axis=1)).flatten()
    badmask_sub2interp = badmask_full[rowsNrs2interp, :]

    for col in range(data_2d.shape[1]):
        rowNrs_with_nodata_at_col = rowsNrs2interp[badmask_sub2interp[:, col]]  # badmask (bool) -> col works as mask

        # TODO remove those rows that also have nodata in other cols

        data_2d_sub_with_nodata_at_col = data_2d[rowNrs_with_nodata_at_col, :]
        badpos = col
        goodpos = np.delete(np.arange(data_2d.shape[1]), badpos)

        data_2d_sub_with_nodata_at_col[:, badpos] = \
            interp1d(goodpos, data_2d_sub_with_nodata_at_col[:, goodpos],
                     axis=1, kind=method, fill_value=fill_value, bounds_error=False)(badpos)

        data_2d[rowNrs_with_nodata_at_col, :] = data_2d_sub_with_nodata_at_col

    return data_2d if axis == 1 else data_2d.T


def interp_nodata_along_axis(data, axis=0, nodata: Union[np.ndarray, Number] = np.nan,
                             method='linear', fill_value='extrapolate', progress=False):
    assert axis <= 2
    if data.ndim not in [2, 3]:
        raise ValueError('Expected a 2D or 3D array. Received a %dD array.' % data.ndim)
    if isinstance(nodata, np.ndarray) and nodata.shape != data.shape:
        raise ValueError('No-data mask and data must have the same shape.')

    kw = dict(method=method, fill_value=fill_value, progress=progress)

    if data.ndim == 2:
        return interp_nodata_along_axis_2d(data, axis=axis, nodata=nodata, **kw)

    else:
        def reshape_input(In):
            R, C, B = In.shape
            return \
                In.reshape(C, R * B) if axis == 0 else \
                np.transpose(In, axes=[1, 0, 2]).reshape(C, R * B).T if axis == 1 else \
                In.reshape(R * C, B)

        def reshape_output(out):
            return \
                out.reshape(data.shape) if axis in [0, 2] else \
                np.transpose(out.T.reshape(data.shape), axes=[1, 0, 2])

        return \
            reshape_output(
                interp_nodata_along_axis_2d(
                    data_2d=reshape_input(data),
                    nodata=reshape_input(nodata) if isinstance(nodata, np.ndarray) else nodata,
                    axis=axis if axis != 2 else 1, **kw))
