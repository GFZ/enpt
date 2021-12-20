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

"""EnPT 'dead pixel correction' module.

Performs the Dead Pixel Correction using a linear interpolation in spectral dimension.
"""
from typing import Union
from numbers import Number  # noqa: F401
import logging

import numpy as np
import numpy_indexed as npi
from multiprocessing import Pool, cpu_count
from scipy.interpolate import griddata, interp1d
from pandas import DataFrame
from geoarray import GeoArray

__author__ = 'Daniel Scheffler'


class Dead_Pixel_Corrector(object):
    """EnPT Dead Pixel Correction class.

    The EnPT dead pixel correction uses the pixel masks provided by DLR and interpolates the EnMAP image
    data at the indicated dead pixel positions. It supports two interpolation algorithms:

    1. spectral interpolation
        * Interpolates the data in the spectral domain.
        * Points outside the data range are extrapolated.
    2. spatial interpolation
        * Interpolates the data spatially.
        * Remaining missing data positions (e.g., outermost columns) are spectrally interpolated.
    """

    def __init__(self,
                 algorithm: str = 'spectral',
                 interp_spectral: str = 'linear',
                 interp_spatial: str = 'linear',
                 CPUs: int = None,
                 logger: logging.Logger = None):
        """Get an instance of Dead_Pixel_Corrector.

        :param algorithm:           algorithm how to correct dead pixels
                                    'spectral': interpolate in the spectral domain
                                    'spatial':  interpolate in the spatial domain
        :param interp_spectral:     spectral interpolation algorithm
                                    (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, etc.)
        :param interp_spatial:      spatial interpolation algorithm ('linear', 'bilinear', 'cubic', 'spline')
        :param CPUs:                number of CPUs to use for interpolation (only relevant if algorithm = 'spatial')
        :param logger:
        """
        self.algorithm = algorithm
        self.interp_alg_spectral = interp_spectral
        self.interp_alg_spatial = interp_spatial
        self.CPUs = CPUs or cpu_count()
        self.logger = logger or logging.getLogger()

    @staticmethod
    def _validate_inputs(image2correct: GeoArray,
                         deadpixel_map: GeoArray):
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

    def _interpolate_nodata_spectrally(self,
                                       image2correct: GeoArray,
                                       deadpixel_map: GeoArray):
        assert deadpixel_map.ndim == 3, "3D dead pixel map expected."
        if deadpixel_map.shape != image2correct.shape:
            raise ValueError("Dead pixel map and image to be corrected must have equal shape.")

        image_corrected = interp_nodata_along_axis(image2correct, axis=2, nodata=deadpixel_map[:],
                                                   method=self.interp_alg_spectral, fill_value='extrapolate')

        return image_corrected

    def _interpolate_nodata_spatially(self,
                                      image2correct: GeoArray,
                                      deadpixel_map: GeoArray):
        assert deadpixel_map.ndim == 3, "3D dead pixel map expected."
        if deadpixel_map.shape != image2correct.shape:
            raise ValueError("Dead pixel map and image to be corrected must have equal shape.")

        band_indices_with_nodata = np.argwhere(np.any(np.any(deadpixel_map, axis=0), axis=0)).flatten()
        image_sub = image2correct[:, :, band_indices_with_nodata]
        deadpixel_map_sub = deadpixel_map[:, :, band_indices_with_nodata]

        kw = dict(method=self.interp_alg_spatial, fill_value=np.nan, implementation='pandas', CPUs=self.CPUs)

        # correct dead columns
        image_sub_interp = interp_nodata_spatially_3d(image_sub, axis=1, nodata=deadpixel_map_sub, **kw)

        # correct dead rows
        if np.isnan(image_sub_interp).any():
            image_sub_interp = interp_nodata_spatially_3d(image_sub_interp, axis=0,
                                                          nodata=np.isnan(image_sub_interp), **kw)

        image2correct[:, :, band_indices_with_nodata] = image_sub_interp

        # correct remaining nodata by spectral interpolation (e.g., outermost columns)
        if np.isnan(image2correct).any():
            image2correct = interp_nodata_along_axis(image2correct, axis=2, nodata=np.isnan(image2correct),
                                                     method=self.interp_alg_spectral, fill_value='extrapolate')

        return image2correct

    def correct(self,
                image2correct: Union[np.ndarray, GeoArray],
                deadpixel_map: Union[np.ndarray, GeoArray]):
        """Run the dead pixel correction.

        :param image2correct:   image to correct
        :param deadpixel_map:   dead pixel map (2D (bands x columns) or 3D (rows x columns x bands)
        :return:    corrected image
        """
        image2correct = GeoArray(image2correct) if not isinstance(image2correct, GeoArray) else image2correct

        self._validate_inputs(image2correct, deadpixel_map)

        if 1 in list(np.unique(deadpixel_map)):
            if deadpixel_map.ndim == 2:
                deadcolumn_map = deadpixel_map

                # compute dead pixel percentage
                prop_dp_anyband = \
                    np.any(deadcolumn_map, axis=0).sum() * image2correct.shape[0] / np.dot(*image2correct.shape[:2])
                prop_dp = deadcolumn_map.sum() * image2correct.shape[0] / image2correct.size

                # convert 2D deadcolumn_map to 3D deadpixel_map
                B, C = deadcolumn_map.shape
                deadpixel_map = np.empty((image2correct.shape[0], C, B), bool)
                deadpixel_map[:, :, :] = deadcolumn_map.T

            else:
                # compute dead pixel percentage
                prop_dp_anyband = np.any(deadpixel_map, axis=2).sum() / np.dot(*image2correct.shape[:2])
                prop_dp = deadpixel_map.sum() / image2correct.size

            self.logger.info('Percentage of defective pixels: %.2f' % (prop_dp * 100))
            self.logger.debug('Percentage of pixels with a defect in any band: %.2f' % (prop_dp_anyband * 100))

            # run correction
            if self.algorithm == 'spectral':
                return self._interpolate_nodata_spectrally(image2correct, deadpixel_map)
            else:
                return self._interpolate_nodata_spatially(image2correct, deadpixel_map)

        else:
            self.logger.info("Dead pixel correction skipped because dead pixel mask labels no pixels as 'defective'.")
            return image2correct


def _get_baddata_mask(data: np.ndarray,
                      nodata: Union[np.ndarray, Number] = np.nan,
                      transpose_inNodata: bool = False):
    if isinstance(nodata, np.ndarray):
        badmask = nodata.T if transpose_inNodata else nodata

        if badmask.shape != data.shape:
            raise ValueError('No-data mask and data must have the same shape.')

    else:
        badmask = ~np.isfinite(data) if ~np.isfinite(nodata) else data == nodata

    return badmask


def interp_nodata_along_axis_2d(data_2d: np.ndarray,
                                axis: int = 0,
                                nodata: Union[np.ndarray, Number] = np.nan,
                                method: str = 'linear',
                                fill_value: Union[float, str] = 'extrapolate'):
    """Interpolate a 2D array along the given axis (based on scipy.interpolate.interp1d).

    :param data_2d:         data to interpolate
    :param axis:            axis to interpolate (0: along columns; 1: along rows)
    :param nodata:          nodata array in the shape of data or nodata value
    :param method:          interpolation method (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, etc.)
    :param fill_value:      value to fill into positions where no interpolation is possible
                            - if 'extrapolate': extrapolate the missing values
    :return:    interpolated array
    """
    if data_2d.ndim != 2:
        raise ValueError('Expected a 2D array. Received a %dD array.' % data_2d.ndim)
    if axis > data_2d.ndim:
        raise ValueError("axis=%d is out of bounds for data with %d dimensions." % (axis, data_2d.ndim))

    data_2d = data_2d if axis == 1 else data_2d.T

    badmask_full = _get_baddata_mask(data_2d, nodata, transpose_inNodata=axis == 0)

    # call 1D interpolation vectorized
    #   => group the dataset by rows that have nodata at the same column position
    #   => remember the row positions, call the interpolation for these rows at once (vectorized)
    #      and substitute the original data  at the previously grouped row positions
    groups_unique_rows = npi.group_by(badmask_full).split(np.arange(len(badmask_full)))

    for indices_unique_rows in groups_unique_rows:
        badmask_grouped_rows = badmask_full[indices_unique_rows, :]

        if np.any(badmask_grouped_rows[0, :]):
            badpos = np.argwhere(badmask_grouped_rows[0, :]).flatten()
            goodpos = np.delete(np.arange(data_2d.shape[1]), badpos)

            if goodpos.size > 1:
                data_2d_grouped_rows = data_2d[indices_unique_rows]

                data_2d_grouped_rows[:, badpos] = \
                    interp1d(goodpos, data_2d_grouped_rows[:, goodpos],
                             axis=1, kind=method, fill_value=fill_value, bounds_error=False)(badpos)

                data_2d[indices_unique_rows, :] = data_2d_grouped_rows

    return data_2d if axis == 1 else data_2d.T


def interp_nodata_along_axis(data,
                             axis=0,
                             nodata: Union[np.ndarray, Number] = np.nan,
                             method: str = 'linear',
                             fill_value: Union[float, str] = 'extrapolate'):
    """Interpolate a 2D or 3D array along the given axis (based on scipy.interpolate.interp1d).

    :param data:            data to interpolate
    :param axis:            axis to interpolate (0: along columns; 1: along rows, 2: along bands)
    :param nodata:          nodata array in the shape of data or nodata value
    :param method:          interpolation method (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, etc.)
    :param fill_value:      value to fill into positions where no interpolation is possible
                            - if 'extrapolate': extrapolate the missing values
    :return:    interpolated array
    """
    assert axis <= 2
    if data.ndim not in [2, 3]:
        raise ValueError('Expected a 2D or 3D array. Received a %dD array.' % data.ndim)
    if isinstance(nodata, np.ndarray) and nodata.shape != data.shape:
        raise ValueError('No-data mask and data must have the same shape.')

    if data.ndim == 2:
        return interp_nodata_along_axis_2d(data, axis=axis, nodata=nodata, method=method, fill_value=fill_value)

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
                    axis=axis if axis != 2 else 1,
                    method=method, fill_value=fill_value))


def interp_nodata_spatially_2d(data_2d: np.ndarray,
                               axis: int = 0,
                               nodata: Union[np.ndarray, Number] = np.nan,
                               method: str = 'linear',
                               fill_value: float = np.nan,
                               implementation: str = 'pandas'
                               ) -> np.ndarray:
    """Interpolate a 2D array spatially.

    NOTE: If data_2d contains NaN values that are unlabelled by a given nodata array,
          they are also overwritten in the pandas implementation.

    :param data_2d:         data to interpolate
    :param axis:            axis to interpolate (0: along columns; 1: along rows)
    :param nodata:          nodata array in the shape of data or nodata value
    :param method:          interpolation method
                            - if implementation='scipy': ‘linear’, ‘nearest’, ‘cubic’
                            - if implementation='pandas': ‘linear’, ‘nearest’, 'slinear’, ‘quadratic’, ‘cubic’, etc.
    :param fill_value:      value to fill into positions where no interpolation is possible
    :param implementation:  'scipy': interpolation based on scipy.interpolate.griddata
                            'pandas': interpolation based on pandas.core.resample.Resampler.interpolate
    :return:    interpolated array
    """
    assert axis < 2
    if data_2d.ndim != 2:
        raise ValueError('Expected a 2D array. Received a %dD array.' % data_2d.ndim)

    badmask_full = _get_baddata_mask(data_2d, nodata)

    if badmask_full.any():
        if implementation == 'scipy':
            if axis == 0:
                y, x = np.indices(data_2d.shape)
            else:
                x, y = np.indices(data_2d.shape)

            data_2d[badmask_full] = \
                griddata(np.array([x[~badmask_full], y[~badmask_full]]).T,  # points we know
                         data_2d[~badmask_full],  # values we know
                         np.array([x[badmask_full], y[badmask_full]]).T,  # points to interpolate
                         method=method, fill_value=fill_value)

        elif implementation == 'pandas':
            data2int = data_2d.astype(float)
            data2int[badmask_full] = np.nan

            data_2d = np.array(DataFrame(data2int)
                               .interpolate(method=method, axis=axis)).astype(data_2d.dtype)

            if fill_value:
                data_2d[np.isnan(data_2d)] = fill_value

        else:
            raise ValueError(implementation, 'Unknown implementation.')

    return data_2d


def interp_nodata_spatially_3d(data_3d: np.ndarray,
                               axis: int = 0,
                               nodata: Union[np.ndarray, Number] = np.nan,
                               method: str = 'linear',
                               fill_value: float = np.nan,
                               implementation: str = 'pandas',
                               CPUs: int = None
                               ) -> np.ndarray:
    """Interpolate a 3D array spatially, band-for-band.

    :param data_3d:         data to interpolate
    :param axis:            axis to interpolate (0: along columns; 1: along rows)
    :param nodata:          nodata array in the shape of data or nodata value
    :param method:          interpolation method
                            - if implementation='scipy': ‘linear’, ‘nearest’, ‘cubic’
                            - if implementation='pandas': ‘linear’, ‘nearest’, 'slinear’, ‘quadratic’, ‘cubic’, etc.
    :param fill_value:      value to fill into positions where no interpolation is possible
    :param implementation:  'scipy': interpolation based on scipy.interpolate.griddata
                            'pandas': interpolation based on pandas.core.resample.Resampler.interpolate
    :param CPUs:            number of CPUs to use
    :return:    interpolated array
    """
    assert axis < 2

    badmask_full = _get_baddata_mask(data_3d, nodata)

    if CPUs > 1:
        with Pool(CPUs or cpu_count()) as pool:
            args = [[data_3d[:, :, band], axis, badmask_full[:, :, band], method, fill_value, implementation]
                    for band in range(data_3d.shape[2])]
            results = pool.starmap(interp_nodata_spatially_2d, args)

            pool.close()  # needed for coverage to work in multiprocessing
            pool.join()

        return np.dstack(results)

    else:
        return \
            np.dstack([interp_nodata_spatially_2d(data_3d[:, :, band], axis=axis,
                                                  nodata=badmask_full[:, :, band], method=method,
                                                  fill_value=fill_value, implementation=implementation)
                       for band in range(data_3d.shape[2])])
