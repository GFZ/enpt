# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2026 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
# (GFZ Potsdam, danschef@gfz.de), Niklas Bohn (GFZ Potsdam, nbohn@gfz.de),
# Stéphane Guillaso (GFZ Potsdam, stephane.guillaso@gfz.de)
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
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""
EnPT de-striping module.

improves the image quality in case of image striping


Input
- L1B data

Output
- L1B data

Memory Budget

Necessary Functions

Process Work Flow
- estimate striping pattern
- correct striping pattern
"""

import numpy as np
from scipy import ndimage
from scipy.signal import savgol_filter, convolve2d
from numpy.fft import fft, ifft

__author__ = ['Maximilian Brell', 'Daniel Scheffler']


def destripe_rad_band_wise(img, detrend_window_length=25, lf=False):
    """correct striping artifacts based on single band L1B EnMAP image.

    :param img:
        L1B (sensor geometry) single band radiance numpy array [ALT, ACT]
        For along-track de-striping [ACT, ALT]

    :param detrend_window_length:
        threshold defining window length for de-trending of cumulated column median.
    """
    img0 = np.copy(img)
    # img = ndimage.uniform_filter1d(img, 25, 0)
    # img01 = ndimage.uniform_filter1d(img, 5, 1)
    img01 = savgol_filter(img, 7, 3, axis=1)

    # calculate gradients in x and y axes
    dx = convolve2d(img, [[1, -1]], boundary='symm', mode='same')
    dy = convolve2d(img, [[1, 1], [-1, -1]], boundary='symm', mode='same')
    dx01 = convolve2d(img01, [[1, -1]], boundary='symm', mode='same')

    # smooth gradients along x- and y-axis (mandatory but eliminates some outliers efficiently)
    smooth_dy = ndimage.uniform_filter1d(dy, 3, 1)
    smooth_dx = ndimage.uniform_filter1d(dx, 3, 0)
    smooth_dx01 = ndimage.uniform_filter1d(dx01, 3, 0)

    # smooth_dx = savgol_filter(dx, 5, 1, axis=1)

    # mask heterogeneous areas (initially a hard-coded threshold (mask_th = 30) made a good job)
    mask_th = np.max(np.percentile(np.abs(dy - smooth_dy), 5, axis=0))
    smooth_dx[np.abs(dy - smooth_dy) > mask_th] = np.nan
    smooth_dx[0, :] = np.nan
    smooth_dx01[0, :] = np.nan

    if lf:
        detrend_window_length = 999
        smooth_dx = ndimage.uniform_filter1d(dx, detrend_window_length, 0)
        smooth_dx0 = np.cumsum(np.nanmedian(smooth_dx, 0))
        smooth_dx1 = np.cumsum(np.nanmedian(smooth_dx01, 0))
        # smooth_dx1 = savgol_filter(smooth_dx0, detrend_window_length, 3)
    else:
        # cumulate column median in x direction
        smooth_dx0 = np.cumsum(np.nanmedian(smooth_dx, 0))
        # smooth_dx1 = np.cumsum(np.nanmedian(smooth_dx01, 0))
        # de-trend based on savgol filter
        smooth_dx1 = savgol_filter(smooth_dx0, detrend_window_length, 3)
    # subtract high-frequent across-track gradient form image
    img0 -= smooth_dx0[None, :] - smooth_dx1
    # img0 -= smooth_dx0 - smooth_dx1
    smooth_dx2 = np.nanmedian(dx, 0)

    return img0, smooth_dx0, smooth_dx2


class Destriper:
    def __init__(self):
        self.high_freq: bool = True
        self.low_freq: bool = False
        self.spatial_domain: bool = True
        self.spectral_domain: bool = False
        self.mode: str = ''
        self.along_track_direction: bool = False

    def destripe(self,
                 array,
                 sensor: str = '',
                 high_freq: bool = True,
                 low_freq: bool = False,
                 spatial_domain: bool = True,
                 spectral_domain: bool = False,
                 mode: str = 'stripes',
                 along_track_direction: bool = False
                 ):
        """Run destriping on L1B EnMAP image.

        :type array:
            L1B (sensor geometry) radiance array

        :param sensor:
            name supplement, e.g., 'vnir', 'swir', or 'vswir'

        :param high_freq:
            correct high frequent features

        :param low_freq:
            correct low frequent features

        :param spatial_domain:
            correct features in the spatial domain

        :param spectral_domain:
            correct features in the spectral domain

        :param mode:
            correct features which are stable over entire along-track dimension = 'stripes'
            correct features which are variable over along-track dimension = 'flicker'

        :param along_track_direction:
            correct pattern in along_track_direction

        :return:    destriped array, difference between original and destriped array
        """
        if along_track_direction:
            array = np.rot90(array, k=1)

        # initialize parameters
        self.high_freq = high_freq
        self.low_freq = low_freq
        self.spatial_domain = spatial_domain
        self.spectral_domain = spectral_domain
        self.mode = mode
        modes = ['stripes', 'flicker']
        if mode not in modes:
            raise ValueError("Invalid mode. Expected one of: %s" % modes)
        self.along_track_direction = along_track_direction

        destriped_data = []

        if self.spatial_domain:
            # iterate through the bands (spatial domain)
            for band in range(0, array.shape[2]):
                img = array[:, :, band].astype(np.float32)
                if np.all(img == 0.0):
                    img.fill(np.nan)
                if self.checking_overall_integrity(img):
                    img, smooth_dx0, smooth_dx1 = self.stripe_detection(img)
                destriped_data.append(img[:, :, None])
            destriped_data = np.concatenate(destriped_data, axis=2)

        if self.spectral_domain:
            # iterating trough columns spectral_domain
            for col in range(0, array.shape[1]):
                img = array[:, col, :]
                if self.checking_overall_integrity(img):
                    img = self.stripe_detection(img)
                destriped_data.append(img[:, :, None])
            destriped_data = np.concatenate(destriped_data, axis=2)

        if along_track_direction:
            if sensor == 'swir':
                destriped_data = self.swir_alt_threshold(corrected_array=destriped_data, original_array=array)

            destriped_data = np.rot90(destriped_data, k=-1)
            array = np.rot90(array, k=-1)

        return destriped_data, destriped_data - array

    def stripe_detection(self, img):
        if self.high_freq:
            img, smooth_dx0, smooth_dx1 = destripe_rad_band_wise(img, detrend_window_length=25)
        if self.low_freq:
            img, smooth_dx0, smooth_dx1 = destripe_rad_band_wise(img, lf=True)
            # img, smooth_dx0, smooth_dx1 = destripe_rad_band_wise(img, detrend_window_length=1000, lf=True)
        return img, smooth_dx0, smooth_dx1

    @staticmethod
    def checking_overall_integrity(img):
        return np.sum(img) != 0.0

    @staticmethod
    def swir_alt_threshold(corrected_array, original_array, threshold=0.60):
        # get rid of values < 0
        wh_zero = np.where(np.logical_or(original_array <= 0.0, corrected_array <= 0.0))
        original_array[wh_zero] = np.nan
        corrected_array[wh_zero] = np.nan

        # calculate along-track percentage miscalibration per band
        mis_cal = np.nanmedian((original_array - corrected_array) / (original_array / 100), axis=0)
        mis_cal = mis_cal - savgol_filter(mis_cal, 15, 3, axis=0)

        # cross-correlation along-track per band
        fft_mis_cal = fft(mis_cal, axis=0)
        cross_mis_cal = ifft(fft_mis_cal * np.conjugate(fft_mis_cal), axis=0).real

        # normalize and cross-correlation max identification
        norm_cross_mis_cal = cross_mis_cal / np.nanmax(cross_mis_cal, axis=0)
        max_cross_mis_cal = np.mean(np.sort(norm_cross_mis_cal[1:, :], axis=0)[-200:, :], axis=0)
        wh_band_mis_cal = max_cross_mis_cal < threshold
        corrected_array[:, :, wh_band_mis_cal] = original_array[:, :, wh_band_mis_cal]
        return corrected_array
