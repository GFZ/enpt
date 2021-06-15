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

"""EnPT module for handling spectral response functions."""

from typing import Union, List

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

__author__ = 'Daniel Scheffler'


class SRF(object):
    def __init__(self, wvl_unit: str = 'nanometers', wvl_min: float = 400, wvl_max: float = 2500, specres_nm: float = 1,
                 format_bandnames: bool = False, v: bool = False):
        """SRF instance provides SRF functions, wavelength positions, etc..

        :param wvl_unit:            the wavelengths unit to be used within SRF instance ('nanometers' or 'micrometers)
        :param wvl_min:
        :param wvl_max:
        :param specres_nm:          output spectral resolution of SRFs in nanometers
        :param format_bandnames:    whether to format default strings from LayerBandsAssignment as 'B01', 'B02' etc..
        :param v:                   verbose mode
        """
        if wvl_unit not in ['micrometers', 'nanometers']:
            raise ValueError('Unknown wavelength unit %s.' % wvl_unit)

        self.srfs_wvl = []  # wavelength positions with 1 nm precision
        self.srfs = {}  # srf values with 1 nm precision
        self.srfs_norm01 = {}  # srf values with 1 nm precision
        self.bands = []
        self.wvl = None
        self.wvl_unit = wvl_unit
        self.wvl_min = wvl_min
        self.wvl_max = wvl_max
        self.specres_nm = specres_nm
        self.format_bandnames = format_bandnames
        self.conv = {}
        self.v = v

    @staticmethod
    def compute_gaussian_srf(cwl: float, fwhm: float, wvl_min: float, wvl_max: float, wvl_res: float,
                             normalize: bool = True) -> np.ndarray:
        """Compute a spectral response function based on center wavelength and band width using on a gaussian curve.

        :param cwl:         target center wavelength position
        :param fwhm:        target band width (full width half maximum)
        :param wvl_min:     minimum wavelength to compute spectral response for
        :param wvl_max:     maximum wavelength to compute spectral response for
        :param wvl_res:     spectral resolution at which spectral response is to be computed
        :param normalize:   whether to normalize the output spectral response to values between 0 and 1
        :return:            2D numpy.ndarray: rows: response per wavelength; columns: wavelength/response
        """
        x = np.arange(wvl_min, wvl_max, wvl_res)
        dist = stats.norm(cwl, fwhm)
        y = dist.pdf(x)

        if normalize:
            y *= (1.0 / y.max())

        rsp = np.empty((x.size, 2), dtype=float)
        rsp[:, 0] = x
        rsp[:, 1] = y

        return rsp

    @classmethod
    def from_cwl_fwhm(cls, cwls: Union[list, np.ndarray], fwhms: Union[list, np.ndarray], **kwargs: dict) -> 'SRF':
        """Create an instance of SRF based on center wavelength positions and band widths (using gaussian responses).

        :param cwls:    center wavelength positions
        :param fwhms:   band widths
        :param kwargs:  Keyword arguments to be passed to SRF.__init__().
        :return:        SRF instance
        """
        srf = cls(**kwargs)

        if srf.wvl_unit != 'nanometers':
            cwls, fwhms = cwls * 1000, fwhms * 1000

        srf.bands = [str(i + 1) for i in range(len(cwls))]

        # compute SRFs and set attributes
        for bN, cwl, fwhm in zip(srf.bands, cwls, fwhms):
            gaussian_srf = cls.compute_gaussian_srf(cwl, fwhm, srf.wvl_min, srf.wvl_max, srf.specres_nm)
            srf.srfs_wvl = gaussian_srf[:, 0].flatten()
            srf_norm01 = gaussian_srf[:, 1].flatten()
            srf.srfs_norm01[bN] = srf_norm01
            srf.srfs[bN] = srf_norm01 / np.trapz(x=srf.srfs_wvl, y=srf_norm01)

        srf.wvl = np.array(cwls)

        srf.conv.update({key: value for key, value in zip(srf.bands, srf.wvl)})
        srf.conv.update({value: key for key, value in zip(srf.bands, srf.wvl)})

        return srf

    def instrument(self, bands):
        return {
            'rspf': np.vstack([self[band] for band in bands]),
            'wvl_rsp': np.copy(self.srfs_wvl),
            'wvl_inst': np.copy(self.wvl),
            'sol_irr': None
        }

    def convert_wvl_unit(self):
        """Convert the wavelength unit to nanometers if they are in micrometers or vice versa."""
        factor = 1/1000 if self.wvl_unit == 'nanometers' else 1000
        self.srfs_wvl = self.srfs_wvl * factor
        self.wvl = self.wvl * factor
        self.wvl_unit = 'nanometers' if self.wvl_unit == 'micrometers' else 'micrometers'

    def __call__(self, band):
        return self.srfs[band]

    def __getitem__(self, band):
        return self.srfs[band]

    def __iter__(self):
        for band in self.bands:
            yield self[band]

    def plot_srfs(self, figsize: tuple = (15, 5), band: Union[str, List[str]] = None, normalize: bool = True):
        """Show a plot of all spectral response functions.

        :param figsize: figure size of the plot
        :param band:    band key to plot a single band instead of all bands
        :param normalize:    normalize SRFs to 0-1 (default: True)
        """
        if band and band not in self.bands:
            raise ValueError("Parameter 'band' must be a string out of those: %s."
                             % ', '.join(self.bands))

        plt.figure(figsize=figsize)
        bands2plot = [band] if band else self.bands
        for band in bands2plot:
            srfs = list(self.srfs_norm01[band]) if normalize else list(self.srfs[band])
            plt.plot(self.srfs_wvl, srfs, '-', label='Band %s' % band)
        plt.title('EnMAP spectral response functions')
        plt.xlabel('wavelength [%s]' % self.wvl_unit)
        plt.ylabel('spectral response')
        plt.legend(loc='upper right')
