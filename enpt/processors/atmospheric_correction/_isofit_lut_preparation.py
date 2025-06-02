# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2024 Karl Segl (GFZ Potsdam, segl@gfz.de), Daniel Scheffler
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

"""Module to prepare the pre-built LUT passed to ISOFIT."""

import logging
from typing import Union

import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d

from ...utils.logging import EnPT_Logger


class LUTTransformer(object):
    def __init__(self, path_lut: str, sza_scene: float, logger: Union[EnPT_Logger, logging.Logger] = None):
        """Get an instance of LUTTransformer.

        :param path_lut:    atmospheric look-up-table in binary format as created by Luis
        """
        self.p_lut_bin = path_lut
        self._offset = 0
        self.sza_scene = sza_scene
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger()
            self.logger.setLevel(logging.INFO)

    def read_binary_modtran_lut(self, path_out_nc: str):
        """
        Read MODTRAN® LUT.

        :param path_out_nc: path to output LUT
        :return:         LUT of atmospheric functions, x and y axes grid points, LUT wavelengths
        """
        if self.logger:
            self.logger.debug("Reading binary LUT")

        def read_int16(array: np.ndarray, count):
            val = np.array(array[self._offset:self._offset + count * 2].view('int16'))
            self._offset += count * 2
            return val

        def read_float32(array: np.ndarray, count):
            val = np.array(array[self._offset:self._offset + count * 4].view('f4'))
            self._offset += count * 4
            return val

        with open(self.p_lut_bin, 'rb') as fd:
            data = np.frombuffer(fd.read(), dtype=np.uint8)  # Read all data as bytes

        wvl, vza, sza, alt, aot, raa, cwv = [read_float32(data, count=read_int16(data, count=1)[0]) for _ in range(7)]
        npar1, npar2 = [read_int16(data, count=1)[0] for _ in range(2)]

        # LUT1 containing path radiance
        # NOTE: The rhoatm values do not vary much with changing CVW contents,
        #       therefore the LUT was created with a single CWV value.
        # lut1_axnames = ['VZA', 'SZA', 'ALT', 'AOT', 'RAA', 'CWV', 'WVL', 'PAR1']
        lut1 = read_float32(data, count=10584000).reshape(
            (len(vza), len(sza), len(alt), len(aot), len(raa), 1, len(wvl), npar1))
        # LUT2 containing downwelling direct and diffuse surface solar irradiance, spherical albedo and transmittance
        # NOTE: The values in LUT2 do not vary much with changing RAAs,
        #       therefore the LUT was created with a single RAA value.
        # lut2_axnames = ['VZA', 'SZA', 'ALT', 'AOT', 'RAA', 'CWV', 'WVL', 'PAR2']
        lut2 = read_float32(data, count=42336000).reshape(
            (len(vza), len(sza), len(alt), len(aot), 1, len(cwv), len(wvl), npar2))

        ##############################
        # TRANSFORM TO ISOFIT FORMAT #
        ##############################

        if self.logger:
            self.logger.debug("Transforming LUT to ISOFIT format")

        # Get data from LUT1 and process it
        fwhm = np.zeros_like(wvl)
        fwhm[:len(wvl) - 1] = np.diff(wvl)
        fwhm[len(wvl) - 1] = fwhm[len(wvl) - 2]
        cos_sza = np.cos(np.deg2rad(sza))[None, :, None, None, None, None, None]

        # NOTE: ISOFIT expects a LUT with all parameters simulated for multiple CWV values.
        #       Since the binary MODTRAN LUT1 is created with a single CWV value (see above),
        #       LUT1 has to be replicated 7 times for the CWV axis (axis=5). Same applies
        #       to LUT2 and RAA (axis=4).
        _lut1 = np.repeat(lut1, 7, axis=5)
        _lut2 = np.repeat(lut2, 7, axis=4)

        # TODO: reduce the LUT according to the state vector of the current EnMAP image

        # extract all the parameters from the 2 LUTs
        rhoatm = _lut1[..., 0]  # rhoatm
        edir = _lut2[..., 0]  # transm_down_dir
        edif = _lut2[..., 1]  # transm_down_dif
        sab = _lut2[..., 2]  # spherical albedo
        # The binary MODTRAN LUT has values in _lut2[..., 3], which should be transm_up_dir + transm_up_dif.
        # However, this is already contained in the Edir/Edif variables, thus we need to set it to 1 here.
        var5 = 1  # transm_up_dir + transm_up_dif
        var6 = 0  # not needed ... already 0 and the transm_up_dir includes the transm_up_dif already

        # convert parameters to what ISOFIT expects
        _unit_scalefactor = 1000  # adapt MODTRAN LUT unit to ISOFIT
        rhoatm *= _unit_scalefactor  # path radiance
        sphalb = sab  # atmospheric spherical albedo
        # tg = tdir + tdif
        # transm_down_dir  = tg * cos(sza) * edir /  PI
        # transm_down_dif  = tg * edif /  PI
        _constant = (var5 + var6) / np.pi * _unit_scalefactor
        transm_down_dir = edir * cos_sza * _constant  # direct downward transmittance
        transm_down_dif = edif * _constant  # (diffuse downward transmittance)

        del lut1, lut2, _lut1, _lut2, edir, edif, sab, data

        # extrapolate data at 8 km elevation due to a bug in the MODTRAN-LUT (data at 8km = data at 0 km)
        for arr in [rhoatm, sphalb, transm_down_dir, transm_down_dif]:
            self.extrapolate_8km(arr, alt)

        transm_up_dir = np.zeros_like(rhoatm)  # direct upward transmittance
        transm_up_dif = np.zeros_like(rhoatm)  # diffuse upward transmittance
        thermal_upwelling = np.zeros_like(rhoatm)  # thermal up-welling
        thermal_downwelling = np.zeros_like(rhoatm)  # thermal down-welling

        if self.logger:
            self.logger.debug(f"Writing LUT to {path_out_nc}")

        # write LUT to NetCDF file
        with nc.Dataset(path_out_nc, 'w', format='NETCDF4') as nc_out:
            nc_out.setncattr("RT_mode", "rdn")

            var = nc_out.createVariable('coszen', 'f4', ())
            var.assignValue(np.cos(np.deg2rad(self.sza_scene)))

            # Add dimensions
            for t, v in zip(('H2OSTR', 'AERFRAC_2', 'surface_elevation_km', 'solar_zenith',
                             'observer_zenith', 'relative_azimuth', 'wl'),
                            (cwv, aot, alt, sza, vza, raa, wvl)):
                nc_out.createDimension(t, v.size)

            # Add data
            dims = ('observer_zenith', 'solar_zenith', 'surface_elevation_km',
                    'AERFRAC_2', 'relative_azimuth', 'H2OSTR', 'wl')
            for t, d, v in [
                # 1D data
                ('H2OSTR', ('H2OSTR',), cwv),
                ('AERFRAC_2', ('AERFRAC_2',), aot),
                ('surface_elevation_km', ('surface_elevation_km',), alt),
                ('solar_zenith', ('solar_zenith',), sza),
                ('observer_zenith', ('observer_zenith',), vza),
                ('relative_azimuth', ('relative_azimuth',), raa),
                ('wl', ('wl',), wvl),
                ('fwhm', ('wl',), fwhm),
                ('solar_irr', ('wl',), np.zeros_like(wvl)),  # not used by ISOFIT)
                # 7D data
                ('rhoatm', dims, rhoatm),
                ('sphalb', dims, sphalb),
                ('transm_down_dir', dims, transm_down_dir),
                ('transm_down_dif', dims, transm_down_dif),
                ('transm_up_dir', dims, transm_up_dir),
                ('transm_up_dif', dims, transm_up_dif),
                ('thermal_upwelling', dims, thermal_upwelling),
                ('thermal_downwelling', dims, thermal_downwelling),
            ]:
                var = nc_out.createVariable(t, "f4", d)
                var[:] = v

        if self.logger:
            self.logger.debug("LUT successfully saved")

    @staticmethod
    def extrapolate_8km(var: np.ndarray, alt_grid: np.ndarray):
        """Replace data at 8km altitude by linearly extrapolating with 0.7km and 2.5km as grid points.

        :param var:         variable to extrapolate
        :param alt_grid:    altitude grid points where to extrapolate
        :return:
        """
        var[:, :, 3, :, :, :, :] = \
            interp1d(x=alt_grid[1:3],
                     y=var[:, :, 1:3, :, :, :, :],
                     axis=2,
                     fill_value='extrapolate',
                     bounds_error=False)(
                np.array([8.])
            )[:, :, 0, :, :, :, :]

        return var
