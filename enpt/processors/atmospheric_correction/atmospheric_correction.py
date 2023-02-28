# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2018-2023 Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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
# with this program.  If not, see <https://www.gnu.org/licenses/>.

"""EnPT atmospheric correction module.

Performs the atmospheric correction of EnMAP L1B data.
"""
import pprint
import numpy as np
from multiprocessing import cpu_count
from logging import Logger

from packaging.version import parse as parse_version
from sicor.sicor_enmap import sicor_ac_enmap
from sicor.options import get_options as get_ac_options

from ...model.images import EnMAPL1Product_SensorGeo
from ...options.config import EnPTConfig
from ...utils.path_generator import get_path_ac_options

__author__ = 'Daniel Scheffler'


class AtmosphericCorrector(object):
    """Class for performing atmospheric correction of EnMAP L1 images using SICOR."""

    def __init__(self, config: EnPTConfig = None):
        """Create an instance of AtmosphericCorrector."""
        self.cfg = config

    def _get_sicor_options(self, enmap_ImageL1: EnMAPL1Product_SensorGeo, land_only=False) -> dict:
        """Get a dictionary containing the SICOR options.

        :param enmap_ImageL1:   EnMAP Level 1 product in sensor geometry
        :param land_only:       True: SICOR is applied to land only; False: SICOR is applied to all pixels
        :return:                dictionary of SICOR options
        """
        path_opts = get_path_ac_options()

        try:
            options = get_ac_options(path_opts, validation=True)

            # adjust options
            if enmap_ImageL1.meta.aot is not None:
                options["retrieval"]["default_aot_value"] = enmap_ImageL1.meta.aot

            options["retrieval"]["cpu"] = self.cfg.CPUs or cpu_count()
            options["retrieval"]["disable_progressbars"] = self.cfg.disable_progress_bars

            # TODO: issue is closed -> revise
            # temporarily disable uncertainty measures to avoid https://git.gfz-potsdam.de/EnMAP/sicor/-/issues/86
            # if set to False, uncertainty values are not contained in the additional output of SICOR
            options["retrieval"]["inversion"]["full"] = False

            # set land_only mode
            options["retrieval"]["land_only"] = land_only

            # disable first guess water vapor retrieval for now
            options["retrieval"]["state_vector"]["water_vapor"]["use_prior_mean"] = True
            options["retrieval"]["state_vector"]["water_vapor"]["prior_mean"] = \
                enmap_ImageL1.meta.water_vapour  # = default = 2.5

            # disable first guess liquid water retrieval for now
            options["retrieval"]["state_vector"]["liquid_water"]["use_prior_mean"] = True

            # disable first guess ice retrieval for now
            options["retrieval"]["state_vector"]["ice"]["use_prior_mean"] = True

        except FileNotFoundError:
            raise FileNotFoundError(f'Could not locate options file for atmospheric correction at {path_opts}')

        enmap_ImageL1.logger.debug('SICOR AC configuration: \n' +
                                   pprint.pformat(options))

        return options

    def _is_acwater_operable(self, logger: Logger):
        """Return True if ACWater/Polymer is operable, else raise a warning and return False."""
        try:
            import acwater as _acwater  # noqa: F401

            if parse_version(_acwater.__version__) < parse_version('0.3.0'):
                if self.cfg.mode_ac in ['water', 'combined']:
                    logger.warning(f"The installed version of ACwater (v{_acwater.__version__}) is too old. "
                                   f"At least version 0.3.0 is required. Instead of ACwater, SICOR is applied to water "
                                   f"surfaces as a workaround.")

                return False

        except ImportError as e:
            if self.cfg.mode_ac in ['water', 'combined']:
                logger.warning(f"The atmospheric correction mode was set to '{self.cfg.mode_ac}' but "
                               f"ACwater cannot be imported (Error was: {e.msg}). "
                               f"As a fallback, SICOR is applied to water surfaces instead.")
            return False

        try:
            from acwater.acwater import polymer_ac_enmap as _polymer_ac_enmap
            if not _polymer_ac_enmap:
                logger.warning("Polymer is not callable. "
                               "As a fallback, SICOR is applied to water surfaces instead.")
                return False
        except ImportError as e:
            if self.cfg.mode_ac in ['water', 'combined']:
                logger.warning(f"The atmospheric correction mode was set to '{self.cfg.mode_ac}' but "
                               f"Polymer cannot be imported (Error was: {e.msg}). "
                               f"As a fallback, SICOR is applied to water surfaces instead.")
            return False

        return True

    def _run_AC__land_mode(self,
                           enmap_ImageL1: EnMAPL1Product_SensorGeo
                           ) -> (np.ndarray, np.ndarray, dict):
        """Run atmospheric correction in 'land' mode, i.e., use SICOR for all surfaces."""
        if 2 in enmap_ImageL1.vnir.mask_landwater[:]:
            enmap_ImageL1.logger.info(
                "Running atmospheric correction in 'land' mode, i.e., SICOR is applied to ALL surfaces. "
                "Uncertainty is expected for water surfaces because SICOR is designed for land only.")

        # run SICOR
        # NOTE: - boa_ref_vnir, boa_ref_swir: reflectance between 0 and 1
        #       - res: a dictionary containing retrieval maps with path lengths of the three water phases
        #              and several retrieval uncertainty measures
        #              -> cwv_model, liq_model, ice_model, intercept_model, slope_model, toa_model,
        #                 sx (posterior predictive uncertainty matrix), scem (correlation error matrix),
        #                 srem (relative error matrix)
        #                 optional (if options["retrieval"]["inversion"]["full"] = True):
        #                 jacobian, convergence, iterations, gain, averaging_kernel, cost_function,
        #                 dof (degrees of freedom), information_content, retrieval_noise, smoothing_error
        #              -> SWIR geometry (?)  # FIXME
        boa_ref_vnir, boa_ref_swir, land_additional_results = \
            sicor_ac_enmap(enmap_l1b=enmap_ImageL1,
                           options=self._get_sicor_options(enmap_ImageL1, land_only=False),
                           unknowns=True,
                           logger=enmap_ImageL1.logger)

        return boa_ref_vnir, boa_ref_swir, land_additional_results

    def _run_AC__water_mode(self, enmap_ImageL1: EnMAPL1Product_SensorGeo
                            ) -> (np.ndarray, np.ndarray):
        """Run atmospheric correction in 'water' mode, i.e., use ACWater/Polymer for water surfaces only.

        NOTE:
            - Land surfaces are NOT included in the EnMAP L2A product.
            - The output radiometric unit for water surfaces is 'water leaving reflectance'.
        """
        from acwater.acwater import polymer_ac_enmap

        if 1 in enmap_ImageL1.vnir.mask_landwater[:]:
            enmap_ImageL1.logger.info(
                "Running atmospheric correction in 'water' mode, i.e., ACWater/Polymer is applied to water "
                "surfaces only. Note that land surfaces will NOT be included in the EnMAP L2A product.")

        # run ACWater/Polymer for water surfaces only
        # NOTE: polymer_ac_enmap() returns masked (nan) values for land
        #       - res: a dictionary containing retrieval maps with several additional retrieval measures
        #              -> chla, bitmask, logfb, Rnir, Rgli
        try:
            wl_ref_vnir, wl_ref_swir, water_additional_results = \
                polymer_ac_enmap(enmap_l1b=enmap_ImageL1,
                                 config=self.cfg,
                                 detector='vnir')
        except:  # noqa
            enmap_ImageL1.logger.error(
                "The atmospheric correction for water surfaces based on ACwater/Polymer failed (issue tracker at "
                "https://gitlab.awi.de/phytooptics/acwater/-/issues).\n"
                "Alternatively, you may run EnPT in the 'land' atmospheric correction mode based on SICOR.\n"
                "The error message is now raised:"
            )
            raise

        return wl_ref_vnir, wl_ref_swir, water_additional_results

    def _run_AC__combined_mode(self,
                               enmap_ImageL1: EnMAPL1Product_SensorGeo
                               ) -> (np.ndarray, np.ndarray, dict):
        """Run atmospheric corr. in 'combined' mode, i.e., use SICOR for land and ACWater/Polymer for water surfaces.

        NOTE:
            - The output radiometric units are:
                - 'surface reflectance' for land surfaces
                - 'water leaving reflectance' for water surfaces
            - There might be visible edge effects, e.g., at coastlines.
        """
        from acwater.acwater import polymer_ac_enmap

        only = 'water' if 1 not in enmap_ImageL1.vnir.mask_landwater[:] else 'land'
        if 1 not in enmap_ImageL1.vnir.mask_landwater[:] or \
           2 not in enmap_ImageL1.vnir.mask_landwater[:]:
            enmap_ImageL1.logger.info(
                f"Running atmospheric correction in 'combined' mode, i.e., SICOR is applied to land and "
                f"ACWater/Polymer is applied to water surfaces. But the input image only contains {only} surfaces.")

        # run SICOR for land surfaces only
        boa_ref_vnir_land, boa_ref_swir_land, land_additional_results = \
            sicor_ac_enmap(enmap_l1b=enmap_ImageL1,
                           options=self._get_sicor_options(enmap_ImageL1, land_only=True),
                           unknowns=True,
                           logger=enmap_ImageL1.logger)

        # run ACWater/Polymer for water surfaces only
        # NOTE: polymer_ac_enmap() returns masked (nan) values for land
        #       - res: a dictionary containing retrieval maps with several additional retrieval measures
        #              -> chla, bitmask, logfb, Rnir, Rgli
        try:
            wl_ref_vnir_water, wl_ref_swir_water, water_additional_results = \
                polymer_ac_enmap(enmap_l1b=enmap_ImageL1,
                                 config=self.cfg,
                                 detector='vnir')
        except:  # noqa
            enmap_ImageL1.logger.error(
                "The atmospheric correction for water surfaces based on ACwater/Polymer failed (issue tracker at "
                "https://gitlab.awi.de/phytooptics/acwater/-/issues).\n"
                "Alternatively, you may run EnPT in the 'land' atmospheric correction mode based on SICOR.\n"
                "The error message is now raised:"
            )
            raise

        # Overwrite the SICOR output at water positions with the output from ACwater/Polymer
        # -> output contains Water-leaving-reflectance over water and BOA-reflectance over land
        water_mask_vnir_3D = (enmap_ImageL1.vnir.mask_landwater[:] == 2)[:, :, None]  # 2 = water
        wlboa_ref_vnir = np.where(water_mask_vnir_3D, wl_ref_vnir_water, boa_ref_vnir_land)
        wlboa_ref_swir = np.where(water_mask_vnir_3D, wl_ref_swir_water, boa_ref_swir_land)

        return wlboa_ref_vnir, wlboa_ref_swir, water_additional_results, land_additional_results

    @staticmethod
    def _validate_AC_results(reflectance_vnir: np.ndarray,
                             reflectance_swir: np.ndarray,
                             logger: Logger):
        for detectordata, detectorname in zip([reflectance_vnir, reflectance_swir],
                                              ['VNIR', 'SWIR']):
            mean0 = np.nanmean(detectordata[:, :, 0])
            std0 = np.nanstd(detectordata[:, :, 0])

            if np.isnan(mean0) or \
               mean0 == 0 or \
               std0 == 0:
                logger.warning(f'The atmospheric correction returned empty {detectorname} bands!')

    def run_ac(self,
               enmap_ImageL1: EnMAPL1Product_SensorGeo
               ) -> EnMAPL1Product_SensorGeo:
        """Run atmospheric correction according to the specified 'mode_ac' parameter.

        :param enmap_ImageL1:   input EnMAP image containing TOA reflectance (an instance EnMAPL1Product_SensorGeo)
        :return:    atmospherically corrected output EnMAP image containing BOA reflectance / water leaving reflectance
                    (an instance EnMAPL1Product_SensorGeo)
        """
        enmap_ImageL1.set_SWIRattr_with_transformedVNIRattr('mask_landwater')

        enmap_ImageL1.logger.info(
            f"Starting atmospheric correction for VNIR and SWIR detector in '{self.cfg.mode_ac}' mode. "
            f"Source radiometric unit code is '{enmap_ImageL1.meta.vnir.unitcode}'.")

        # set initial values for land_additional_results and water_additional_results
        land_additional_results = None
        water_additional_results = None

        # run the AC
        if self.cfg.mode_ac in ['water', 'combined'] and not self._is_acwater_operable(enmap_ImageL1.logger):
            # use SICOR as fallback AC for water surfaces if ACWater/Polymer is not installed
            reflectance_vnir, reflectance_swir, land_additional_results = \
                self._run_AC__land_mode(enmap_ImageL1)

        else:
            if self.cfg.mode_ac == 'combined':
                reflectance_vnir, reflectance_swir, water_additional_results, land_additional_results = \
                    self._run_AC__combined_mode(enmap_ImageL1)

            elif self.cfg.mode_ac == 'water':
                reflectance_vnir, reflectance_swir, water_additional_results = \
                    self._run_AC__water_mode(enmap_ImageL1)

            elif self.cfg.mode_ac == 'land':
                reflectance_vnir, reflectance_swir, land_additional_results = \
                    self._run_AC__land_mode(enmap_ImageL1)

            else:
                raise ValueError(self.cfg.mode_ac,
                                 "Unexpected 'mode_ac' parameter. "
                                 "Choose one out of 'land', 'water', 'combined'.")

        # validate outputs
        self._validate_AC_results(reflectance_vnir, reflectance_swir, logger=enmap_ImageL1.logger)

        # join results
        enmap_ImageL1.logger.info('Joining results of atmospheric correction.')

        for in_detector, out_detector in zip([enmap_ImageL1.vnir, enmap_ImageL1.swir],
                                             [reflectance_vnir, reflectance_swir]):
            in_detector.data = (out_detector * self.cfg.scale_factor_boa_ref).astype(np.int16)
            # NOTE: geotransform and projection are missing due to sensor geometry

            in_detector.detector_meta.unit = '0-%d' % self.cfg.scale_factor_boa_ref
            in_detector.detector_meta.unitcode = 'BOARef'

        # FIXME: Consider to also join SICOR's land_additional_results
        #  (contains three phases of water maps and several retrieval uncertainty measures)

        # join additional results from ACwater/Polymer
        if water_additional_results and self.cfg.polymer_additional_results:

            water_mask = enmap_ImageL1.vnir.mask_landwater[:] == 2
            for k in water_additional_results.keys():
                if k == 'polymer_bitmask':
                    # the bitmask already explicitly indicates land pixels with "1"
                    continue
                else:
                    v = water_additional_results[k]
                    v[~water_mask] = -9999
                    v[np.isnan(v)] = -9999

                    water_additional_results[k] = v

            enmap_ImageL1.vnir.polymer_logchl = water_additional_results['polymer_logchl']
            enmap_ImageL1.vnir.polymer_logfb = water_additional_results['polymer_logfb']
            enmap_ImageL1.vnir.polymer_rgli = water_additional_results['polymer_rgli']
            enmap_ImageL1.vnir.polymer_rnir = water_additional_results['polymer_rnir']
            enmap_ImageL1.vnir.polymer_bitmask = water_additional_results['polymer_bitmask']

        return enmap_ImageL1
