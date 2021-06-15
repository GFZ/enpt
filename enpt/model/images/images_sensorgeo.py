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

"""EnPT EnMAP objects in sensor geometry."""

import logging
from types import SimpleNamespace
from typing import Tuple, Optional, List  # noqa: F401
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import numpy as np
from os import path, makedirs
from glob import glob
import utm
from scipy.interpolate import interp2d
from geoarray import GeoArray

from ...utils.logging import EnPT_Logger
from .image_baseclasses import _EnMAP_Image
from ...model.metadata import EnMAP_Metadata_L1B_SensorGeo, EnMAP_Metadata_L1B_Detector_SensorGeo
from ...model.metadata import EnMAP_Metadata_L2A_MapGeo  # noqa: F401  # only used for type hint
from ...options.config import EnPTConfig
from ...processors.dead_pixel_correction import Dead_Pixel_Corrector
from ...processors.dem_preprocessing import DEM_Processor
from ...processors.spatial_transform import compute_mapCoords_within_sensorGeoDims

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class EnMAP_Detector_SensorGeo(_EnMAP_Image):
    """Class representing a single detector of an EnMAP image (in sensor geometry).

    NOTE:
        - Inherits all attributes from _EnMAP_Image class.
        - All functionality that VNIR and SWIR detectors (sensor geometry) have in common is to be implemented here.

    Attributes:
        - to be listed here. Check help(_EnMAP_Detector_SensorGeo) in the meanwhile!

    """

    def __init__(self, detector_name: str, root_dir: str, config: EnPTConfig, logger=None, meta=None) -> None:
        """Get an instance of _EnMAP_Detector_SensorGeo.

        :param detector_name:   'VNIR' or 'SWIR'
        :param root_dir:
        :param logger:
        :param meta: import meta if already loaded
        """
        if detector_name not in ['VNIR', 'SWIR']:
            raise ValueError("'detector_name' must be 'VNIR' or 'SWIR'. Received %s!" % detector_name)

        self.cfg = config
        self._root_dir = root_dir
        self.detector_name = detector_name
        self.logger = logger or logging.getLogger()

        # get all attributes of base class "_EnMAP_Image"
        super(EnMAP_Detector_SensorGeo, self).__init__()

        # instance an empty metadata object
        if meta is None:
            self.detector_meta: EnMAP_Metadata_L1B_Detector_SensorGeo = \
                EnMAP_Metadata_L1B_Detector_SensorGeo(self.detector_name, config=self.cfg, logger=self.logger)
        else:
            self.detector_meta: EnMAP_Metadata_L1B_Detector_SensorGeo = meta

    def get_paths(self) -> SimpleNamespace:
        """Get all file paths associated with the current instance of EnMAP_Detector_SensorGeo.

        NOTE: This information is read from the detector_meta.

        :return: paths
        """
        self.paths.root_dir = self._root_dir
        self.paths.data = path.join(self._root_dir, self.detector_meta.filename_data)

        def path_or_None(filename):
            return path.join(self._root_dir, filename) if filename else None

        self.paths.mask_landwater = path_or_None(self.detector_meta.filename_mask_landwater)
        self.paths.mask_clouds = path_or_None(self.detector_meta.filename_mask_clouds)
        self.paths.mask_cloudshadow = path_or_None(self.detector_meta.filename_mask_cloudshadow)
        self.paths.mask_haze = path_or_None(self.detector_meta.filename_mask_haze)
        self.paths.mask_snow = path_or_None(self.detector_meta.filename_mask_snow)
        self.paths.mask_cirrus = path_or_None(self.detector_meta.filename_mask_cirrus)
        self.paths.deadpixelmap = path_or_None(self.detector_meta.filename_deadpixelmap)
        self.paths.quicklook = path_or_None(self.detector_meta.filename_quicklook)

        return self.paths

    def correct_dead_pixels(self):
        """Correct dead pixels with respect to the dead pixel mask."""
        algo = self.cfg.deadpix_P_algorithm
        method_spectral, method_spatial = self.cfg.deadpix_P_interp_spectral, self.cfg.deadpix_P_interp_spatial
        self.logger.info("Correcting dead pixels of %s detector...\n"
                         "Used algorithm: %s interpolation in the %s domain"
                         % (self.detector_name, method_spectral, algo if algo == 'spectral' else method_spatial))

        self.data = \
            Dead_Pixel_Corrector(algorithm=algo,
                                 interp_spectral=method_spectral,
                                 interp_spatial=method_spatial,
                                 logger=self.logger)\
            .correct(self.data, self.deadpixelmap)

    def get_preprocessed_dem(self) -> GeoArray:
        """Get a digital elevation model in EnMAP sensor geometry of the current detector."""
        if self.cfg.path_dem:
            self.logger.info('Pre-processing DEM for %s...' % self.detector_name)
            DP = DEM_Processor(self.cfg.path_dem, enmapIm_cornerCoords=tuple(zip(self.detector_meta.lon_UL_UR_LL_LR,
                                                                                 self.detector_meta.lat_UL_UR_LL_LR)),
                               CPUs=self.cfg.CPUs)
            DP.fill_gaps()  # FIXME this will also be needed at other places

            R, C = self.data.shape[:2]
            if DP.dem.is_map_geo:
                lons = self.detector_meta.lons
                lats = self.detector_meta.lats

                if not (lons.ndim == 2 and lats.ndim == 2) and not (lons.ndim == 3 and lats.ndim == 3):
                    raise ValueError((lons.ndim, lats.ndim), 'Geolayer must be either 2- or 3-dimensional.')

                msg_bandinfo = ''
                if lons.ndim == 3:
                    # 3D geolayer (the usual case for EnMAP data provided by DLR)
                    lons = lons[:, :, 0]
                    lats = lats[:, :, 0]
                    msg_bandinfo = ' (using first band of %s geolayer)' % self.detector_name
                else:
                    # 2D geolayer (GFZ-internal test case)
                    # FIXME replace linear interpolation by native geolayers
                    if lons.shape != self.data.shape:
                        lons = self.detector_meta.interpolate_corners(*self.detector_meta.lon_UL_UR_LL_LR, nx=C, ny=R)
                    if lats.shape != self.data.shape:
                        lats = self.detector_meta.interpolate_corners(*self.detector_meta.lat_UL_UR_LL_LR, nx=C, ny=R)

                self.logger.info(('Transforming DEM to %s sensor geometry%s...' % (self.detector_name, msg_bandinfo)))
                self.dem = DP.to_sensor_geometry(lons=lons, lats=lats)
            else:
                self.dem = DP.dem

        else:
            self.logger.info('No DEM for the %s detector provided. Falling back to an average elevation of %d meters.'
                             % (self.detector_name, self.cfg.average_elevation))
            self.dem = GeoArray(np.full(self.data.shape[:2], self.cfg.average_elevation).astype(np.int16))

        return self.dem

    def append_new_image(self, img2: 'EnMAP_Detector_SensorGeo', n_lines: int = None) -> None:
        # TODO convert method to function?
        """Check if a second image matches with the first image and if so, append the given number of lines below.

        In this version we assume that the image will be added below. If it is the case, the method will create
        temporary files that will be used in the following.

        :param img2:
        :param n_lines: number of lines to be added from the new image
        :return: None
        """
        if not self.cfg.is_dummy_dataformat:
            basename_img1 = self.detector_meta.filename_data.split('-SPECTRAL_IMAGE')[0] + '::%s' % self.detector_name
            basename_img2 = img2.detector_meta.filename_data.split('-SPECTRAL_IMAGE')[0] + '::%s' % img2.detector_name
        else:
            basename_img1 = path.basename(self._root_dir)
            basename_img2 = path.basename(img2._root_dir)

        self.logger.info("Check new image for %s: %s " % (self.detector_name, basename_img2))

        distance_min = 27.0
        distance_max = 34.0

        # compute distance between image1 LL and image2 UL
        x1, y1, _, _ = utm.from_latlon(self.detector_meta.lat_UL_UR_LL_LR[2], self.detector_meta.lon_UL_UR_LL_LR[2])
        x2, y2, _, _ = utm.from_latlon(img2.detector_meta.lat_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[0])
        distance_left = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        # compute distance between image1 LR and image2 UR
        x1, y1, _, _ = utm.from_latlon(self.detector_meta.lat_UL_UR_LL_LR[3], self.detector_meta.lon_UL_UR_LL_LR[3])
        x2, y2, _, _ = utm.from_latlon(img2.detector_meta.lat_UL_UR_LL_LR[1], img2.detector_meta.lon_UL_UR_LL_LR[1])
        distance_right = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

        if distance_min < distance_left < distance_max and distance_min < distance_right < distance_max:
            self.logger.info("Append new image to %s: %s" % (self.detector_name, basename_img2))
        else:
            self.logger.warning("%s and %s don't fit to be appended." % (basename_img1, basename_img2))
            return

        # set new number of line
        n_lines = n_lines or img2.detector_meta.nrows

        if n_lines > img2.detector_meta.nrows:
            self.logger.warning("n_lines (%s) exceeds the total number of line of second image" % n_lines)
            self.logger.warning("Set to the image number of line")
            n_lines = img2.detector_meta.nrows

        if n_lines < 50:  # TODO: determine these values
            self.logger.warning("A minimum of 50 lines is required, only %s were selected" % n_lines)
            self.logger.warning("Set the number of line to 50")
            n_lines = 50

        self.detector_meta.nrows += n_lines

        # Compute new lower coordinates
        if not self.cfg.is_dummy_dataformat:
            img2_cornerCoords = tuple(zip(img2.detector_meta.lon_UL_UR_LL_LR,
                                          img2.detector_meta.lat_UL_UR_LL_LR))
            elevation = DEM_Processor(img2.cfg.path_dem,
                                      enmapIm_cornerCoords=img2_cornerCoords).dem \
                if img2.cfg.path_dem else self.cfg.average_elevation

            LL, LR = compute_mapCoords_within_sensorGeoDims(
                sensorgeoCoords_YX=[(n_lines - 1, 0),  # LL
                                    (n_lines - 1, img2.detector_meta.ncols - 1)],  # LR
                rpc_coeffs=list(img2.detector_meta.rpc_coeffs.values())[0],  # RPC coeffs of first band of the detector
                elevation=elevation,
                enmapIm_cornerCoords=img2_cornerCoords,
                enmapIm_dims_sensorgeo=(img2.detector_meta.nrows, img2.detector_meta.ncols)
            )

            self.detector_meta.lon_UL_UR_LL_LR[2], self.detector_meta.lat_UL_UR_LL_LR[2] = LL
            self.detector_meta.lon_UL_UR_LL_LR[3], self.detector_meta.lat_UL_UR_LL_LR[3] = LR
        else:
            # lats
            interp_lats = interp2d(x=[0, 1],
                                   y=[0, 1],
                                   z=[[img2.detector_meta.lat_UL_UR_LL_LR[0], img2.detector_meta.lat_UL_UR_LL_LR[1]],
                                      [img2.detector_meta.lat_UL_UR_LL_LR[2], img2.detector_meta.lat_UL_UR_LL_LR[3]]],
                                   kind='linear')
            self.detector_meta.lat_UL_UR_LL_LR[2] = np.array(interp_lats(0, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lat_UL_UR_LL_LR[3] = np.array(interp_lats(1, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lats = self.detector_meta.interpolate_corners(*self.detector_meta.lat_UL_UR_LL_LR,
                                                                             self.detector_meta.ncols,
                                                                             self.detector_meta.nrows)
            # lons
            interp_lons = interp2d(x=[0, 1],
                                   y=[0, 1],
                                   z=[[img2.detector_meta.lon_UL_UR_LL_LR[0], img2.detector_meta.lon_UL_UR_LL_LR[1]],
                                      [img2.detector_meta.lon_UL_UR_LL_LR[2], img2.detector_meta.lon_UL_UR_LL_LR[3]]],
                                   kind='linear')
            self.detector_meta.lon_UL_UR_LL_LR[2] = np.array(interp_lons(0, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lon_UL_UR_LL_LR[3] = np.array(interp_lons(1, int(n_lines / img2.detector_meta.nrows)))[0]
            self.detector_meta.lons = self.detector_meta.interpolate_corners(*self.detector_meta.lon_UL_UR_LL_LR,
                                                                             self.detector_meta.ncols,
                                                                             self.detector_meta.nrows)

        # append the raster data
        self.data = np.append(self.data, img2.data[0:n_lines, :, :], axis=0)

        # only append masks for the VNIR as they are only provided in VNIR sensor geometry
        if self.detector_name == 'VNIR':
            for attrName in ['mask_landwater', 'mask_clouds', 'mask_cloudshadow',
                             'mask_haze', 'mask_snow', 'mask_cirrus']:
                arr_img1 = getattr(self, attrName)
                arr_img2 = getattr(img2, attrName)

                if arr_img1 is not None and arr_img2 is not None:
                    setattr(self, attrName, np.append(arr_img1, arr_img2[0:n_lines, :], axis=0))
                else:
                    # this mainly applies to the dummy data format that does not have all mask files
                    self.logger.warning("Could not append the '%s' attribute "
                                        "as it does not exist in the current image." % attrName)

        if not self.cfg.is_dummy_dataformat:
            self.deadpixelmap = np.append(self.deadpixelmap, img2.deadpixelmap[0:n_lines, :], axis=0)

        # NOTE: We leave the quicklook out here because merging the quicklook of adjacent scenes might cause a
        #       brightness jump that can be avoided by recomputing the quicklook after DN/radiance conversion.

    def DN2TOARadiance(self):
        """Convert DNs to TOA radiance.

        Convert the radiometric unit of _EnMAP_Detector_SensorGeo.data from digital numbers to top-of-atmosphere
        radiance.

        :return: None
        """
        # TODO move to processors.radiometric_transform?
        if self.detector_meta.unitcode == 'DN':
            self.logger.info('Converting DN values to radiance [mW/m^2/sr/nm] for %s detector...' % self.detector_name)

            if self.detector_meta.l_min is not None and self.detector_meta.l_max is not None:
                # Lλ = (LMINλ + ((LMAXλ - LMINλ)/(QCALMAX-QCALMIN)) * (QCAL-QCALMIN))
                # FIXME this asserts LMIN and LMAX in mW m^-2 sr^-1 nm^-1

                QCALMIN = 1
                QCALMAX = 2 ** 16  # 65535 (16-bit maximum value)
                LMIN = self.detector_meta.l_min
                LMAX = self.detector_meta.l_max
                QCAL = self.data[:]

                radiance = ((LMAX - LMIN)/(QCALMAX - QCALMIN)) * (QCAL - QCALMIN) + LMIN

            elif self.detector_meta.gains is not None and self.detector_meta.offsets is not None:
                # Lλ = QCAL * GAIN + OFFSET
                # NOTE: - DLR provides gains between 2000 and 10000, so we have to DEVIDE by gains
                #       - DLR gains / offsets are provided in W/m2/sr/nm, so we have to multiply by 1000 to get
                #         mW/m2/sr/nm as needed later
                radiance = 1000 * (self.data[:] * self.detector_meta.gains + self.detector_meta.offsets)

            else:
                raise ValueError("Neighter 'l_min'/'l_max' nor 'gains'/'offsets' "
                                 "are available for radiance computation.")

            self.data = radiance.astype(np.float32)
            self.detector_meta.unit = "mW m^-2 sr^-1 nm^-1"
            self.detector_meta.unitcode = "TOARad"
        else:
            self.logger.warning(
                "No DN to Radiance conversion is performed because unitcode is not DN (found: {code}).".format(
                    code=self.detector_meta.unitcode)
            )

    def save_raster_attributes(self, attrNames: List[str], outputDir: str):
        """Save the specified raster attributes to the given output directory.

        :param attrNames:   list of attribute names to be saved
        :param outputDir:   output directory
        """
        for attrName in attrNames:
            attr: GeoArray = getattr(self, attrName)
            fN = getattr(self.detector_meta, 'filename_%s' % attrName)

            if attr is not None:
                attr.save(path.join(outputDir, fN), fmt="GTiff")
            else:
                self.logger.warning("Could not save the %s attribute '%s' as it does not exist in the current image."
                                    % (self.detector_name, attrName))

    def _transform_raster_geometry_from_other_detector(self,
                                                       array: np.ndarray,
                                                       src_lons: np.ndarray,
                                                       src_lats: np.ndarray,
                                                       src_epsg: int,
                                                       resamp_alg: str = 'nearest',
                                                       respect_keystone: bool = False
                                                       ) -> np.ndarray:
        """Transform the given input raster from SWIR to VNIR or from SWIR to VNIR sensor geometry.

        NOTE:

        - The transformation target is always the EnMAP_Detector_SensorGeo instance sensor geometry
          (e.g., VNIR sensorgeo if self.detector_name == 'VNIR').
        - In case a 3D array is given and the array has the exact dimensions of the source detector,
          the full geolayer is only used if 'respect_keystone' is set to True. This saves computation time
          for input arrays where the keystone uncertainty does not matter.


        :param array:               input array to be transformed (2- or 3-dimensional)
        :param src_lons:            geolayer longitudes corresponding to the input array
                                    (same nbands like the source detector)
        :param src_lats:            geolayer latitudes corresponding to the input array
                                    (same nbands like the source detector)
        :param src_epsg:            projection EPSG code of the source array
        :param resamp_alg:          resampling algorithm ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case a 3D array
                                    in the dimension of the source detector is passed (default: False)
        :return:
        """
        detN = self.detector_name
        if self.detector_meta.lons is None or self.detector_meta.lats is None:
            raise RuntimeError(f"The {detN} geolayer must be computed first "
                               f"to transform arrays from {'SWIR' if detN == 'VNIR' else 'VNIR'} "
                               f"to {detN} sensor geometry.")

        if src_lons.shape != src_lats.shape:
            raise ValueError("'src_lons' must have the same shape as 'src_lats'.")

        vnir_lons = self.detector_meta.lons if detN == 'VNIR' else src_lons
        vnir_lats = self.detector_meta.lats if detN == 'VNIR' else src_lats
        swir_lons = self.detector_meta.lons if detN == 'SWIR' else src_lons
        swir_lats = self.detector_meta.lats if detN == 'SWIR' else src_lats
        prj_vnir = self.detector_meta.epsg_ortho if detN == 'VNIR' else src_epsg
        prj_swir = self.detector_meta.epsg_ortho if detN == 'SWIR' else src_epsg

        # use first geolayer band if the input array has only one band
        if array.ndim == 2 or \
           array.shape[2] != src_lons.shape[2] or \
           not respect_keystone:
            vnir_lons = vnir_lons[:, :, 0]
            vnir_lats = vnir_lats[:, :, 0]
            swir_lons = swir_lons[:, :, 0]
            swir_lats = swir_lats[:, :, 0]

        from ...processors.spatial_transform import VNIR_SWIR_SensorGeometryTransformer
        VS_SGT = VNIR_SWIR_SensorGeometryTransformer(lons_vnir=vnir_lons,
                                                     lats_vnir=vnir_lats,
                                                     lons_swir=swir_lons,
                                                     lats_swir=swir_lats,
                                                     prj_vnir=prj_vnir,
                                                     prj_swir=prj_swir,
                                                     res_vnir=(30, 30),
                                                     res_swir=(30, 30),
                                                     resamp_alg=resamp_alg,
                                                     # radius_of_influence=45,
                                                     nprocs=self.cfg.CPUs
                                                     )
        if detN == 'VNIR':
            return VS_SGT.transform_sensorgeo_SWIR_to_VNIR(array)
        else:
            return VS_SGT.transform_sensorgeo_VNIR_to_SWIR(array)


class EnMAP_VNIR_SensorGeo(EnMAP_Detector_SensorGeo):
    def __init__(self, root_dir: str, config: EnPTConfig, logger=None, meta=None) -> None:
        super().__init__(detector_name='VNIR', root_dir=root_dir, config=config, logger=logger, meta=meta)

    def read_masks(self):
        """Read the L1B masks."""
        self.logger.info('Reading image masks in VNIR sensor geometry.')

        # water mask (0=backgr.; 1=land; 2=water)
        if self.paths.mask_landwater:
            self.mask_landwater = GeoArray(self.paths.mask_landwater)[:]

        # cloud mask (0=none; 1=cloud)
        if self.paths.mask_clouds:
            self.mask_clouds = GeoArray(self.paths.mask_clouds)[:] == 1

        # cloud shadow mask (0=none; 1=cloud shadow)
        if self.paths.mask_cloudshadow:
            self.mask_cloudshadow = GeoArray(self.paths.mask_cloudshadow)[:] == 1

        # haze mask (0=none; 1=haze)
        if self.paths.mask_landwater:
            self.mask_haze = GeoArray(self.paths.mask_haze)[:] == 1

        # snow mask (0=none; 1=snow)
        if self.paths.mask_snow:
            self.mask_snow = GeoArray(self.paths.mask_snow)[:] == 1

        # cirrus mask (0=none; 1=thin, 2=medium, 3=thick)
        if self.paths.mask_cirrus:
            self.mask_cirrus = GeoArray(self.paths.mask_cirrus)[:]

    def transform_swir_to_vnir_raster(self,
                                      array_swirsensorgeo: np.ndarray,
                                      swir_lons: np.ndarray,
                                      swir_lats: np.ndarray,
                                      swir_epsg: int,
                                      resamp_alg: str = 'nearest',
                                      respect_keystone: bool = False
                                      ) -> np.ndarray:
        """Transform the given SWIR sensor-geometry raster array into VNIR sensor geometry.

        :param array_swirsensorgeo: source array in SWIR sensor geometry to be transformed
        :param swir_lons:           longitude geolayer array of the SWIR
        :param swir_lats:           latitude geolayer array of the SWIR
        :param swir_epsg:           EPSG code of the SWIR when transformed to map geometry
        :param resamp_alg:          resampling algorith ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case a 3D array
                                    in the dimension of the SWIR detector is passed (default: False)
        """
        return self._transform_raster_geometry_from_other_detector(
            array_swirsensorgeo, swir_lons, swir_lats, swir_epsg, resamp_alg, respect_keystone)


class EnMAP_SWIR_SensorGeo(EnMAP_Detector_SensorGeo):
    def __init__(self, root_dir: str, config: EnPTConfig, logger=None, meta=None) -> None:
        super().__init__(detector_name='SWIR', root_dir=root_dir, config=config, logger=logger, meta=meta)

    def __getattribute__(self, item):  # called whenever an instance attribute is accessed
        if item in ['mask_landwater', 'mask_clouds', 'mask_cloudshadow',
                    'mask_haze', 'mask_snow', 'mask_cirrus'] \
                and getattr(self, '_%s' % item) is None:
            self.logger.warning('The %s is not yet available in SWIR sensor geometry. '
                                'Use EnMAP_SWIR_SensorGeo.transform_vnir_to_swir_raster() to set it with a '
                                'transformed version of the one provided in VNIR sensor geometry.' % item)
        return super().__getattribute__(item)

    def transform_vnir_to_swir_raster(self,
                                      array_vnirsensorgeo: np.ndarray,
                                      vnir_lons: np.ndarray,
                                      vnir_lats: np.ndarray,
                                      vnir_epsg: int,
                                      resamp_alg: str = 'nearest',
                                      respect_keystone: bool = False
                                      ) -> np.ndarray:
        """Transform the given VNIR sensor-geometry raster array into SWIR sensor geometry.

        :param array_vnirsensorgeo: source array in VNIR sensor geometry to be transformed
        :param vnir_lons:           longitude geolayer array of the VNIR
        :param vnir_lats:           latitude geolayer array of the VNIR
        :param vnir_epsg:           EPSG code of the VNIR when transformed to map geometry
        :param resamp_alg:          resampling algorith ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case a 3D array
                                    in the dimension of the VNIR detector is passed (default: False)
        """
        return self._transform_raster_geometry_from_other_detector(
            array_vnirsensorgeo, vnir_lons, vnir_lats, vnir_epsg, resamp_alg, respect_keystone)


class EnMAPL1Product_SensorGeo(object):
    """Class for EnPT EnMAP object in sensor geometry.

    Attributes:
        - logger:
            - logging.Logger instance or subclass instance
        - vnir
            - instance of EnMAP_VNIR_SensorGeo class
        - swir
            - instance of EnMAP_SWIR_SensorGeo class
        - paths:
            - paths belonging to the EnMAP product
        - meta:
            - instance of EnMAP_Metadata_SensorGeo class
        - detector_attrNames:
            - list of attribute names for VNIR and SWIR detectors,

    """

    def __init__(self, root_dir: str, config: EnPTConfig, logger=None):
        """Get instance of EnPT EnMAP object in sensor geometry.

        :param root_dir:    Root directory of EnMAP Level-1B product
        :param config:      EnPT configuration object
        :param logger:      None or logging instance to be appended to EnMAPL1Product_SensorGeo instance
                            (If None, a default EnPT_Logger is used).
        """
        # protected attributes
        self._logger = None

        # populate attributes
        self.cfg = config
        if logger:
            self.logger = logger

        # Read metadata here in order to get all information needed by to get paths.
        if not self.cfg.is_dummy_dataformat:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*METADATA.XML"))[0],
                                                     config=self.cfg, logger=self.logger)
        else:
            self.meta = EnMAP_Metadata_L1B_SensorGeo(glob(path.join(root_dir, "*_header.xml"))[0],
                                                     config=self.cfg, logger=self.logger)
        self.meta.read_metadata()

        # define VNIR and SWIR detector
        self.detector_attrNames = ['vnir', 'swir']
        self.vnir = EnMAP_VNIR_SensorGeo(root_dir, config=self.cfg, logger=self.logger, meta=self.meta.vnir)
        self.swir = EnMAP_SWIR_SensorGeo(root_dir, config=self.cfg, logger=self.logger, meta=self.meta.swir)

        # Get the paths according information delivered in the metadata
        self.paths = self.get_paths()

        # associate raster attributes with file links (raster data is read lazily / on demand)
        # or directly read here in case the user does not want to include all L1B bands into the processing
        self.vnir.data = self.paths.vnir.data
        self.swir.data = self.paths.swir.data
        self.vnir.read_masks()

        if self.cfg.drop_bad_bands:
            self.vnir.data = self.vnir.data[:, :, self.meta.vnir.goodbands_inds]
            self.swir.data = self.swir.data[:, :, self.meta.swir.goodbands_inds]

        try:
            vnir_dpm = GeoArray(self.paths.vnir.deadpixelmap)
            swir_dpm = GeoArray(self.paths.swir.deadpixelmap)

            if self.cfg.drop_bad_bands:
                if vnir_dpm.ndim == 3:
                    self.vnir.deadpixelmap = vnir_dpm[:, :, self.meta.vnir.goodbands_inds]
                    self.swir.deadpixelmap = swir_dpm[:, :, self.meta.swir.goodbands_inds]
                else:  # 2D
                    self.vnir.deadpixelmap = vnir_dpm[self.meta.vnir.goodbands_inds, :]
                    self.swir.deadpixelmap = swir_dpm[self.meta.swir.goodbands_inds, :]
            else:
                self.vnir.deadpixelmap = vnir_dpm
                self.swir.deadpixelmap = swir_dpm

        except ValueError:
            self.logger.warning("Unexpected dimensions of dead pixel mask. Setting all pixels to 'normal'.")
            self.vnir.deadpixelmap = np.zeros(self.vnir.data.shape)
            self.swir.deadpixelmap = np.zeros(self.swir.data.shape)

        # NOTE: We leave the quicklook out here because merging the quicklook of adjacent scenes might cause a
        #       brightness jump that can be avoided by recomputing the quicklook after DN/radiance conversion.

        # compute radiance
        self.DN2TOARadiance()

    def get_paths(self) -> SimpleNamespace:
        """Get all file paths associated with the current instance of EnMAPL1Product_SensorGeo.

        :return: paths.SimpleNamespace()
        """
        paths = SimpleNamespace()
        paths.vnir = self.vnir.get_paths()
        paths.swir = self.swir.get_paths()
        paths.root_dir = self.meta.rootdir
        paths.metaxml = self.meta.path_xml

        return paths

    @property
    def logger(self) -> EnPT_Logger:
        """Get a an instance of enpt.utils.logging.EnPT_Logger.

        NOTE:
            - The logging level will be set according to the user inputs of EnPT.
            - The path of the log file is directly derived from the attributes of the _EnMAP_Image instance.

        Usage:
            - get the logger:
                logger = self.logger
            - set the logger
                self.logger = logging.getLogger()  # NOTE: only instances of logging.Logger are allowed here
            - delete the logger:
                del self.logger  # or "self.logger = None"

        :return: EnPT_Logger instance
        """
        if self._logger and self._logger.handlers[:]:
            return self._logger
        else:
            basename = path.splitext(path.basename(self.cfg.path_l1b_enmap_image))[0]
            path_logfile = path.join(self.cfg.output_dir, basename + '.log') \
                if self.cfg.create_logfile and self.cfg.output_dir else ''
            self._logger = EnPT_Logger('log__' + basename, fmt_suffix=None, path_logfile=path_logfile,
                                       log_level=self.cfg.log_level, append=False)

            return self._logger

    @logger.setter
    def logger(self, logger: logging.Logger):
        assert isinstance(logger, logging.Logger) or logger in ['not set', None], \
            "%s.logger can not be set to %s." % (self.__class__.__name__, logger)

        # save prior logs
        # if logger is None and self._logger is not None:
        #     self.log += self.logger.captured_stream
        self._logger = logger

    @property
    def log(self) -> str:
        """Return a string of all logged messages until now.

        NOTE: self.log can also be set to a string.
        """
        return self.logger.captured_stream

    @log.setter
    def log(self, string: str):
        assert isinstance(string, str), "'log' can only be set to a string. Got %s." % type(string)
        self.logger.captured_stream = string

    # @classmethod
    # def from_L1B_provider_data(cls, path_enmap_image: str, config: EnPTConfig=None) -> EnMAPL1Product_SensorGeo:
    #     """
    #
    #     :param path_enmap_image:
    #     :param config:
    #     :return:
    #     """
    #     # input validation
    #     if not path.isdir(path_enmap_image) and \
    #        not (path.exists(path_enmap_image) and path_enmap_image.endswith('.zip')):
    #         raise ValueError("The parameter 'path_enmap_image' must be a directory or the path to an existing zip "
    #                          "archive.")
    #
    #     # extract L1B image archive if needed
    #     if path_enmap_image.endswith('.zip'):
    #         path_enmap_image = self.extract_zip_archive(path_enmap_image)
    #         if not path.isdir(path_enmap_image):
    #             raise NotADirectoryError(path_enmap_image)
    #
    #     # run the reader
    #     from ..io.reader import L1B_Reader
    #     RD = L1B_Reader(config=config)
    #     L1_obj = RD.read_inputdata(path_enmap_image, observation_time=datetime(2015, 12, 7, 10))
    #
    #     return L1_obj

    def DN2TOARadiance(self):
        """Convert self.vnir.data and self.swir.data from digital numbers to top-of-atmosphere radiance.

        :return: None
        """
        if self.vnir.detector_meta.unitcode != 'TOARad':
            self.vnir.DN2TOARadiance()
            self.meta.vnir.unitcode = self.vnir.detector_meta.unitcode  # FIXME possible duplicate?
            self.meta.vnir.unit = self.vnir.detector_meta.unit  # FIXME possible duplicate?
        if self.swir.detector_meta.unitcode != 'TOARad':
            self.swir.DN2TOARadiance()
            self.meta.swir.unitcode = self.swir.detector_meta.unitcode
            self.meta.swir.unit = self.swir.detector_meta.unit

    def append_new_image(self, root_dir, n_line_ext):
        """Append a second EnMAP Level-1B product below the last line of the current L1B product.

        NOTE:   We create new files that will be saved into a temporary directory and their path will be stored in the
                paths of self. We assume that the dead pixel map does not change between two adjacent images.

        :param root_dir:    root directory of EnMAP Level-1B product to be appended
        :param n_line_ext:  number of lines to be added from the new image
        :return:
        """
        l1b_ext_obj = EnMAPL1Product_SensorGeo(root_dir, config=self.cfg)

        self.vnir.append_new_image(l1b_ext_obj.vnir, n_line_ext)
        self.swir.append_new_image(l1b_ext_obj.swir, n_line_ext)

    def calc_snr_from_radiance(self):
        """Compute EnMAP SNR from radiance data.

        SNR equation:    SNR = p0 + p1 * LTOA + p2 * LTOA ^ 2   [W/(m^2 sr nm)].
        """
        with TemporaryDirectory(dir=self.cfg.working_dir) as tmpDir, \
                ZipFile(self.cfg.path_l1b_snr_model, "r") as zf:

            zf.extractall(tmpDir)

            if self.meta.vnir.unitcode == 'TOARad':
                self.vnir.detector_meta.calc_snr_from_radiance(rad_data=self.vnir.data, dir_snr_models=tmpDir)

            if self.meta.swir.unitcode == 'TOARad':
                self.swir.detector_meta.calc_snr_from_radiance(rad_data=self.swir.data, dir_snr_models=tmpDir)

    def correct_dead_pixels(self):
        """Correct dead pixels of both detectors."""
        self.vnir.correct_dead_pixels()
        self.swir.correct_dead_pixels()

    # def correct_VNIR_SWIR_shift(self):
    #     # use first geolayer bands for VNIR and SWIR
    #     Vlons, Vlats = self.vnir.detector_meta.lons[:, :, 0], self.vnir.detector_meta.lats[:, :, 0]
    #     Slons, Slats = self.swir.detector_meta.lons[:, :, 0], self.swir.detector_meta.lats[:, :, 0]
    #
    #     # get corner coordinates of VNIR and SWIR according to geolayer
    #     def get_coords(lons, lats):
    #         return tuple([(lons[Y, X], lats[Y, X]) for Y, X in [(0, 0), (0, -1), (-1, 0), (-1, -1)]])
    #
    #     VUL, VUR, VLL, VLR = get_coords(Vlons, Vlats)
    #     SUL, SUR, SLL, SLR = get_coords(Slons, Slats)
    #
    #     # get map coordinates of VNIR/SWIR overlap
    #     ovUL = max(VUL[0], SUL[0]), min(VUL[1], SUL[1])
    #     ovUR = min(VUR[0], SUR[0]), min(VUR[1], SUR[1])
    #     ovLL = max(VLL[0], SLL[0]), max(VLL[1], SLL[1])
    #     ovLR = min(VLR[0], SLR[0]), max(VLR[1], SLR[1])
    #
    #     # find nearest image positions for VNIR and SWIR to the map coordinates of the VNIR/SWIR overlap
    #     def nearest_imCoord(lons_arr, lats_arr, lon, lat):
    #         dists = np.sqrt((lons_arr - lon) ** 2 + (lats_arr - lat) ** 2)
    #         row, col = np.unravel_index(dists.argmin(), dists.shape)
    #
    #         return row, col
    #
    #     overlapImVNIR = tuple([nearest_imCoord(Vlons, Vlats, *ovCoords) for ovCoords in [ovUL, ovUR, ovLL, ovLR]])
    #     overlapImSWIR = tuple([nearest_imCoord(Slons, Slats, *ovCoords) for ovCoords in [ovUL, ovUR, ovLL, ovLR]])
    #
    #     raise NotImplementedError()  # FIXME
    #     self.vnir.data.get_mapPos()  # FIXME

    def get_preprocessed_dem(self):
        self.vnir.get_preprocessed_dem()
        self.swir.get_preprocessed_dem()

    def transform_vnir_to_swir_raster(self,
                                      array_vnirsensorgeo: np.ndarray,
                                      resamp_alg: str = 'nearest',
                                      respect_keystone: bool = False
                                      ) -> np.ndarray:
        """Transform the given array from VNIR into SWIR sensor geometry.

        :param array_vnirsensorgeo: raster array in VNIR sensor geometry to be transformed into SWIR sensor geometry
        :param resamp_alg:          resampling algorithm ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case a 3D array
                                    in the dimension of the VNIR detector is passed (default: False)
        """
        if self.meta.vnir.lons is None or self.meta.vnir.lats is None or \
           self.meta.swir.lons is None or self.meta.swir.lats is None:
            raise RuntimeError('The VNIR/SWIR geolayers must be computed first '
                               'to transform arrays from VNIR to SWIR sensor geometry.')

        return self.swir.transform_vnir_to_swir_raster(array_vnirsensorgeo=array_vnirsensorgeo,
                                                       vnir_lons=self.meta.vnir.lons,
                                                       vnir_lats=self.meta.vnir.lats,
                                                       vnir_epsg=self.meta.vnir.epsg_ortho,
                                                       resamp_alg=resamp_alg,
                                                       respect_keystone=respect_keystone)

    def transform_swir_to_vnir_raster(self, array_swirsensorgeo: np.ndarray,
                                      resamp_alg: str = 'nearest',
                                      respect_keystone: bool = False
                                      ) -> np.ndarray:
        """Transform the given array from SWIR into VNIR sensor geometry.

        :param array_swirsensorgeo: raster array in SWIR sensor geometry to be transformed into VNIR sensor geometry
        :param resamp_alg:          resampling algorithm ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case a 3D array
                                    in the dimension of the VNIR detector is passed (default: False)
        """
        if self.meta.vnir.lons is None or self.meta.vnir.lats is None or \
           self.meta.swir.lons is None or self.meta.swir.lats is None:
            raise RuntimeError('The VNIR/SWIR geolayers must be computed first '
                               'to transform arrays from VNIR to SWIR sensor geometry.')

        return self.vnir.transform_swir_to_vnir_raster(array_swirsensorgeo=array_swirsensorgeo,
                                                       swir_lons=self.meta.swir.lons,
                                                       swir_lats=self.meta.swir.lats,
                                                       swir_epsg=self.meta.swir.epsg_ortho,
                                                       resamp_alg=resamp_alg,
                                                       respect_keystone=respect_keystone)

    def set_SWIRattr_with_transformedVNIRattr(self, attrName: str,
                                              resamp_alg: str = 'nearest',
                                              respect_keystone: bool = False
                                              ) -> None:
        """Set the specified SWIR raster attribute with a VNIR attribute transformed to SWIR sensor geometry.

        :param attrName:            name of the attribute to be set
        :param resamp_alg:          resampling algorithm ('nearest', 'bilinear', 'gauss', 'custom')
        :param respect_keystone:    whether to use the full geoarray (all bands) in case the attribute
                                    to be transformed is 'data' (default: False)
        """
        self.logger.info("Transforming the '%s' attribute from VNIR to SWIR sensor geometry." % attrName)

        vnir_rasterAttr = getattr(self.vnir, attrName)

        if vnir_rasterAttr is None:
            raise RuntimeError("%s.vnir.%s has not yet been set." % (self.__class__.__name__, attrName))

        attr_transformed = self.transform_vnir_to_swir_raster(array_vnirsensorgeo=np.array(vnir_rasterAttr),
                                                              resamp_alg=resamp_alg,
                                                              respect_keystone=respect_keystone)
        setattr(self.swir, attrName, attr_transformed)

    def run_AC(self):
        from ...processors.atmospheric_correction import AtmosphericCorrector
        AC = AtmosphericCorrector(config=self.cfg)
        AC.run_ac(self)

    def save(self, outdir: str, suffix="") -> str:
        """Save the product to disk using almost the same input format.

        NOTE: Radiance is saved instead of DNs to avoid a radiometric jump between concatenated images.

        :param outdir: path to the output directory
        :param suffix: suffix to be appended to the output filename (???)
        :return: root path (root directory) where products were written
        """
        product_dir = path.join(path.abspath(outdir),
                                "{name}{suffix}".format(name=self.meta.scene_basename, suffix=suffix))

        self.logger.info("Write product to: %s" % product_dir)
        makedirs(product_dir, exist_ok=True)

        # write the VNIR
        self.vnir.save_raster_attributes(['data', 'mask_landwater', 'mask_clouds', 'mask_cloudshadow',
                                          'mask_haze', 'mask_snow', 'mask_cirrus', 'deadpixelmap'],
                                         product_dir)

        # FIXME we could also write the quicklook included in DLR L1B format instead of generating an own one
        self.vnir.generate_quicklook(bands2use=self.vnir.detector_meta.preview_bands) \
            .save(path.join(product_dir, path.splitext(self.meta.vnir.filename_quicklook)[0] + '.png'), fmt='PNG')

        # write the SWIR
        # NOTE: The masks only exist in VNIR sensor geometry (would have to be transformed before to match for the SWIR)
        self.swir.save_raster_attributes(['data', 'deadpixelmap'], product_dir)

        self.swir.generate_quicklook(bands2use=self.swir.detector_meta.preview_bands) \
            .save(path.join(product_dir, path.splitext(self.meta.swir.filename_quicklook)[0] + '.png'), fmt='PNG')

        # Update the xml file
        metadata_string = self.meta.to_XML()
        with open(product_dir + path.sep + path.basename(self.paths.metaxml), "w") as xml_file:
            xml_file.write(metadata_string)

        self.logger.info("L1B product successfully written!")

        return product_dir
