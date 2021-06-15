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

"""EnPT metadata objects for EnMAP data in sensor geometry."""

from datetime import datetime
from lxml import etree as ElementTree
import logging
import os
import fnmatch
from typing import Union, List, Tuple, Optional  # noqa: F401
from collections import OrderedDict
import numpy as np
from py_tools_ds.geo.vector.topology import Polygon, get_footprint_polygon  # noqa: F401  # flake8 issue
from geoarray import GeoArray

from ...options.config import EnPTConfig
from ..srf import SRF
from ...processors.spatial_transform import RPC_3D_Geolayer_Generator

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


class EnMAP_Metadata_L1B_Detector_SensorGeo(object):
    """Class for all EnMAP metadata associated with a single EnMAP detector in sensor geometry.

    NOTE:
        - All metadata that have VNIR and SWIR detector in sensor geometry in common should be included here.

    """

    def __init__(self, detector_name: str, config: EnPTConfig, logger: logging.Logger = None):
        """Get a metadata object containing the metadata of a single EnMAP detector in sensor geometry.

        :param detector_name:   Name of the detector (VNIR or SWIR)
        :param config:          EnPT configuration object
        :param logger:          instance of logging.logger or subclassed
        """
        from . import L1B_product_props, L1B_product_props_DLR
        self.cfg = config
        self.detector_name: str = detector_name
        if not self.cfg.is_dummy_dataformat:
            self.detector_label = L1B_product_props_DLR['xml_detector_label'][detector_name]
        else:
            self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.logger = logger or logging.getLogger()

        # These lines are used to load path information
        self.filename_data: Optional[str] = ''  # detector data filename
        self.scene_basename: Optional[str] = ''  # basename of the EnMAP image
        self.filename_deadpixelmap: Optional[str] = ''  # filename of the dead pixel file
        self.filename_quicklook: Optional[str] = ''  # filename of the quicklook file
        # FIXME masks of BOTH detectors
        self.filename_mask_landwater: Optional[str] = ''  # filename of the land/water mask file
        self.filename_mask_snow: Optional[str] = ''  # filename of the snow mask file
        self.filename_mask_cloudshadow: Optional[str] = ''  # filename of the cloud shadow mask file
        self.filename_mask_clouds: Optional[str] = ''  # filename of the cloud mask file
        self.filename_mask_haze: Optional[str] = ''  # filename of the haze mask file
        self.filename_mask_cirrus: Optional[str] = ''  # filename of the cirrus mask file

        self.wvl_center: Optional[np.ndarray] = None  # Center wavelengths for each EnMAP band
        self.fwhm: Optional[np.ndarray] = None  # Full width half maximmum for each EnMAP band
        self.srf: Optional[SRF] = None  # SRF object holding the spectral response functions for each EnMAP band
        self.solar_irrad: Optional[np.ndarray] = None  # solar irradiance in [W/m2/nm] for each band
        self.nwvl: Optional[int] = None  # Number of wave bands
        self.nrows: Optional[int] = None  # number of rows
        self.ncols: Optional[int] = None  # number of columns
        self.smile_coef: Optional[np.ndarray] = None  # smile coefficients needed for smile computation
        self.nsmile_coef: Optional[int] = None  # number of smile coefficients
        self.smile: Optional[np.ndarray] = None  # smile for each EnMAP image column
        self.gains: Optional[np.ndarray] = None  # band-wise gains for computing radiance from DNs
        self.offsets: Optional[np.ndarray] = None   # band-wise offsets for computing radiance from DNs
        self.l_min: Optional[np.ndarray] = None  # band-wise l-min for computing radiance from DNs
        self.l_max: Optional[np.ndarray] = None  # band-wise l-max for computing radiance from DNs
        self.goodbands_inds: Optional[List] = None  # list of band indices included in the processing (all other bands are removed)  # noqa
        self.lat_UL_UR_LL_LR: Optional[List[float, float, float, float]] = None  # latitude coords for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR: Optional[List[float, float, float, float]] = None  # longitude coords for UL, UR, LL, LR
        self.epsg_ortho: Optional[int] = None  # EPSG code of the orthorectified image
        self.rpc_coeffs: OrderedDict = OrderedDict()  # RPC coefficients for geolayer computation
        self.ll_mapPoly: Optional[Polygon] = None  # footprint polygon in longitude/latitude map coordinates
        self.lats: Optional[np.ndarray] = None  # 2D array of latitude coordinates according to given lon/lat sampling
        self.lons: Optional[np.ndarray] = None  # 2D array of longitude coordinates according to given lon/lat sampling
        self.unit: str = ''  # radiometric unit of pixel values
        self.unitcode: str = ''  # code of radiometric unit
        self.preview_wvls: Optional[List[float]] = None  # wavelengths to be used for quicklook images
        self.preview_bands: Optional[List[int]] = None  # band indices to be used for quicklook images
        self.snr: Optional[np.ndarray] = None   # Signal to noise ratio as computed from radiance data

    def read_metadata(self, path_xml):
        """Read the metadata of a specific EnMAP detector in sensor geometry.

        :param path_xml: file path of the metadata file
        :return: None
        """
        xml = ElementTree.parse(path_xml).getroot()

        if not self.cfg.is_dummy_dataformat:
            lbl = self.detector_label
            self.logger.info("Reading metadata for %s detector..." % self.detector_name)

            # read data filenames
            all_filenames = [ele.text for ele in xml.findall("product/productFileInformation/file/name")]

            def get_filename(matching_exp: str):
                matches = fnmatch.filter(all_filenames, '%s.GEOTIFF' % matching_exp)
                if not matches:
                    matches = fnmatch.filter(all_filenames, '%s.TIF' % matching_exp)
                if not matches:
                    raise FileNotFoundError("Could not find a file that matches the expression '%s'." % matching_exp)

                return matches[0]

            self.filename_data = xml.find("product/image/%s/name" % lbl).text
            self.scene_basename = self.filename_data.split('-SPECTRAL_IMAGE')[0]
            self.filename_quicklook = xml.find("product/quicklook/%s/name" % lbl).text
            self.filename_deadpixelmap = get_filename('*QL_PIXELMASK_%s' % self.detector_name)
            self.filename_mask_landwater = get_filename('*QL_QUALITY_CLASSES')
            self.filename_mask_snow = get_filename('*QL_QUALITY_SNOW')
            self.filename_mask_cloudshadow = get_filename('*QL_QUALITY_CLOUDSHADOW')
            self.filename_mask_clouds = get_filename('*-QL_QUALITY_CLOUD')
            self.filename_mask_haze = get_filename('*QL_QUALITY_HAZE')
            self.filename_mask_cirrus = get_filename('*QL_QUALITY_CIRRUS')

            # FIXME combine different cloud masks?
            # TODO: Add test flags layer.

            # read some basic information concerning the detector
            self.nrows = int(xml.find("product/image/%s/dimension/rows" % lbl).text)
            self.ncols = int(xml.find("product/image/%s/dimension/columns" % lbl).text)
            self.unitcode = 'DN'
            self.unit = ''

            # Read image coordinates
            # FIXME why do we have the same corner coordinates for VNIR and SWIR?
            points = xml.findall("base/spatialCoverage/boundingPolygon/point")
            coords_xy = {p.find('frame').text: (float(p.find('longitude').text),
                                                float(p.find('latitude').text))
                         for p in points}
            coord_tags = ['upper_left', 'upper_right', 'lower_left', 'lower_right']
            self.lon_UL_UR_LL_LR = [coords_xy[ct][0] for ct in coord_tags]
            self.lat_UL_UR_LL_LR = [coords_xy[ct][1] for ct in coord_tags]

            # read the band related information: wavelength, fwhm
            self.nwvl = int(xml.find("product/image/%s/channels" % lbl).text)
            # FIXME hardcoded - DLR does not provide any smile information
            #   => smile coefficients are set to zero until we have valid ones
            self.nsmile_coef = 5
            self.smile_coef = np.zeros((self.nwvl, self.nsmile_coef), dtype=float)

            startband = 0 if self.detector_name == 'VNIR' else int(xml.find("product/image/vnir/channels").text)
            subset = slice(startband, startband + self.nwvl)
            bi = "specific/bandCharacterisation/bandID/"
            self.wvl_center = np.array([float(ele.text) for ele in xml.findall(bi + "wavelengthCenterOfBand")[subset]])
            self.fwhm = np.array([float(ele.text) for ele in xml.findall(bi + "FWHMOfBand")[subset]])
            self.gains = np.array([float(ele.text) for ele in xml.findall(bi + "GainOfBand")[subset]])
            self.offsets = np.array([float(ele.text) for ele in xml.findall(bi + "OffsetOfBand")[subset]])

            # read preview bands
            wvl_red = float(xml.find("product/image/%s/qlChannels/red" % lbl).text)
            wvl_green = float(xml.find("product/image/%s/qlChannels/green" % lbl).text)
            wvl_blue = float(xml.find("product/image/%s/qlChannels/blue" % lbl).text)
            self.preview_wvls = [wvl_red, wvl_green, wvl_blue]
            self.preview_bands = np.array([np.argmin(np.abs(self.wvl_center - wvl)) for wvl in self.preview_wvls])

            # read RPC coefficients
            for bID in xml.findall("product/navigation/RPC/bandID")[subset]:
                bN = 'band_%d' % np.int64(bID.attrib['number'])

                keys2combine = ('row_num', 'row_den', 'col_num', 'col_den')

                tmp = OrderedDict([(ele.tag.lower(), float(ele.text)) for ele in bID.findall('./')])
                self.rpc_coeffs[bN] = {k: v for k, v in tmp.items() if not k.startswith(keys2combine)}

                for n in keys2combine:
                    self.rpc_coeffs[bN]['%s_coeffs' % n.lower()] = \
                        np.array([v for k, v in tmp.items() if k.startswith(n)])

        else:
            lbl = self.detector_label
            self.logger.info("Reading metadata for %s detector..." % self.detector_name)

            # read data filenames
            self.filename_data = xml.findall("ProductComponent/%s/Data/Filename" % lbl)[0].text
            self.scene_basename = os.path.splitext(self.filename_data)[0]
            self.filename_deadpixelmap = xml.findall("ProductComponent/%s/Sensor/DeadPixel/Filename" % lbl)[0].text
            self.filename_quicklook = xml.findall("ProductComponent/%s/Preview/Filename" % lbl)[0].text
            self.filename_mask_clouds = xml.findall("ProductComponent/%s/Data/CloudMaskMap/Filename" % lbl)[0].text

            # read preview bands
            self.preview_bands = np.zeros(3, dtype=int)
            self.preview_bands[0] = int(xml.findall("ProductComponent/%s/Preview/Bands/Red" % lbl)[0].text)
            self.preview_bands[1] = int(xml.findall("ProductComponent/%s/Preview/Bands/Green" % lbl)[0].text)
            self.preview_bands[2] = int(xml.findall("ProductComponent/%s/Preview/Bands/Blue" % lbl)[0].text)

            # read some basic information concerning the detector
            self.nrows = int(xml.findall("ProductComponent/%s/Data/Size/NRows" % lbl)[0].text)
            self.ncols = int(xml.findall("ProductComponent/%s/Data/Size/NCols" % lbl)[0].text)
            self.unitcode = xml.findall("ProductComponent/%s/Data/Type/UnitCode" % lbl)[0].text
            self.unit = xml.findall("ProductComponent/%s/Data/Type/Unit" % lbl)[0].text

            # Read image coordinates
            scene_corner_coordinates = xml.findall("ProductComponent/%s/Data/SceneInformation/"
                                                   "SceneCornerCoordinates" % lbl)
            self.lat_UL_UR_LL_LR = [
                float(scene_corner_coordinates[0].findall("Latitude")[0].text),
                float(scene_corner_coordinates[1].findall("Latitude")[0].text),
                float(scene_corner_coordinates[2].findall("Latitude")[0].text),
                float(scene_corner_coordinates[3].findall("Latitude")[0].text)
            ]
            self.lon_UL_UR_LL_LR = [
                float(scene_corner_coordinates[0].findall("Longitude")[0].text),
                float(scene_corner_coordinates[1].findall("Longitude")[0].text),
                float(scene_corner_coordinates[2].findall("Longitude")[0].text),
                float(scene_corner_coordinates[3].findall("Longitude")[0].text)
            ]

            # read the band related information: wavelength, fwhm
            self.nwvl = int(xml.findall("ProductComponent/%s/Data/BandInformationList/NumberOfBands" % lbl)[0].text)
            self.nsmile_coef = int(xml.findall(
                "ProductComponent/%s/Data/BandInformationList/SmileInformation/NumberOfCoefficients" % lbl)[0].text)
            self.fwhm = np.zeros(self.nwvl, dtype=float)
            self.wvl_center = np.zeros(self.nwvl, dtype=float)
            self.smile_coef = np.zeros((self.nwvl, self.nsmile_coef), dtype=float)
            self.l_min = np.zeros(self.nwvl, dtype=float)
            self.l_max = np.zeros(self.nwvl, dtype=float)
            band_informations = xml.findall("ProductComponent/%s/Data/BandInformationList/BandInformation" % lbl)
            for bi in band_informations:
                k = np.int64(bi.attrib['Id']) - 1
                self.wvl_center[k] = float(bi.findall("CenterWavelength")[0].text)
                self.fwhm[k] = float(bi.findall("FullWidthHalfMaximum")[0].text)
                self.l_min[k] = float(bi.findall("L_min")[0].text)
                self.l_max[k] = float(bi.findall("L_max")[0].text)
                scl = bi.findall("Smile/Coefficient")
                for sc in scl:
                    self.smile_coef[k, np.int64(sc.attrib['exponent'])] = float(sc.text)

            self.lats = self.interpolate_corners(*self.lat_UL_UR_LL_LR, self.ncols, self.nrows)
            self.lons = self.interpolate_corners(*self.lon_UL_UR_LL_LR, self.ncols, self.nrows)
            self.preview_wvls = np.array([self.wvl_center[i] for i in self.preview_bands])

        # drop bad bands from metadata if desired
        if self.cfg.drop_bad_bands:
            wvls = list(self.wvl_center)
            self.goodbands_inds = gbl = [wvls.index(wvl) for wvl in wvls
                                         if not (1358 < wvl < 1453 or
                                                 1814 < wvl < 1961)]

            if len(gbl) < len(wvls):
                for attrN in ['wvl_center', 'fwhm', 'offsets', 'gains', 'l_min', 'l_max']:
                    if getattr(self, attrN) is not None:
                        setattr(self, attrN, getattr(self, attrN)[gbl])

                self.nwvl = len(gbl)
                self.smile_coef = self.smile_coef[gbl, :]
                self.rpc_coeffs = OrderedDict([(k, v) for i, (k, v) in enumerate(self.rpc_coeffs.items())
                                               if i in gbl])

        # compute metadata derived from read data
        self.smile = self.calc_smile()
        self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
        self.solar_irrad = self.calc_solar_irradiance_CWL_FWHM_per_band()
        self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                          self.lat_UL_UR_LL_LR)), fix_invalid=True)
        from ...processors.spatial_transform import get_UTMEPSG_from_LonLat_cornersXY
        # NOTE:   self.cfg.target_epsg is set if user-provided or in case of Lon/Lat coordinates
        self.epsg_ortho = \
            self.cfg.target_epsg or \
            get_UTMEPSG_from_LonLat_cornersXY(lons=self.lon_UL_UR_LL_LR, lats=self.lat_UL_UR_LL_LR)

    def calc_smile(self):
        """Compute smile for each EnMAP column.

        The sum in (1) is expressed as inner product of two arrays with inner dimension as the polynomial smile
        coefficients shape = (ncols, nsmile_coef) of polynomial x

        :return:
        """
        # smile(icol, iwvl) = sum_(p=0)^(nsmile_coef-1) smile_coef[iwvl, p] * icol**p (1)
        return np.inner(
            np.array([[icol ** p for p in np.arange(self.nsmile_coef)] for icol in np.arange(self.ncols)]),
            self.smile_coef  # shape = (nwvl, nsmile_coef)
        )  # shape = (ncols, nwvl)

    def calc_snr_from_radiance(self, rad_data: Union[GeoArray, np.ndarray], dir_snr_models: str):
        """Compute EnMAP SNR from radiance data for the given detector.

        SNR equation:    SNR = p0 + p1 * LTOA + p2 * LTOA ^ 2   [W/(m^2 sr nm)].

        NOTE:   The SNR model files (SNR_D1.bsq/SNR_D2.bsq) contain polynomial coefficients needed to compute SNR.

                SNR_D1.bsq: SNR model for EnMAP SWIR detector (contains high gain and low gain model coefficients)

                - 1000 columns (for 1000 EnMAP columns)
                - 88 bands for 88 EnMAP VNIR bands
                - 7 lines:  - 3 lines: high gain coefficients
                            - 3 lines: low gain coefficients
                            - 1 line: threshold needed to decide about high gain or low gain

                SNR_D2.bsq: SNR model for EnMAP SWIR detector

                - 1000 columns (for 1000 EnMAP columns)
                - x bands for x EnMAP SWIR bands
                - 3 lines for 3 coefficients

        :param rad_data:        image radiance data of EnMAP_Detector_SensorGeo
        :param dir_snr_models:  root directory where SNR model data is stored (must contain SNR_D1.bsq/SNR_D2.bsq)
        """
        path_snr_model = os.path.join(dir_snr_models, "SNR_D1.bsq" if self.detector_name == 'VNIR' else "SNR_D2.bsq")
        rad_data = np.array(rad_data)

        assert self.unitcode == 'TOARad'
        self.logger.info("Computing SNR for %s using %s" % (self.detector_name, path_snr_model))

        if self.detector_name == 'VNIR':
            gA = GeoArray(path_snr_model)
            coeffs_highgain = gA[0:3, :, :]  # [3 x ncols x nbands]
            coeffs_lowgain = gA[3:6, :, :]  # [3 x ncols x nbands]
            gain_threshold = np.squeeze(gA[6, :, :])

            self.snr = np.zeros(rad_data.shape)
            for irow in range(self.nrows):
                highgain_mask = rad_data[irow, :, :] < gain_threshold  # a single row
                rad_highgain = rad_data[irow, highgain_mask]
                self.snr[irow, :, :][highgain_mask] = \
                    coeffs_highgain[0, highgain_mask] + \
                    coeffs_highgain[1, highgain_mask] * rad_highgain + \
                    coeffs_highgain[2, highgain_mask] * rad_highgain ** 2

                lowgain_mask = rad_data[irow, :, :] >= gain_threshold
                rad_lowgain = rad_data[irow, lowgain_mask]
                self.snr[irow, :, :][lowgain_mask] = \
                    coeffs_lowgain[0, lowgain_mask] + \
                    coeffs_lowgain[1, lowgain_mask] * rad_lowgain + \
                    coeffs_lowgain[2, lowgain_mask] * rad_lowgain ** 2

        else:
            gA = GeoArray(path_snr_model)

            # some SWIR bands may be missing -> find indices of bands we need for SNR computation
            cwls_gA = gA.metadata.band_meta['wavelength']
            cwls_needed = self.wvl_center

            idxs_needed = [np.argmin(np.abs(cwls_gA - cwl)) for cwl in cwls_needed]
            if not len(set(idxs_needed)) == len(idxs_needed) or len(idxs_needed) != self.nwvl:
                raise RuntimeError('Unclear band allocation during SWIR SNR computation.')

            coeffs = gA[:, :, idxs_needed]  # [3 x ncols x nbands]
            self.snr = coeffs[0, :, :] + coeffs[1, :, :] * rad_data[:, :, :] + coeffs[2, :, :] * rad_data[:, :, :] ** 2

    @staticmethod
    def interpolate_corners(ul: float, ur: float, ll: float, lr: float, nx: int, ny: int):
        """Compute interpolated field from corner values of a scalar field given at: ul, ur, ll, lr.

        :param ul:  tbd
        :param ur:  tbd
        :param ll:  tbd
        :param lr:  tbd
        :param nx: final shape (x-axis direction)
        :param ny: final shape (y-axis direction)
        """
        # FIXME this method must later be replaced by the geolayer provided by the ground segment
        #       => a linear interpolation between the EnMAP corner coordinates is NOT sufficient for modelling the
        #          geometry of VNIR and SWIR
        #       - especially at off-nadir acquisitions and with some terrain present, a linear interpolation leads
        #         to large deviations (> 180 m y-coordinate offset for the EnPT test dataset)

        # TODO ensure that lons/lats represent UL coordinates not pixel coordinate centers (as given by Karl / DLR(?))

        corner_coords = np.array([[ul, ur],
                                  [ll, lr]])
        rowpos, colpos = [0, 1], [0, 1]

        from scipy.interpolate import RegularGridInterpolator
        rgi = RegularGridInterpolator([rowpos, colpos], corner_coords, method='linear')
        out_rows_grid, out_cols_grid = np.meshgrid(np.linspace(0, 1, ny),
                                                   np.linspace(0, 1, nx),
                                                   indexing='ij')

        coords = rgi(np.dstack([out_rows_grid, out_cols_grid]))

        return coords

    def compute_geolayer_for_cube(self):
        self.logger.info('Computing %s geolayer...' % self.detector_name)
        lons, lats = \
            RPC_3D_Geolayer_Generator(rpc_coeffs_per_band=self.rpc_coeffs,
                                      elevation=self.cfg.path_dem if self.cfg.path_dem else self.cfg.average_elevation,
                                      enmapIm_cornerCoords=tuple(zip(self.lon_UL_UR_LL_LR, self.lat_UL_UR_LL_LR)),
                                      enmapIm_dims_sensorgeo=(self.nrows, self.ncols),
                                      CPUs=self.cfg.CPUs)\
            .compute_geolayer()

        return lons, lats

    def calc_solar_irradiance_CWL_FWHM_per_band(self) -> np.array:
        from ...io.reader import Solar_Irradiance_reader

        self.logger.debug('Calculating solar irradiance...')

        sol_irr = Solar_Irradiance_reader(path_solar_irr_model=self.cfg.path_solar_irr, wvl_min_nm=350, wvl_max_nm=2500)

        irr_bands = []
        for band in self.srf.bands:
            WVL_band = self.srf.srfs_wvl if self.srf.wvl_unit == 'nanometers' else self.srf.srfs_wvl * 1000
            RSP_band = self.srf.srfs_norm01[band]
            sol_irr_at_WVL = np.interp(WVL_band, sol_irr[:, 0], sol_irr[:, 1], left=0, right=0)

            irr_bands.append(np.round(np.sum(sol_irr_at_WVL * RSP_band) / np.sum(RSP_band), 2))

        return np.array(irr_bands)


class EnMAP_Metadata_L1B_SensorGeo(object):
    """EnMAP Metadata class for the metadata of the complete EnMAP L1B product in sensor geometry incl. VNIR and SWIR.

    Attributes:
        - logger(logging.Logger):  None or logging instance
        - observation_datetime(datetime.datetime):  datetime of observation time
        - geom_view_zenith: viewing zenith angle
        - geom_view_azimuth: viewing azimuth angle
        - geom_sun_zenith: sun zenith angle
        - geom_sun_azimuth: sun azimuth angle
        - mu_sun: needed by SICOR for TOARad > TOARef conversion
        - vnir(EnMAP_Metadata_VNIR_SensorGeo)
        - swir(EnMAP_Metadata_SWIR_SensorGeo)
        - detector_attrNames: attribute names of the detector objects
    """

    def __init__(self, path_metaxml, config: EnPTConfig, logger=None):
        """Get a metadata object instance for the given EnMAP L1B product in sensor geometry.

        :param path_metaxml:    file path of the EnMAP L1B metadata XML file
        :para, config:          EnPT configuration object
        :param logger:          instance of logging.logger or subclassed
        """
        self.cfg = config
        self.logger = logger or logging.getLogger()
        self.path_xml = path_metaxml
        self.rootdir = os.path.dirname(path_metaxml)

        # defaults - Common
        self.proc_level: Optional[str] = None   # Dataset processing level
        self.observation_datetime: Optional[datetime] = None  # Date and Time of image observation
        self.geom_view_zenith: Optional[float] = None  # viewing zenith angle
        self.geom_view_azimuth: Optional[float] = None  # viewing azimuth angle
        self.geom_sun_zenith: Optional[float] = None  # sun zenith angle
        self.geom_sun_azimuth: Optional[float] = None   # sun azimuth angle
        self.mu_sun: Optional[float] = None   # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist: Optional[float] = None  # earth-sun distance
        self.aot: Optional[float] = None  # scene aerosol optical thickness
        self.water_vapour: Optional[float] = None  # scene water vapour [cm]
        self.vnir: Optional[EnMAP_Metadata_L1B_Detector_SensorGeo] = None  # metadata of VNIR only
        self.swir: Optional[EnMAP_Metadata_L1B_Detector_SensorGeo] = None  # metadata of SWIR only
        self.detector_attrNames: list = ['vnir', 'swir']  # attribute names of the detector objects
        self.filename_metaxml: Optional[str] = None  # filename of XML metadata file

        self._scene_basename: Optional[str] = None  # basename of the EnMAP image

    @property
    def scene_basename(self):
        if self.vnir:
            self._scene_basename = self.vnir.scene_basename

        return self._scene_basename

    def read_common_meta(self, path_xml):
        """Read the common metadata, principally stored in General Info.

        - the acquisition time
        - the geometrical observation and illumination

        :param path_xml: path to the main xml file
        :return: None
        """
        # load the metadata xml file
        xml = ElementTree.parse(path_xml).getroot()

        self.filename_metaxml = os.path.basename(path_xml)

        if not self.cfg.is_dummy_dataformat:
            # read processing level
            self.proc_level = xml.find("base/level").text
            if self.proc_level != 'L1B':
                raise RuntimeError(self.proc_level, "Unexpected input data processing level. Expected 'L1B'.")

            # read the acquisition time
            self.observation_datetime = \
                datetime.strptime(xml.find("base/temporalCoverage/startTime").text, '%Y-%m-%dT%H:%M:%S.%fZ')

            # get the distance earth sun from the acquisition date
            self.earthSunDist = self.get_earth_sun_distance(self.observation_datetime)

            # read Geometry (observation/illumination) angle
            # NOTE: EnMAP metadata provide also the angles for the image corners
            #       -> would allow even more precise computation (e.g., specific/sunElevationAngle/upper_left)
            # NOTE: alongOffNadirAngle is always near 0 and therefore ignored here (not relevant for AC)
            # FIXME VZA may be negative in DLR L1B data -> correct to always use the absolute value for SICOR?
            self.geom_view_zenith = np.abs(float(xml.find("specific/acrossOffNadirAngle/center").text))
            # FIXME correct to directly use sceneAzimuthAngle (14.24 (DLR) vs. 101.1 (AlpineTest)
            self.geom_view_azimuth = float(xml.find("specific/sceneAzimuthAngle/center").text)
            self.geom_sun_zenith = 90 - float(xml.find("specific/sunElevationAngle/center").text)
            self.geom_sun_azimuth = float(xml.find("specific/sunAzimuthAngle/center").text)
            self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))
            self.aot = float(xml.find("specific/qualityFlag/sceneAOT").text) / 1000  # scale factor is 1000
            self.water_vapour = float(xml.find("specific/qualityFlag/sceneWV").text) / 1000  # scale factor is 1000
        else:
            # read the acquisition time
            self.observation_datetime = \
                datetime.strptime(xml.findall("GeneralInfo/ProductInfo/ProductStartTime")[0].text,
                                  '%Y-%m-%dT%H:%M:%S.%fZ')

            # get the distance earth sun from the acquisition date
            self.earthSunDist = self.get_earth_sun_distance(self.observation_datetime)

            # read Geometry (observation/illumination) angle
            self.geom_view_zenith = float(xml.findall("GeneralInfo/Geometry/Observation/ZenithAngle")[0].text)
            self.geom_view_azimuth = float(xml.findall("GeneralInfo/Geometry/Observation/AzimuthAngle")[0].text)
            self.geom_sun_zenith = float(xml.findall("GeneralInfo/Geometry/Illumination/ZenithAngle")[0].text)
            self.geom_sun_azimuth = float(xml.findall("GeneralInfo/Geometry/Illumination/AzimuthAngle")[0].text)
            self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))

    def get_earth_sun_distance(self, acqDate: datetime):
        """Get earth sun distance (requires file of pre calculated earth sun distance per day).

        :param acqDate:
        """
        if not os.path.exists(self.cfg.path_earthSunDist):
            self.logger.warning("\n\t WARNING: Earth Sun Distance is assumed to be "
                                "1.0 because no database can be found at %s.""" % self.cfg.path_earthSunDist)
            return 1.0
        if not acqDate:
            self.logger.warning("\n\t WARNING: Earth Sun Distance is assumed to be 1.0 because "
                                "acquisition date could not be read from metadata.")
            return 1.0

        with open(self.cfg.path_earthSunDist, "r") as EA_dist_f:
            EA_dist_dict = {}
            for line in EA_dist_f:
                date, EA = [item.strip() for item in line.split(",")]
                EA_dist_dict[date] = EA

        return float(EA_dist_dict[acqDate.strftime('%Y-%m-%d')])

    def read_metadata(self):
        """Read the metadata of the entire EnMAP L1B product in sensor geometry."""
        # first read common metadata
        self.read_common_meta(self.path_xml)

        # define and read the VNIR metadata
        self.vnir = EnMAP_Metadata_L1B_Detector_SensorGeo('VNIR', config=self.cfg, logger=self.logger)
        self.vnir.read_metadata(self.path_xml)

        # define and read the SWIR metadata
        self.swir = EnMAP_Metadata_L1B_Detector_SensorGeo('SWIR', config=self.cfg, logger=self.logger)
        self.swir.read_metadata(self.path_xml)

    def to_XML(self) -> str:
        """Generate an XML metadata string from the L1B metadata."""
        from . import L1B_product_props, L1B_product_props_DLR
        xml = ElementTree.parse(self.path_xml).getroot()

        if not self.cfg.is_dummy_dataformat:
            for detName, detMeta in zip(['VNIR', 'SWIR'], [self.vnir, self.swir]):
                lbl = L1B_product_props_DLR['xml_detector_label'][detName]
                xml.find("product/image/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
                xml.find("product/image/%s/dimension/columns" % lbl).text = str(detMeta.ncols)
                xml.find("product/quicklook/%s/dimension/rows" % lbl).text = str(detMeta.nrows)
                xml.find("product/quicklook/%s/dimension/columns" % lbl).text = str(detMeta.ncols)

        else:
            for detName, detMeta in zip(['VNIR', 'SWIR'], [self.vnir, self.swir]):
                lbl = L1B_product_props['xml_detector_label'][detName]
                xml.find("ProductComponent/%s/Data/Size/NRows" % lbl).text = str(detMeta.nrows)
                xml.find("ProductComponent/%s/Data/Type/UnitCode" % lbl).text = detMeta.unitcode
                xml.find("ProductComponent/%s/Data/Type/Unit" % lbl).text = detMeta.unit

        xml_string = ElementTree.tostring(xml, encoding='unicode', pretty_print=True)

        return xml_string
