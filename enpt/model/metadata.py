# -*- coding: utf-8 -*-

# EnPT, EnMAP Processing Tool - A Python package for pre-processing of EnMAP Level-1B data
#
# Copyright (C) 2019  Karl Segl (GFZ Potsdam, segl@gfz-potsdam.de), Daniel Scheffler
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

"""EnPT metadata modules. All object and functions regarding EnMAP metadata are implemented here."""

from datetime import datetime
from lxml import etree as ElementTree
import logging
import os
import fnmatch
from typing import Union, List, Tuple  # noqa: F401
from collections import OrderedDict
import numpy as np
from py_tools_ds.geo.vector.topology import Polygon, get_footprint_polygon  # noqa: F401  # flake8 issue
from geoarray import GeoArray

from ..options.config import EnPTConfig, enmap_xres
from .srf import SRF
from ..processors.spatial_transform import RPC_3D_Geolayer_Generator

__author__ = ['Daniel Scheffler', 'Stéphane Guillaso', 'André Hollstein']


# Define L1B_product_props
L1B_product_props = dict(
    xml_detector_label=dict(
        VNIR='VNIRDetector',
        SWIR='SWIRDetector'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)


L1B_product_props_DLR = dict(
    xml_detector_label=dict(
        VNIR='vnir',
        SWIR='swir'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)


# Define L1B_product_props
L2A_product_props_DLR = dict(
    xml_detector_label=dict(
        VNIR='vnir',
        SWIR='swir'
    ),
    fn_detector_suffix=dict(
        VNIR='D1',
        SWIR='D2'
    )
)


#########################################################
# EnPT metadata objects for EnMAP data in sensor geometry
#########################################################


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
        self.cfg = config
        self.detector_name = detector_name  # type: str
        if not self.cfg.is_dummy_dataformat:
            self.detector_label = L1B_product_props_DLR['xml_detector_label'][detector_name]
        else:
            self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.logger = logger or logging.getLogger()

        # These lines are used to load path information
        self.data_filename = None  # type: str # detector data filename
        self.scene_basename = None  # type: str # basename of the EnMAP image
        self.dead_pixel_filename = None  # type: str # filename of the dead pixel file
        self.quicklook_filename = None  # type: str # filename of the quicklook file
        # FIXME cloud mask of BOTH detectors
        self.cloud_mask_filename = None  # type: str # filename of the cloud mask file

        self.wvl_center = None  # type: np.ndarray  # Center wavelengths for each EnMAP band
        self.fwhm = None  # type: np.ndarray  # Full width half maximmum for each EnMAP band
        self.srf = None  # type: SRF  # SRF object holding the spectral response functions for each EnMAP band
        self.solar_irrad = None  # type: np.array  # solar irradiance in [W/m2/nm] for each band
        self.nwvl = None  # type: int  # Number of wave bands
        self.nrows = None  # type: int  # number of rows
        self.ncols = None  # type: int  # number of columns
        self.smile_coef = None  # type: np.ndarray  # smile coefficients needed for smile computation
        self.nsmile_coef = None  # type: int  # number of smile coefficients
        self.smile = None  # type: np.ndarray  # smile for each EnMAP image column
        self.gains = None  # type: np.ndarray  # band-wise gains for computing radiance from DNs
        self.offsets = None  # type: np.ndarray  # band-wise offsets for computing radiance from DNs
        self.l_min = None  # type: np.ndarray  # band-wise l-min for computing radiance from DNs
        self.l_max = None  # type: np.ndarray  # band-wise l-max for computing radiance from DNs
        self.lat_UL_UR_LL_LR = None  # type:  List[float, float, float, float]  # latitude coords for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR = None  # type:  List[float, float, float, float]  # longitude coords for UL, UR, LL, LR
        self.rpc_coeffs = OrderedDict()  # type: OrderedDict  # RPC coefficients for geolayer computation
        self.ll_mapPoly = None  # type: Polygon  # footprint polygon in longitude/latitude map coordinates
        self.lats = None  # type: np.ndarray  # 2D array of latitude coordinates according to given lon/lat sampling
        self.lons = None  # type: np.ndarray  # 2D array of longitude coordinates according to given lon/lat sampling
        self.unit = ''  # type: str  # radiometric unit of pixel values
        self.unitcode = ''  # type: str  # code of radiometric unit
        self.preview_bands = None
        self.snr = None  # type: np.ndarray  # Signal to noise ratio as computed from radiance data

    def read_metadata(self, path_xml):
        """
        Read the metadata of a specific EnMAP detector in sensor geometry

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

            self.data_filename = xml.find("product/image/%s/name" % lbl).text
            self.scene_basename = self.data_filename.split('-SPECTRAL_IMAGE')[0]
            self.dead_pixel_filename = get_filename('*QL_PIXELMASK_%s' % self.detector_name)
            self.quicklook_filename = xml.find("product/quicklook/%s/name" % lbl).text
            # FIXME multiple cloud masks provided. QL_QUALITY_CLASSES as combined product?
            #   - QL_QUALITY_CLOUD
            #   - QL_QUALITY_CIRRUS
            #   - QL_QUALITY_SNOW
            #   - QL_QUALITY_CLOUDSHADOW
            #   - QL_QUALITY_HAZE
            self.cloud_mask_filename = get_filename('*-QL_QUALITY_CLOUD')
            self.logger.warning('DLR test data provide multiple cloud masks. Added only *%s!'
                                % self.cloud_mask_filename.split(self.scene_basename)[1])

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
            self.smile_coef = np.zeros((self.nwvl, self.nsmile_coef), dtype=np.float)

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
            self.preview_bands = np.array([np.argmin(np.abs(self.wvl_center - wvl))
                                           for wvl in [wvl_red, wvl_green, wvl_blue]])

            # read RPC coefficients
            for bID in xml.findall("product/navigation/RPC/bandID")[subset]:
                bN = 'band_%d' % np.int64(bID.attrib['number'])

                keys2combine = ('row_num', 'row_den', 'col_num', 'col_den')

                tmp = OrderedDict([(ele.tag.lower(), np.float(ele.text)) for ele in bID.findall('./')])
                self.rpc_coeffs[bN] = {k: v for k, v in tmp.items() if not k.startswith(keys2combine)}

                for n in keys2combine:
                    self.rpc_coeffs[bN]['%s_coeffs' % n.lower()] = \
                        np.array([v for k, v in tmp.items() if k.startswith(n)])

            # compute metadata derived from read data
            self.smile = self.calc_smile()
            self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
            self.solar_irrad = self.calc_solar_irradiance_CWL_FWHM_per_band()
            self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                              self.lat_UL_UR_LL_LR)), fix_invalid=True)

        else:
            lbl = self.detector_label
            self.logger.info("Reading metadata for %s detector..." % self.detector_name)

            # read data filenames
            self.data_filename = xml.findall("ProductComponent/%s/Data/Filename" % lbl)[0].text
            self.scene_basename = os.path.splitext(self.data_filename)[0]
            self.dead_pixel_filename = xml.findall("ProductComponent/%s/Sensor/DeadPixel/Filename" % lbl)[0].text
            self.quicklook_filename = xml.findall("ProductComponent/%s/Preview/Filename" % lbl)[0].text
            self.cloud_mask_filename = xml.findall("ProductComponent/%s/Data/CloudMaskMap/Filename" % lbl)[0].text

            # read preview bands
            self.preview_bands = np.zeros(3, dtype=np.int)
            self.preview_bands[0] = np.int(xml.findall("ProductComponent/%s/Preview/Bands/Red" % lbl)[0].text)
            self.preview_bands[1] = np.int(xml.findall("ProductComponent/%s/Preview/Bands/Green" % lbl)[0].text)
            self.preview_bands[2] = np.int(xml.findall("ProductComponent/%s/Preview/Bands/Blue" % lbl)[0].text)

            # read some basic information concerning the detector
            self.nrows = np.int(xml.findall("ProductComponent/%s/Data/Size/NRows" % lbl)[0].text)
            self.ncols = np.int(xml.findall("ProductComponent/%s/Data/Size/NCols" % lbl)[0].text)
            self.unitcode = xml.findall("ProductComponent/%s/Data/Type/UnitCode" % lbl)[0].text
            self.unit = xml.findall("ProductComponent/%s/Data/Type/Unit" % lbl)[0].text

            # Read image coordinates
            scene_corner_coordinates = xml.findall("ProductComponent/%s/Data/SceneInformation/"
                                                   "SceneCornerCoordinates" % lbl)
            self.lat_UL_UR_LL_LR = [
                np.float(scene_corner_coordinates[0].findall("Latitude")[0].text),
                np.float(scene_corner_coordinates[1].findall("Latitude")[0].text),
                np.float(scene_corner_coordinates[2].findall("Latitude")[0].text),
                np.float(scene_corner_coordinates[3].findall("Latitude")[0].text)
            ]
            self.lon_UL_UR_LL_LR = [
                np.float(scene_corner_coordinates[0].findall("Longitude")[0].text),
                np.float(scene_corner_coordinates[1].findall("Longitude")[0].text),
                np.float(scene_corner_coordinates[2].findall("Longitude")[0].text),
                np.float(scene_corner_coordinates[3].findall("Longitude")[0].text)
            ]

            # read the band related information: wavelength, fwhm
            self.nwvl = np.int(xml.findall("ProductComponent/%s/Data/BandInformationList/NumberOfBands" % lbl)[0].text)
            self.nsmile_coef = np.int(xml.findall(
                "ProductComponent/%s/Data/BandInformationList/SmileInformation/NumberOfCoefficients" % lbl)[0].text)
            self.fwhm = np.zeros(self.nwvl, dtype=np.float)
            self.wvl_center = np.zeros(self.nwvl, dtype=np.float)
            self.smile_coef = np.zeros((self.nwvl, self.nsmile_coef), dtype=np.float)
            self.l_min = np.zeros(self.nwvl, dtype=np.float)
            self.l_max = np.zeros(self.nwvl, dtype=np.float)
            band_informations = xml.findall("ProductComponent/%s/Data/BandInformationList/BandInformation" % lbl)
            for bi in band_informations:
                k = np.int64(bi.attrib['Id']) - 1
                self.wvl_center[k] = np.float(bi.findall("CenterWavelength")[0].text)
                self.fwhm[k] = np.float(bi.findall("FullWidthHalfMaximum")[0].text)
                self.l_min[k] = np.float(bi.findall("L_min")[0].text)
                self.l_max[k] = np.float(bi.findall("L_max")[0].text)
                scl = bi.findall("Smile/Coefficient")
                for sc in scl:
                    self.smile_coef[k, np.int64(sc.attrib['exponent'])] = np.float(sc.text)
            self.smile = self.calc_smile()
            self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
            self.solar_irrad = self.calc_solar_irradiance_CWL_FWHM_per_band()
            self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR, self.lat_UL_UR_LL_LR)),
                                                    fix_invalid=True)
            self.lats = self.interpolate_corners(*self.lat_UL_UR_LL_LR, self.ncols, self.nrows)
            self.lons = self.interpolate_corners(*self.lon_UL_UR_LL_LR, self.ncols, self.nrows)

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

        lons, lats = RPC_3D_Geolayer_Generator(rpc_coeffs_per_band=self.rpc_coeffs,
                                               dem=self.cfg.path_dem,
                                               enmapIm_cornerCoords=tuple(zip(self.lon_UL_UR_LL_LR,
                                                                              self.lat_UL_UR_LL_LR)),
                                               enmapIm_dims_sensorgeo=(self.nrows, self.ncols),
                                               CPUs=self.cfg.CPUs)\
            .compute_geolayer()

        return lons, lats

    def calc_solar_irradiance_CWL_FWHM_per_band(self) -> np.array:
        from ..io.reader import Solar_Irradiance_reader

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
        self.proc_level = None  # type: str  # Dataset processing level
        self.observation_datetime = None  # type: datetime  # Date and Time of image observation
        self.geom_view_zenith = None  # type: float  # viewing zenith angle
        self.geom_view_azimuth = None  # type: float  # viewing azimuth angle
        self.geom_sun_zenith = None  # type: float  # sun zenith angle
        self.geom_sun_azimuth = None  # type: float  # sun azimuth angle
        self.mu_sun = None  # type: float  # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist = None  # type: float  # earth-sun distance
        self.aot = None  # type: float  # scene aerosol optical thickness
        self.water_vapour = None  # type: float  # scene water vapour [cm]
        self.vnir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of VNIR only
        self.swir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of SWIR only
        self.detector_attrNames = ['vnir', 'swir']  # type: list # attribute names of the detector objects
        self.metaxml_filename = None  # type: str # filename of XML metadata file

        self._scene_basename = None  # type: str # basename of the EnMAP image

    @property
    def scene_basename(self):
        if self.vnir:
            self._scene_basename = self.vnir.scene_basename

        return self._scene_basename

    # Read common metadata method
    def read_common_meta(self, path_xml):
        """Read the common metadata, principally stored in General Info
        - the acquisition time
        - the geometrical observation and illumination
        :param path_xml: path to the main xml file
        :return: None
        """

        # load the metadata xml file
        xml = ElementTree.parse(path_xml).getroot()

        self.metaxml_filename = os.path.basename(path_xml)

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
            self.geom_view_zenith = np.abs(np.float(xml.find("specific/acrossOffNadirAngle/center").text))
            # FIXME correct to directly use sceneAzimuthAngle (14.24 (DLR) vs. 101.1 (AlpineTest)
            self.geom_view_azimuth = np.float(xml.find("specific/sceneAzimuthAngle/center").text)
            self.geom_sun_zenith = 90 - np.float(xml.find("specific/sunElevationAngle/center").text)
            self.geom_sun_azimuth = np.float(xml.find("specific/sunAzimuthAngle/center").text)
            self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))
            self.aot = np.float(xml.find("specific/qualityFlag/sceneAOT").text) / 1000  # scale factor is 1000
            self.water_vapour = np.float(xml.find("specific/qualityFlag/sceneWV").text) / 1000  # scale factor is 1000
        else:
            # read the acquisition time
            self.observation_datetime = \
                datetime.strptime(xml.findall("GeneralInfo/ProductInfo/ProductStartTime")[0].text,
                                  '%Y-%m-%dT%H:%M:%S.%fZ')

            # get the distance earth sun from the acquisition date
            self.earthSunDist = self.get_earth_sun_distance(self.observation_datetime)

            # read Geometry (observation/illumination) angle
            self.geom_view_zenith = np.float(xml.findall("GeneralInfo/Geometry/Observation/ZenithAngle")[0].text)
            self.geom_view_azimuth = np.float(xml.findall("GeneralInfo/Geometry/Observation/AzimuthAngle")[0].text)
            self.geom_sun_zenith = np.float(xml.findall("GeneralInfo/Geometry/Illumination/ZenithAngle")[0].text)
            self.geom_sun_azimuth = np.float(xml.findall("GeneralInfo/Geometry/Illumination/AzimuthAngle")[0].text)
            self.mu_sun = np.cos(np.deg2rad(self.geom_sun_zenith))

    def get_earth_sun_distance(self, acqDate: datetime):
        """Get earth sun distance (requires file of pre calculated earth sun distance per day)

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
        """
        Read the metadata of the entire EnMAP L1B product in sensor geometry

        :return: None
        """

        # first read common metadata
        self.read_common_meta(self.path_xml)

        # define and read the VNIR metadata
        self.vnir = EnMAP_Metadata_L1B_Detector_SensorGeo('VNIR', config=self.cfg, logger=self.logger)
        self.vnir.read_metadata(self.path_xml)

        # define and read the SWIR metadata
        self.swir = EnMAP_Metadata_L1B_Detector_SensorGeo('SWIR', config=self.cfg, logger=self.logger)
        self.swir.read_metadata(self.path_xml)

    def to_XML(self) -> str:
        """
        Generate an XML metadata string from the L1B metadata.
        """
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


class EnMAP_Metadata_L2A_MapGeo(object):
    def __init__(self,
                 config: EnPTConfig,
                 meta_l1b: EnMAP_Metadata_L1B_SensorGeo,
                 wvls_l2a: Union[List, np.ndarray],
                 dims_mapgeo: Tuple[int, int, int],
                 logger=None):
        """EnMAP Metadata class for the metadata of the complete EnMAP L2A product in map geometry incl. VNIR and SWIR.

        :param config:              EnPT configuration object
        :param meta_l1b:            metadata object of the L1B dataset in sensor geometry
        :param wvls_l2a:            list of center wavelengths included in the L2A product
        :param dims_mapgeo:         dimensions of the EnMAP raster data in map geometry, e.g., (1024, 1000, 218)
        :param logger:              instance of logging.logger or subclassed
        """
        self.cfg = config
        self._meta_l1b = meta_l1b
        self.logger = logger or logging.getLogger()

        # defaults
        self.band_means = None  # type: np.ndarray # band-wise means in unscaled values (percent in case of reflectance)
        self.band_stds = None  # type: np.ndarray # band-wise standard deviations in unscaled values
        self.fileinfos = []  # type: list # file informations for each file beloning to the EnMAP L2A product

        self.proc_level = 'L2A'
        self.observation_datetime = meta_l1b.observation_datetime  # type: datetime  # Date and Time of observation
        # FIXME VZA may be negative in DLR data
        self.geom_view_zenith = meta_l1b.geom_view_zenith  # type: float  # viewing zenith angle
        self.geom_view_azimuth = meta_l1b.geom_view_azimuth  # type: float  # viewing azimuth angle
        self.geom_sun_zenith = meta_l1b.geom_sun_zenith  # type: float  # sun zenith angle
        self.geom_sun_azimuth = meta_l1b.geom_sun_azimuth  # type: float  # sun azimuth angle
        self.mu_sun = meta_l1b.mu_sun  # type: float  # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist = meta_l1b.earthSunDist  # type: float  # earth-sun distance

        # generate file names for L2A output
        if not self.cfg.is_dummy_dataformat:
            self.scene_basename = meta_l1b.vnir.data_filename.split('-SPECTRAL_IMAGE')[0].replace('L1B-', 'L2A-')
        else:
            self.scene_basename = os.path.splitext(meta_l1b.vnir.data_filename)[0]
        self.data_filename = meta_l1b.vnir.data_filename.replace('L1B-', 'L2A-').replace('_VNIR', '')
        self.dead_pixel_filename_vnir = meta_l1b.vnir.dead_pixel_filename.replace('L1B-', 'L2A-')
        self.dead_pixel_filename_swir = meta_l1b.swir.dead_pixel_filename.replace('L1B-', 'L2A-')
        self.quicklook_filename_vnir = meta_l1b.vnir.quicklook_filename.replace('L1B-', 'L2A-')
        self.quicklook_filename_swir = meta_l1b.swir.quicklook_filename.replace('L1B-', 'L2A-')
        self.cloud_mask_filename = meta_l1b.vnir.cloud_mask_filename.replace('L1B-', 'L2A-')
        self.metaxml_filename = meta_l1b.metaxml_filename.replace('L1B-', 'L2A-')

        # fuse band-wise metadata (sort all band-wise metadata by wavelengths but band number keeps as it is)
        # get band index order
        wvls_sorted = np.array(sorted(np.hstack([self._meta_l1b.vnir.wvl_center,
                                                 self._meta_l1b.swir.wvl_center])))
        bandidx_order = np.array([np.argmin(np.abs(wvls_sorted - cwl)) for cwl in wvls_l2a])

        self.wvl_center = np.hstack([meta_l1b.vnir.wvl_center, meta_l1b.swir.wvl_center])[bandidx_order]
        self.fwhm = np.hstack([meta_l1b.vnir.fwhm, meta_l1b.swir.fwhm])[bandidx_order]
        self.gains = np.full((dims_mapgeo[2],), 100)  # implies reflectance scaled between 0 and 10000
        self.offsets = np.zeros((dims_mapgeo[2],))
        self.srf = SRF.from_cwl_fwhm(self.wvl_center, self.fwhm)
        self.solar_irrad = np.hstack([meta_l1b.vnir.solar_irrad, meta_l1b.swir.solar_irrad])[bandidx_order]

        if not meta_l1b.vnir.nsmile_coef == meta_l1b.swir.nsmile_coef:
            raise ValueError('Expected equal number of smile coefficients for VNIR and SWIR. Received %d/%s.'
                             % (meta_l1b.vnir.nsmile_coef, meta_l1b.swir.nsmile_coef))

        self.nsmile_coef = meta_l1b.vnir.nsmile_coef
        self.smile_coef = np.vstack([meta_l1b.vnir.smile_coef, meta_l1b.swir.smile_coef])[bandidx_order, :]
        self.smile = np.hstack([meta_l1b.vnir.smile, meta_l1b.swir.smile])[:, bandidx_order]

        if not self.cfg.is_dummy_dataformat:
            self.rpc_coeffs = OrderedDict(zip(
                ['band_%d' % (i + 1) for i in range(dims_mapgeo[2])],
                [meta_l1b.vnir.rpc_coeffs['band_%d' % (i + 1)] if 'band_%d' % (i + 1) in meta_l1b.vnir.rpc_coeffs else
                 meta_l1b.swir.rpc_coeffs['band_%d' % (i + 1)] for i in bandidx_order]))
        else:
            self.rpc_coeffs = OrderedDict()

        self.nrows = dims_mapgeo[0]
        self.ncols = dims_mapgeo[1]
        self.nwvl = dims_mapgeo[2]
        common_UL_UR_LL_LR = self.get_common_UL_UR_LL_LR()
        self.lon_UL_UR_LL_LR = [lon for lon, lat in common_UL_UR_LL_LR]
        self.lat_UL_UR_LL_LR = [lat for lon, lat in common_UL_UR_LL_LR]
        self.ll_mapPoly = get_footprint_polygon(tuple(zip(self.lon_UL_UR_LL_LR,
                                                          self.lat_UL_UR_LL_LR)), fix_invalid=True)

        if meta_l1b.vnir.unit != meta_l1b.swir.unit or meta_l1b.vnir.unitcode != meta_l1b.swir.unitcode:
            raise RuntimeError('L2A data should have the same radiometric unit for VNIR and SWIR. '
                               'Received %s in %s for VNIR and %s in %s for SWIR.'
                               % (meta_l1b.vnir.unitcode, meta_l1b.vnir.unit,
                                  meta_l1b.swir.unitcode, meta_l1b.swir.unit))

        self.unit = meta_l1b.vnir.unit
        self.unitcode = meta_l1b.vnir.unitcode
        self.preview_bands_vnir = meta_l1b.vnir.preview_bands
        self.preview_bands_swir = meta_l1b.swir.preview_bands

        self.snr = None
        if meta_l1b.vnir.snr is not None:
            assert meta_l1b.swir.snr is not None
            self.snr = np.dstack([meta_l1b.vnir.snr, meta_l1b.swir.snr])[:, :, bandidx_order]

    def get_common_UL_UR_LL_LR(self):
        vnir_ulx, vnir_urx, vnir_llx, vnir_lrx = self._meta_l1b.vnir.lon_UL_UR_LL_LR
        vnir_uly, vnir_ury, vnir_lly, vnir_lry = self._meta_l1b.vnir.lat_UL_UR_LL_LR
        swir_ulx, swir_urx, swir_llx, swir_lrx = self._meta_l1b.swir.lon_UL_UR_LL_LR
        swir_uly, swir_ury, swir_lly, swir_lry = self._meta_l1b.swir.lat_UL_UR_LL_LR

        # use OUTER coordinates
        return ((min([vnir_ulx, swir_ulx]), max([vnir_uly, swir_uly])),
                (max([vnir_urx, swir_urx]), max([vnir_ury, swir_ury])),
                (min([vnir_llx, swir_llx]), min([vnir_lly, swir_lly])),
                (max([vnir_lrx, swir_lrx]), min([vnir_lry, swir_lry])))

    def add_band_statistics(self, datastack_vnir_swir: Union[np.ndarray, GeoArray]):
        R, C, B = datastack_vnir_swir.shape
        # NOTE:  DEVIDE by gains to reflectance in percent
        self.band_means = np.mean(datastack_vnir_swir.reshape(1, R * C, B), axis=1) / self.gains
        self.band_stds = np.mean(datastack_vnir_swir.reshape(1, R * C, B), axis=1) / self.gains

    def add_product_fileinformation(self, filepaths: List[str], sizes: List[int] = None, versions: List[str] = None):
        self.fileinfos = []

        for i, fp in enumerate(filepaths):
            ismeta = fp.endswith('METADATA.XML') or fp.endswith('_header.xml')  # FIXME
            if not os.path.exists(fp):
                if ismeta:
                    pass  # does not yet exist
                else:
                    raise FileNotFoundError(fp)

            ext = os.path.splitext(fp)[1]
            fileinfo_dict = dict(
                name=os.path.basename(fp),
                size=sizes[i] if sizes else int(os.path.getsize(fp) / 1024) if not ismeta else '',
                version=versions[i] if versions else '',
                format='binary' if ext in ['.GEOTIFF',
                                           '.TIF',
                                           '.TIFF',
                                           '.GTIFF',
                                           '.BSQ',
                                           '.BIL',
                                           '.BIP',
                                           '.JPEG2000'] else 'xml' if ext == '.XML' else 'NA'
            )

            self.fileinfos.append(fileinfo_dict)

    def to_XML(self) -> str:
        """
        Generate an XML metadata string from the L2A metadata.
        """
        # use an XML parser that creates properly indented XML files even if new SubElements have been added
        parser = ElementTree.XMLParser(remove_blank_text=True)

        # parse (use L1B metadata as template)
        xml = ElementTree.parse(self._meta_l1b.path_xml, parser).getroot()

        if self.cfg.is_dummy_dataformat:
            self.logger.warning('No XML metadata conversion implemented for datasets different to the DLR format.'
                                'Metadata XML file will be empty.')
            return ''

        self.logger.warning('Currently, the L2A metadata XML file does not contain all relevant keys and contains '
                            'not updated values!')  # FIXME

        ############
        # metadata #
        ############

        xml.find("metadata/schema/processingLevel").text = self.proc_level
        xml.find("metadata/name").text = self.metaxml_filename
        # xml.find("metadata/comment").text = 'EnMAP Level 0 Product of datatake 987'  # FIXME hardcoded

        ##############
        # processing #
        ##############

        xml.find("processing/terrainCorrection").text = 'Yes'  # FIXME hardcoded {Yes, No}
        xml.find("processing/ozoneValue").text = 'NA'  # FIXME {[200-500], NA}
        xml.find("processing/season").text = 'NA'  # FIXME {summer, winter, NA}
        xml.find("processing/productFormat").text = 'GeoTIFF+Metadata'  # FIXME hardcoded
        # {BSQ+Metadata, BIL+Metadata, BIP+Metadata, JPEG2000+Metadata, GeoTiff+Metadata}
        xml.find("processing/mapProjection").text = 'UTM_Zone_of_Scene_Center'  # FIXME hardcoded
        # {UTM_Zone_of_Scene_Center, UTM_Zone_of_Scene_Center(-1), UTM_Zone_of_Scene_Center(+1),
        #  UTM_Zone_of_Datatake_Center, Geographic, European_Projection_LAEA, NA}
        xml.find("processing/DEMDBVersion").text = 'SRTM-C_v4'  # FIXME hardcoded
        # {SRTM-C-X_vv.rr, best-of-DEM_vv.rr, DEM-derivedfrom-Tandem-X_vv.rr, ASTER-GDEM_vv.rr, NA}
        xml.find("processing/correctionType").text = 'NA'  # FIXME hardcoded {Combined, Land_Mode, Water_Mode, NA}
        xml.find("processing/cirrusHazeRemoval").text = 'NA'  # FIXME hardcoded {Yes, No}
        xml.find("processing/bandInterpolation").text = 'NA'  # FIXME hardcoded {Yes, No}
        xml.find("processing/waterType").text = 'NA'  # FIXME hardcoded {Clear, Turbid, Highly_Turbid, NA}

        ########
        # base #
        ########

        # TODO update corner coordinates? DLR just uses the same coordinates like in L1B
        # xml.find("base/spatialCoverage" % lbl).text =
        xml.find("base/format").text = 'ENMAP_%s' % self.proc_level
        xml.find("base/level").text = self.proc_level
        xml.find("base/size").text = 'NA'  # FIXME Size of product. Attribute unit {byte, Kbyte, Mbyte, Gbyte}

        ############
        # specific #
        ############

        xml.find("specific/code").text = self.proc_level
        bi = "specific/bandCharacterisation/bandID/"
        for ele, gain in zip(xml.findall(bi + "GainOfBand"), self.gains):
            ele.text = str(gain)
        for ele, offset in zip(xml.findall(bi + "OffsetOfBand"), self.offsets):
            ele.text = str(offset)

        ###########
        # product #
        ###########

        if not self.fileinfos:
            raise ValueError('Product file informations must be added before writing metadata. '
                             'Call add_product_fileinformation() before!')

        for detName, detMetaL1B in zip(['VNIR', 'SWIR'], [self._meta_l1b.vnir, self._meta_l1b.swir]):
            lbl = L2A_product_props_DLR['xml_detector_label'][detName]
            # FIXME DLR uses L0 filenames for VNIR/SWIR separately?!
            xml.find("product/image/%s/name" % lbl).text = detMetaL1B.data_filename
            # FIXME this is the size of the VNIR/SWIR stack
            size = [F['size'] for F in self.fileinfos if os.path.splitext(F['name'])[0].endswith('-SPECTRAL_IMAGE')][0]
            xml.find("product/image/%s/size" % lbl).text = str(size)
            # FIXME DLR data dimensions equal neither L2A data nor L1B data
            xml.find("product/image/%s/channels" % lbl).text = str(detMetaL1B.nwvl)
            xml.find("product/image/%s/dimension/rows" % lbl).text = str(self.nrows)
            xml.find("product/image/%s/dimension/columns" % lbl).text = str(self.ncols)
            # xml.find("product/image/%s/dimension/dimensionGeographic/longitude" % lbl).text = 'NA'  # TODO
            # xml.find("product/image/%s/dimension/dimensionGeographic/latitude" % lbl).text = 'NA'

            fN_quicklook = self.quicklook_filename_vnir if detName == 'VNIR' else self.quicklook_filename_swir
            size_quicklook = [F['size'] for F in self.fileinfos
                              if os.path.splitext(F['name'])[0].endswith('-QL_%s' % detName)][0]
            xml.find("product/quicklook/%s/name" % lbl).text = fN_quicklook
            xml.find("product/quicklook/%s/size" % lbl).text = str(size_quicklook)
            xml.find("product/quicklook/%s/dimension/rows" % lbl).text = str(self.nrows)
            xml.find("product/quicklook/%s/dimension/columns" % lbl).text = str(self.ncols)
            # xml.find("product/quicklook/%s/dimension/dimensionGeographic/longitude" % lbl).text = 'NA'
            # xml.find("product/quicklook/%s/dimension/dimensionGeographic/latitude" % lbl).text = 'NA'

        # productFileInformation
        ########################

        # get L1B product file information
        l1b_fileinfos = xmlSubtree2dict(xml, 'product/productFileInformation/')

        # clear old L1B file information in XML
        pFI_root = xml.findall('product/productFileInformation')[0]
        pFI_root.clear()

        # recreate sub-elements for productFileInformation according to L2A file information
        for i, fileInfo in enumerate(self.fileinfos):
            fn_l1b_exp = fileInfo['name'].replace('L2A', '*').replace('-SPECTRAL_IMAGE', '-SPECTRAL_IMAGE_VNIR')
            l1b_fileInfo = [fI for fI in l1b_fileinfos.values() if fnmatch.fnmatch(fI['name'], fn_l1b_exp)]

            if l1b_fileInfo:
                # TODO update file size of METADATA.XML (has not been written yet)
                fileInfo['size'] = fileInfo['size'] or l1b_fileInfo[0]['size']
                fileInfo['version'] = fileInfo['version'] or l1b_fileInfo[0]['version']
            else:
                # FIXME if no L1B equivalent is found for the file to be written, the file version will be empty ('')
                pass

            sub = ElementTree.SubElement(pFI_root, 'file', number=str(i))

            for k, kw in zip(['name', 'size', 'version', 'format'], [{}, {'unit': 'kbytes'}, {}, {}]):
                ele = ElementTree.SubElement(sub, k, **kw)
                ele.text = str(fileInfo[k])

        # TODO update product/ortho/projection
        #      {UTM_ZoneX_North, UTM_ZoneX_South (where X in {1..60}), Geographic, LAEA-ETRS89, NA}
        xml.find('product/ortho/resolution').text = str(enmap_xres)
        xml.find('product/ortho/resampling').text = self.cfg.ortho_resampAlg

        # band statistics
        #################

        if self.band_means is None or self.band_stds is None:
            raise ValueError('Band statistics have not yet been computed. Compute them first by calling '
                             'add_band_statistics()!')

        bs = "specific/bandStatistics/bandID/"
        for ele, mean in zip(xml.findall(bs + "meanReflectance"), self.band_means):
            ele.text = str(mean)
        for ele, std in zip(xml.findall(bs + "stdDeviation"), self.band_stds):
            ele.text = str(std)

        xml_string = ElementTree.tostring(xml, encoding='unicode', pretty_print=True)

        return xml_string


def xmlSubtree2dict(xml_root, path_subtree) -> OrderedDict:
    outDict = OrderedDict()
    allEle = xml_root.findall(path_subtree)

    for ele in allEle:
        eleKey = '%s_%s' % (ele.tag, ele.get('number'))
        outDict[eleKey] = dict()
        for subele in ele:
            outDict[eleKey][subele.tag] = subele.text

    return outDict
