# -*- coding: utf-8 -*-
"""EnPT metadata modules. All object and functions regarding EnMAP metadata are implemented here."""

from datetime import datetime
from xml.etree import ElementTree
import logging
import os
from typing import Union, List  # noqa: F401
import numpy as np
import spectral as sp
from py_tools_ds.geo.vector.topology import Polygon, get_footprint_polygon  # noqa: F401  # flake8 issue
from geoarray import GeoArray

from ..options.config import EnPTConfig
from .srf import SRF


# Define L1B_product_props
L1B_product_props = dict(
    xml_detector_label=dict(
        VNIR='VNIR',
        SWIR='SWIR'
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

        :param detector_name: Name of the detector (VNIR or SWIR)
        :param logger:
        """
        self.cfg = config
        self.detector_name = detector_name  # type: str
        self.detector_label = L1B_product_props['xml_detector_label'][detector_name]
        self.logger = logger or logging.getLogger()

        # These lines are used to load path information
        self.data_filename = None  # type: str # detector data filename
        self.dead_pixel_filename = None  # type: str # filename of the dead pixel file
        self.quicklook_filename = None  # type: str # filename of the quicklook file
        self.cloud_mask_filename = None  # type: str # filename of the cloud mask file

        self.fwhm = None  # type: np.ndarray  # Full width half maximmum for each EnMAP band
        self.wvl_center = None  # type: np.ndarray  # Center wavelengths for each EnMAP band
        self.srf = None  # type: SRF  # SRF object holding the spectral response functions for each EnMAP band
        self.solar_irrad = None  # type: dict  # solar irradiance in [W/m2/nm] for eac band
        self.nwvl = None  # type: int  # Number of wave bands
        self.nrows = None  # type: int  # number of rows
        self.ncols = None  # type: int  # number of columns
        self.smile_coef = None  # type: np.ndarray  # smile coefficients needed for smile computation
        self.nsmile_coef = None  # type: int  # number of smile coefficients
        self.smile = None  # type: np.ndarray  # smile for each EnMAP image column
        self.l_min = None  # type: np.ndarray
        self.l_max = None  # type: np.ndarray
        self.lat_UL_UR_LL_LR = None  # type:  List[float, float, float, float]  # latitude coords for UL, UR, LL, LR
        self.lon_UL_UR_LL_LR = None  # type:  List[float, float, float, float]  # longitude coords for UL, UR, LL, LR
        self.ll_mapPoly = None  # type: Polygon  # footprint polygon in longitude/latitude map coordinates
        self.lats = None  # type: np.ndarray  # 2D array of latitude coordinates according to given lon/lat sampling
        self.lons = None  # type: np.ndarray  # 2D array of longitude coordinates according to given lon/lat sampling
        self.unit = ''  # type: str  # radiometric unit of pixel values
        self.unitcode = ''  # type: str  # code of radiometric unit
        self.preview_bands = None
        self.snr = None  # type: np.ndarray  # Signal to noise ratio as computed from radiance data

    # On this new version of read_data, we don't need anymore nsmile_coef (will be read from xml file)
    def read_metadata(self, path_xml, lon_lat_smpl):
        """
        Read the meadata of a specific EnMAP detector in sensor geometry
        :param path_xml: file path of the metadata file
        :param lon_lat_smpl:  number if sampling in lon, lat fields
        :return: None
        """
        xml = ElementTree.parse(path_xml).getroot()
        lbl = self.detector_label + "Detector"
        self.logger.info("Reading metadata for %s detector..." % self.detector_name)

        # read data filenames
        self.data_filename = xml.findall("ProductComponent/%s/Data/Filename" % lbl)[0].text
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
        scene_corner_coordinates = xml.findall("ProductComponent/%s/Data/SceneInformation/SceneCornerCoordinates" % lbl)
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
        self.lats = self.interpolate_corners(*self.lat_UL_UR_LL_LR, *lon_lat_smpl)
        self.lons = self.interpolate_corners(*self.lon_UL_UR_LL_LR, *lon_lat_smpl)

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
        path_snr_model = os.path.join(dir_snr_models, "SNR_D1.hdr" if self.detector_name == 'VNIR' else "SNR_D2.hdr")
        rad_data = np.array(rad_data)

        assert self.unitcode == 'TOARad'
        self.logger.info("Computing SNR for %s using %s" % (self.detector_name, path_snr_model))

        if self.detector_name == 'VNIR':
            gA = sp.open_image(path_snr_model)
            coeffs_highgain = gA[0:3, :, :]
            coeffs_lowgain = gA[3:6, :, :]
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
            coeffs = sp.open_image(path_snr_model)[:, :, :]
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

    def calc_solar_irradiance_CWL_FWHM_per_band(self) -> dict:
        from ..io.reader import Solar_Irradiance_reader

        self.logger.debug('Calculating solar irradiance...')

        sol_irr = Solar_Irradiance_reader(path_solar_irr_model=self.cfg.path_solar_irr, wvl_min_nm=350, wvl_max_nm=2500)

        irr_bands = {}
        for band in self.srf.bands:
            if self.srf[band] is None:
                irr_bands[band] = None
            else:
                WVL_band = self.srf.srfs_wvl if self.srf.wvl_unit == 'nanometers' else self.srf.srfs_wvl * 1000
                RSP_band = self.srf.srfs_norm01[band]
                sol_irr_at_WVL = np.interp(WVL_band, sol_irr[:, 0], sol_irr[:, 1], left=0, right=0)

                irr_bands[band] = round(np.sum(sol_irr_at_WVL * RSP_band) / np.sum(RSP_band), 2)

        return irr_bands


class EnMAP_Metadata_L1B_SensorGeo(object):
    """EnMAP Metadata class holding the metadata of the complete EnMAP product in sensor geometry incl. VNIR and SWIR.

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
    """

    def __init__(self, path_metaxml, config: EnPTConfig, logger=None):
        """Get a metadata object instance for the given EnMAP L1B product in sensor geometry.

        :param path_metaxml:  file path of the EnMAP L1B metadata XML file
        :param logger:  instance of logging.logger or subclassed
        """
        self.cfg = config
        self.logger = logger or logging.getLogger()
        self._path_xml = path_metaxml

        # defaults - Common
        self.observation_datetime = None  # type: datetime  # Date and Time of image observation
        self.geom_view_zenith = None  # type: float  # viewing zenith angle
        self.geom_view_azimuth = None  # type: float  # viewing azimuth angle
        self.geom_sun_zenith = None  # type: float  # sun zenith angle
        self.geom_sun_azimuth = None  # type: float  # sun azimuth angle
        self.mu_sun = None  # type: float  # needed by SICOR for TOARad > TOARef conversion
        self.earthSunDist = None  # type: float  # earth-sun distance # TODO doc correct?
        self.vnir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of VNIR only
        self.swir = None  # type: EnMAP_Metadata_L1B_Detector_SensorGeo # metadata of SWIR only
        self.detector_attrNames = ['vnir', 'swir']

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

        # read the acquisition time
        self.observation_datetime = \
            datetime.strptime(xml.findall("GeneralInfo/ProductInfo/ProductStartTime")[0].text, '%Y-%m-%dT%H:%M:%S.%fZ')

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

    def read_metadata(self, lon_lat_smpl):
        """
        Read the metadata of the entire EnMAP L1B product in sensor geometry
        :param lon_lat_smpl:  number if sampling point in lon, lat fields
        :return: None
        """

        # first read common metadata
        self.read_common_meta(self._path_xml)

        # define and read the VNIR metadata
        self.vnir = EnMAP_Metadata_L1B_Detector_SensorGeo('VNIR', config=self.cfg, logger=self.logger)
        self.vnir.read_metadata(self._path_xml, lon_lat_smpl)

        # define and read the SWIR metadata
        self.swir = EnMAP_Metadata_L1B_Detector_SensorGeo('SWIR', config=self.cfg, logger=self.logger)
        self.swir.read_metadata(self._path_xml, lon_lat_smpl)
